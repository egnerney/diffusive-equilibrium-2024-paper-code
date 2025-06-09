# -*- coding: utf-8 -*-
"""
Parallel 5‑Maxwellian fit for 4‑D product‑κ grid (flattened task list)

Grid sizes      :  A(10) × T⊥(14) × λ(9) × κ⊥(17)  = 21420 cells
Energy grid     : 150 log points, 0.01 eV … Emax(κ)
Parallelism     : one task per cell, distributed to all CPU cores
Output          : T_i[5,…], f_i[5,…]  in same shape as original script
"""

import numpy as np, time, multiprocessing as mp
from   scipy.integrate import quad
from   scipy.optimize  import least_squares
from   scipy.special   import gammaln

# ----- constants --------------------------------------------------------
EPS, SIG = 1e-80, 1e-60
LN10     = np.log(1.10)
EPS_TAIL = 1e-7

# ----- product‑κ definition --------------------------------------------
def kappa_product(E, Tpar, Tperp, kpar, kperp):
    a,b = E/(kpar*Tpar), E/(kperp*Tperp)
    log_pref = (gammaln(kpar+1.) - gammaln(kpar+0.5)
                - np.log(Tperp) - 0.5*np.log(Tpar*kpar))
    def integrand(t):
        return np.exp(-(1+kpar )*np.log1p(a*t*t)
                      -(1+kperp)*np.log1p(b*(1-t*t)))
    I, err = quad(integrand, 0.0, 1.0, limit=100,
                  epsabs=0.0, epsrel=1.49e-8)
    return np.exp(log_pref)*I

# ----- stick‑breaking helpers ------------------------------------------
def softplus(x): return np.logaddexp(0.0, x)
def q_to_f(q):
    d = softplus(q); b = np.concatenate(([0.0], np.cumsum(d)))
    e = np.exp(-b - (-b).max()); return e/e.sum()
def f_to_q(f, eps=1e-12):
    f = np.asarray(f,float)
    d = np.log(np.clip(f[:-1]/f[1:], eps, None))
    return np.log(np.clip(np.exp(d)-1.0, eps, None))
def logistic(x): return 1/(1+np.exp(-x))
def logistic_inv(u): return np.log(u/(1-u))
def inv_stick_T(T,Tmin,Tmax):
    a=(np.log(T)-np.log(Tmin))/(np.log(Tmax)-np.log(Tmin))
    u=np.empty(5);u[0]=a[0]
    for i in range(1,5): u[i]=(a[i]-a[i-1])/max(1-a[i-1],1e-12)
    return logistic_inv(u)
def p_to_T(p,Tmin,Tmax):
    u=logistic(p);a=np.empty(5);a[0]=u[0]
    for i in range(1,5): a[i]=a[i-1]+u[i]*(1-a[i-1])
    return np.exp(np.log(Tmin)+a*(np.log(Tmax)-np.log(Tmin)))

# ----- 5‑Maxwellian spectrum -------------------------------------------
def five_maxwellians(E,params,Tmin,Tmax):
    fi=q_to_f(params[:4]);Ti=p_to_T(params[4:],Tmin,Tmax)
    return (fi/Ti**1.5*np.exp(-E[:,None]/Ti[None,:])).sum(axis=1)

# ----- residual & helpers ----------------------------------------------
def residual(params,E,data,w,Tmin,Tmax):
    model=five_maxwellians(E,params,Tmin,Tmax)
    raw=(np.log(np.maximum(model,EPS))-np.log(np.maximum(data,EPS)))/LN10
    bad=np.abs(raw)>1; raw[bad]=np.sign(raw[bad])*(1+10*(np.abs(raw[bad])-1))
    return w*raw

def percent_err(E,data,params,Tmin,Tmax):
    model=five_maxwellians(E,params,Tmin,Tmax)
    mask=data>SIG; err=np.zeros_like(model)
    if np.any(mask):
        ratio=np.maximum(model[mask],EPS)/data[mask]
        err[mask]=100*np.abs(ratio-1)
    return err

def lsq_pass(E,data,w,x0,bounds,Tmin,Tmax,nfev=30000):
    return least_squares(residual,x0,
        args=(E,data,w,Tmin,Tmax),
        bounds=bounds,method='trf',jac='3-point',x_scale='jac',
        ftol=1e-12,xtol=1e-12,gtol=1e-12,max_nfev=nfev,verbose=0).x

# ----- ONE cell task ----------------------------------------------------
def _cell_task(idxA,idxT,idxL,idxK, A,Tperp,lmbd,kperp,
               nE,Emin,Tmin_fit,Tmax_fit):
    Tpar=Tperp/A; kpar=kperp/lmbd
    # physics‑based E_max -------------------------------------------------
    kmin=min(kperp,kpar)
    tail_Emax=kmin*max(Tperp,Tpar)*(EPS_TAIL**(-1/(1+kmin))-1)
    E_max = min(500.,100*max(Tperp,Tpar),tail_Emax)
    E_max = max(E_max,2*max(Tperp,Tpar))
    E = np.logspace(np.log10(Emin), np.log10(E_max), nE)

    data = np.array([kappa_product(e,Tpar,Tperp,kpar,kperp) for e in E])
    w    = E/E.max()
    bounds=(np.r_[-10*np.ones(4),-15*np.ones(5)],
            np.r_[ 10*np.ones(4), 15*np.ones(5)])

    # initial guess
    T_init=np.logspace(np.log10(Tperp), np.log10(min(500,50*Tperp)),5)
    p0=inv_stick_T(T_init,Tmin_fit,Tmax_fit)
    q0=f_to_q([0.4,0.25,0.18,0.12,0.05])
    x =np.concatenate([q0,p0])

    # adaptive re‑polish
    for _ in range(8):
        x = lsq_pass(E,data,w,x,bounds,Tmin_fit,Tmax_fit,20000)
        err = percent_err(E,data,x,Tmin_fit,Tmax_fit)
        if err.max()<=10: break
        w *= 1+np.minimum(10,err/10)
    
    if (idxT==0) and (idxL==0) and (idxK == 0):
        print(f"A[{idxA}/10]={A:.3f}",flush=True)
        
    return (idxA,idxT,idxL,idxK,
            p_to_T(x[4:],Tmin_fit,Tmax_fit),
            q_to_f(x[:4]), err.max())

# =======================================================================
if __name__=="__main__":
    t0=time.time()

    # -------- explicit grids (10×14×9×17) -------------------------------
    A_grid      = np.array([0.28,0.5,0.7,0.9,1.0,1.1,1.3,1.5,1.7,1.9])
    Tperp_grid  = np.array([0.97,1.0,2.5,5.,7.5,10.,15.,20.,30.,50.,
                            75.,100.,125.,144.])
    lambda_grid = np.array([1.,1.15,1.35,1.55,1.75,1.95,2.15,2.35,2.6])
    kperp_grid  = np.array([1.01,1.1,1.3,1.7,2.,2.5,3.,4.,5.,7.,10.,15.,
                            20.,30.,50.,150.,300.])

    NA,NT,NL,NK = len(A_grid),len(Tperp_grid),len(lambda_grid),len(kperp_grid)
    Emin,nE = 0.01,150

    T_i = np.full((5,NA,NT,NL,NK),np.nan)
    f_i = np.full_like(T_i,np.nan)
    mErr= np.full((NA,NT,NL,NK),np.nan)

    print(f"Total cells = {NA*NT*NL*NK}  —  using {mp.cpu_count()} cores")

    # -------- build flat task list --------------------------------------
    tasks=[]
    for ia,A in enumerate(A_grid):
        for it,Tp in enumerate(Tperp_grid):
            for il,lmb in enumerate(lambda_grid):
                for ik,kp in enumerate(kperp_grid):
                    tasks.append((ia,it,il,ik, A,Tp,lmb,kp,
                                  nE,Emin,0.01,750.0))

    # -------- run pool ---------------------------------------------------
    with mp.Pool() as pool:
        for (ia,it,il,ik,Topt,fopt,err) in pool.starmap(_cell_task,tasks):
            T_i[:,ia,it,il,ik]=Topt
            f_i[:,ia,it,il,ik]=fopt
            mErr[ia,it,il,ik] = err
            if err>10:
                print(f"  Warn {err:.1f}%  A={A_grid[ia]:.2f} "
                      f"T⊥={Tperp_grid[it]:.1f} λ={lambda_grid[il]:.2f} "
                      f"κ⊥={kperp_grid[ik]:.2g}", flush=True)

    np.savez_compressed("product_kappa_to_5Mxw_phys_cutoff_A10_T14_L9_K17.npz",
                        A=A_grid,Tperp=Tperp_grid,lambda_=lambda_grid,
                        kappa_perp=kperp_grid,
                        T_i=T_i,f_i=f_i,max_err=mErr)

    print("nan‑max significant‑bin error =",np.nanmax(mErr))
    print(f"--- {time.time()-t0:.1f} s ---",flush=True)
