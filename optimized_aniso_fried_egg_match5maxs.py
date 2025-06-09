# -*- coding: utf-8 -*-
"""
5‑Maxwellian surrogate for an anisotropic Fried‑Egg κ‑distribution
• log‑ratio residual with 10 % penalty‑wall
• √E base‑weight, adaptive re‑polishing (≤ 8 passes)
• %‑error evaluated only on significant bins (data > SIG)
• physics‑based κ‑dependent energy cut‑off:
      E_max(kappa) = kappa*T_perp*((ε)^(-1/(1+κ)) − 1)
      final E_max  = min(500 eV, 100·max(T_perp,T_par), E_max(kappa))
"""

import numpy as np, multiprocessing as mp, time
from   scipy.integrate import quad
from   scipy.optimize  import least_squares

# ---------- global numerical constants ---------------------------------
EPS       = 1e-80     # floor for log arguments
SIG       = 1e-60     # bins below this are ignored in %‑error
LN10      = np.log(1.10)
EPS_TAIL  = 1e-7      # target fractional tail contribution (<10⁻⁶ emissivity)

# ---------- Fried‑Egg spectrum -----------------------------------------
def fried_egg(E, Tpar, Tperp, kappa):
    a, b = E/Tpar, E/(kappa*Tperp)
    pref = 1.0/(Tperp*np.sqrt(Tpar))
    func = lambda t: np.exp(-a*t*t - (1+kappa)*np.log1p(b*(1-t*t)))
    I, _ = quad(func, 0.0, 1.0, limit=100, epsabs=0.0, epsrel=1.49e-8)
    return pref*I

# ---------- stick‑breaking helpers -------------------------------------
def softplus(x): return np.logaddexp(0.0, x)

def q_to_f(q):
    d = softplus(q)
    b = np.concatenate(([0.0], np.cumsum(d)))
    e = np.exp(-b - (-b).max())
    return e/e.sum()

def f_to_q(f, eps=1e-12):
    f = np.asarray(f, float)
    d = np.log(np.clip(f[:-1]/f[1:], eps, None))
    return np.log(np.clip(np.exp(d)-1.0, eps, None))

def logistic(x):      return 1/(1+np.exp(-x))
def logistic_inv(u):  return np.log(u/(1-u))

def inv_stick_T(T, T_min, T_max):
    alpha = (np.log(T)-np.log(T_min))/(np.log(T_max)-np.log(T_min))
    u = np.empty(5); u[0]=alpha[0]
    for i in range(1,5):
        u[i]=(alpha[i]-alpha[i-1])/max(1-alpha[i-1],1e-12)
    return logistic_inv(u)

def p_to_T(p, T_min, T_max):
    u  = logistic(p)
    alpha = np.empty(5); alpha[0]=u[0]
    for i in range(1,5):
        alpha[i]=alpha[i-1]+u[i]*(1-alpha[i-1])
    return np.exp(np.log(T_min)+alpha*(np.log(T_max)-np.log(T_min)))

# ---------- 5‑Maxwellian spectrum --------------------------------------
def five_maxwellians(E, params, T_min, T_max):
    fi = q_to_f(params[:4])
    Ti = p_to_T(params[4:], T_min, T_max)
    E  = np.atleast_1d(E)
    return (fi/Ti**1.5 * np.exp(-E[:,None]/Ti[None,:])).sum(axis=1)

# ---------- residual with penalty‑wall ---------------------------------
def residual(params, E, data, w, T_min, T_max):
    model = five_maxwellians(E, params, T_min, T_max)
    raw   = (np.log(np.maximum(model,EPS))-np.log(np.maximum(data,EPS)))/LN10
    bad   = np.abs(raw) > 1
    raw[bad] = np.sign(raw[bad])*(1 + 10*(np.abs(raw[bad])-1))
    return w*raw

# ---------- %‑error on significant bins --------------------------------
def percent_err(E, data, params, T_min, T_max):
    model = five_maxwellians(E, params, T_min, T_max)
    mask  = data > SIG
    err   = np.zeros_like(model)
    if np.any(mask):
        ratio     = np.maximum(model[mask], EPS) / data[mask]
        err[mask] = 100*np.abs(ratio-1)
    return err

# ---------- LSQ helper --------------------------------------------------
def lsq_pass(E, data, w, x0, bounds, T_min, T_max, nfev=30000):
    return least_squares(residual, x0,
                         args=(E, data, w, T_min, T_max),
                         bounds=bounds, method='trf', jac='3-point',
                         x_scale='jac', ftol=1e-12, xtol=1e-12, gtol=1e-12,
                         max_nfev=nfev, verbose=0).x

# ---------- worker for one κ -------------------------------------------
def _proc_kappa(j, E, Tpar, Tperp, kappa, T_min_fit, T_max_fit):
    data = np.array([fried_egg(e, Tpar, Tperp, kappa) for e in E])

    w    = E/E.max()                                   # base weight
    bounds = (np.r_[-10*np.ones(4), -15*np.ones(5)],
              np.r_[ 10*np.ones(4),  15*np.ones(5)])

    # initial parameters
    T_ref  = max(Tpar, Tperp)
    T_init = np.logspace(np.log10(T_ref),
                         np.log10(min(500,50*T_ref)),5)
    p0 = inv_stick_T(T_init, T_min_fit, T_max_fit)
    q0 = f_to_q([0.4,0.25,0.18,0.12,0.05])
    x  = np.concatenate([q0,p0])

    # adaptive re‑polish
    for _ in range(8):
        x   = lsq_pass(E, data, w, x, bounds, T_min_fit, T_max_fit, 20000)
        err = percent_err(E, data, x, T_min_fit, T_max_fit)
        if err.max() <= 10: break
        w *= 1 + np.minimum(10, err/10)                # boost bad bins

    return j, p_to_T(x[4:], T_min_fit, T_max_fit), q_to_f(x[:4]), err.max()

# ================================================================= MAIN
if __name__ == "__main__":
    t0 = time.time()

    # ---------------- grids --------------------
    #A_grid      = np.linspace(0.84,1.12,5)
    #Tperp_grid  = np.linspace(0.999,12.47,5)
    #kappa_grid  = np.logspace(np.log10(1.23),np.log10(300),5)
    
    A_grid = np.array([0.84,0.9,0.95,1.0,1.05,1.1,1.12]) # 7 elements
    Tperp_grid  = np.array([0.99,1.0,2.,2.5,3.,4.,5.0,6.,7.5,10.,11.,12.,12.5]) #13 elements
    # 23 elements
    kappa_grid = np.array([1.23,1.3,1.5,1.7,1.8,1.9,2.,2.5,3.,4.,5.,7.,10.,15.,20.,25.,30.,50.,70., 100., 150.,200., 300.])


    E_min_global, nE = 0.01, 150

    T_i   = np.full((5,len(A_grid),len(Tperp_grid),len(kappa_grid)),np.nan)
    f_i   = np.full_like(T_i, np.nan)
    m_err = np.full((len(A_grid),len(Tperp_grid),len(kappa_grid)),np.nan)

    pool = mp.Pool(mp.cpu_count())

    for ia,A in enumerate(A_grid):
        for ip,Tperp in enumerate(Tperp_grid):
            Tpar = Tperp/A
            print(f"A[{ia}/6]={A:.3f}, Tperp[{ip}/12]={Tperp:.3f}", flush=True)

            tasks=[]
            for j,kappa in enumerate(kappa_grid):

                # physics‑based tail cutoff ---------------------------
                tail_Emax = kappa*Tperp*((EPS_TAIL)**(-1.0/(1.0+kappa))-1.0)

                # final cap (lab limit and legacy 100*T) --------------
                E_max_dyn = min(500.0,
                                100.0*max(Tperp,Tpar),
                                tail_Emax)

                # ensure at least 2*T so grid non‑empty ---------------
                E_max_dyn = max(E_max_dyn, 2.0*max(Tperp,Tpar))

                E_vals = np.logspace(np.log10(E_min_global),
                                     np.log10(E_max_dyn), nE)

                tasks.append((j,E_vals,Tpar,Tperp,kappa,0.01,750.0))

            for j,T_opt,f_opt,err in pool.starmap(_proc_kappa, tasks):
                T_i[:,ia,ip,j] = T_opt
                f_i[:,ia,ip,j] = f_opt
                m_err[ia,ip,j] = err
                if err > 10:
                    print(f"  Warning: {err:.1f}% κ={kappa_grid[j]:.3g}", flush=True)

    pool.close(); pool.join()

    np.savez_compressed("fried_egg_5Mxw_fit_phys_cutoff_7x13x23_AxTperpxkappa.npz",
                        A=A_grid, Tperp=Tperp_grid, kappa=kappa_grid,
                        T_i=T_i, f_i=f_i, max_err=m_err)

    print("nan‑max significant‑bin error =", np.nanmax(m_err))
    print(f"--- {time.time()-t0:.1f} s ---", flush=True)
