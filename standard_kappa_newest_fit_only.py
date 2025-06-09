# -*- coding: utf-8 -*-
"""
Created on Mon Apr 21 07:47:12 2025

@author: Owner
"""

# -*- coding: utf-8 -*-
"""
5‑Maxwellian surrogate for a STANDARD ISOTROPIC κ‑distribution
 • physics‑based E_max(κ) from tail tolerance ε = 1e‑7
 • log‑ratio residual with 10 % penalty‑wall, √E base weight
 • adaptive re‑polish (≤8 passes, per‑bin weight boost)
 • %‑error checked only where data > SIG
 • multiprocessing over κ
"""

import numpy as np, time, multiprocessing as mp
from   scipy.integrate import quad
from   scipy.optimize  import least_squares
from   scipy.special   import gammaln

# ---- global numerical constants ---------------------------------------
EPS       = 1e-80     # floor for log arguments
SIG       = 1e-60     # significance threshold
LN10      = np.log(1.10)
EPS_TAIL  = 1e-7      # tail tolerance ⇒ emissivity error < 1e‑6

# ---- isotropic κ distribution -----------------------------------------
def kappa_iso(E, T, kappa):
    a = E/(kappa*T)
    log_pref = (gammaln(kappa) - gammaln(kappa-0.5)
                - 1.5*np.log(T) - 0.5*np.log(kappa))
    return np.exp(log_pref - (1.0+kappa)*np.log1p(a))

# ---- stick‑breaking helpers -------------------------------------------
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
    u = np.empty(5);  u[0] = alpha[0]
    for i in range(1,5):
        u[i]=(alpha[i]-alpha[i-1])/max(1-alpha[i-1],1e-12)
    return logistic_inv(u)

def p_to_T(p, T_min, T_max):
    u = logistic(p)
    alpha = np.empty(5); alpha[0]=u[0]
    for i in range(1,5):
        alpha[i] = alpha[i-1] + u[i]*(1-alpha[i-1])
    return np.exp(np.log(T_min)+alpha*(np.log(T_max)-np.log(T_min)))

# ---- five‑Maxwellian spectrum -----------------------------------------
def five_maxwellians(E, params, T_min, T_max):
    fi = q_to_f(params[:4])
    Ti = p_to_T(params[4:], T_min, T_max)
    E  = np.atleast_1d(E)
    return (fi/Ti**1.5 * np.exp(-E[:,None]/Ti[None,:])).sum(axis=1)

# ---- residual with penalty‑wall ---------------------------------------
def residual(params, E, data, w, T_min, T_max):
    model = five_maxwellians(E, params, T_min, T_max)
    raw   = (np.log(np.maximum(model,EPS))-np.log(np.maximum(data,EPS)))/LN10
    bad   = np.abs(raw)>1
    raw[bad]=np.sign(raw[bad])*(1+10*(np.abs(raw[bad])-1))
    return w*raw

# ---- %‑error on significant bins --------------------------------------
def percent_err(E, data, params, T_min, T_max):
    model = five_maxwellians(E, params, T_min, T_max)
    mask  = data > SIG
    err   = np.zeros_like(model)
    if np.any(mask):
        ratio     = np.maximum(model[mask], EPS)/data[mask]
        err[mask] = 100*np.abs(ratio-1)
    return err

# ---- single LSQ pass ---------------------------------------------------
def lsq_pass(E, data, w, x0, bounds, T_min, T_max, nfev=30000):
    return least_squares(residual, x0,
                         args=(E,data,w,T_min,T_max),
                         bounds=bounds, method='trf', jac='3-point',
                         x_scale='jac', ftol=1e-12, xtol=1e-12, gtol=1e-12,
                         max_nfev=nfev, verbose=0).x

# ---- worker for one κ --------------------------------------------------
def _proc_kappa(j, E, T, kappa, T_min_fit, T_max_fit):
    data = np.array([kappa_iso(e, T, kappa) for e in E])

    w     = E/E.max()                         # base √E weight ≡ E/E_max
    bnds  = (np.r_[-10*np.ones(4), -15*np.ones(5)],
             np.r_[ 10*np.ones(4),  15*np.ones(5)])

    # initial parameters
    T_init = np.logspace(np.log10(T), np.log10(min(500,50*T)), 5)
    p0 = inv_stick_T(T_init, T_min_fit, T_max_fit)
    q0 = f_to_q([0.4,0.25,0.18,0.12,0.05])
    x  = np.concatenate([q0,p0])

    # adaptive re‑polish (≤8 passes)
    for _ in range(8):
        x   = lsq_pass(E,data,w,x,bnds,T_min_fit,T_max_fit,20000)
        err = percent_err(E,data,x,T_min_fit,T_max_fit)
        if err.max()<=10: break
        w *= 1 + np.minimum(10, err/10)

    return j, p_to_T(x[4:],T_min_fit,T_max_fit), q_to_f(x[:4]), err.max()

# =======================================================================
if __name__=="__main__":
    t0 = time.time()

    # -------- grids (T × κ) -------------------------------------------
    # T_grid     = np.logspace(np.log10(0.97), np.log10(69.0), 5)   # eV
    # 14 elements
    T_grid = np.array([0.97,1.,2.5,5.,7.5,10.,15.,20.,30.,40.,50.,60.,65.,69.])
    
    # kappa_grid = np.logspace(np.log10(1.51), np.log10(300.), 5)
    # 24  elements
    kappa_grid  = np.array([1.51,1.6,1.7,1.8,1.9,2.,2.5,3.,4.,5.,6.,7.,10.,15.,20.,25.,30.,40.,50.,75.,100.,150.,200.,300.])
    E_min_global, nE = 0.01, 150

    T_i  = np.full((5,len(T_grid),len(kappa_grid)), np.nan)
    f_i  = np.full_like(T_i,np.nan)
    mErr = np.full((len(T_grid),len(kappa_grid)), np.nan)

    pool = mp.Pool(mp.cpu_count())

    for iT, T in enumerate(T_grid):
        print(f"T[{iT}/9] = {T:.3f} eV", flush=True)

        # kappa loop ----------------------------------------------------
        tasks=[]
        for j,kappa in enumerate(kappa_grid):

            # physics‑based tail cutoff  E_max = κT((ε)^(-1/(1+κ))−1)
            tail_Emax = kappa*T*((EPS_TAIL)**(-1.0/(1.0+kappa))-1.0)
            E_max_dyn = min(500.0, 100.0*T, tail_Emax)
            E_max_dyn = max(E_max_dyn, 2.0*T)          # grid non‑empty

            E_vals=np.logspace(np.log10(E_min_global),
                               np.log10(E_max_dyn), nE)

            tasks.append((j,E_vals,T,kappa,0.01,750.0))

        for j,T_opt,f_opt,err in pool.starmap(_proc_kappa,tasks):
            T_i[:,iT,j] = T_opt
            f_i[:,iT,j] = f_opt
            mErr[iT,j]  = err
            if err>10:
                print(f"  Warning: {err:.1f}%  κ={kappa_grid[j]:.3g}", flush=True)

    pool.close(); pool.join()

    np.savez_compressed("iso_std_kappa_to_5Mxw_phys_cutoff_Txkappa_14x24.npz",
                        T=T_grid, kappa=kappa_grid,
                        T_i=T_i, f_i=f_i, max_err=mErr)

    print("nan‑max significant‑bin error =", np.nanmax(mErr))
    print(f"--- {time.time()-t0:.1f} s ---", flush=True)
