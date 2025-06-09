# -*- coding: utf-8 -*-
"""
Created on Sun Apr 20 19:51:25 2025

@author: Owner
"""

import numpy as np
# import matplotlib.pyplot as plt   # still commented out
from scipy.special import gammaln, hyp2f1
from scipy.integrate import quad
from scipy.optimize import least_squares
from scipy import stats
import multiprocessing as mp          # <<< NEW
import time

# -------------------------------------------------------------------
# 1D E collapse distributions 
# -------------------------------------------------------------------
def standard_kappa(E, T_par, T_perp, kappa):
    ask = E/(kappa*T_par)
    bsk = E/(kappa*T_perp)
    result = np.exp(gammaln(kappa) - np.log(T_perp) - 0.5*np.log(T_par*kappa) - gammaln(kappa - 0.5)
                    - (1. + kappa)*np.log1p(bsk) )*hyp2f1(0.5, 1. + kappa, 1.5, (bsk - ask)/(bsk + 1.))
    return result

def product_kappa(E, T_par, T_perp, kappa_par, kappa_perp):
    apk   = E / (T_par * kappa_par)
    bpk   = E / (T_perp * kappa_perp)
    pref = np.exp(gammaln(kappa_par + 1.0) - np.log(T_perp)- 0.5*np.log(T_par*kappa_par) - gammaln(kappa_par + 0.5))
                  

    def integrand_pk_lv(t):
        return np.exp(-(1.0 + kappa_par)* np.log1p(apk*t*t) -(1.0 + kappa_perp)* np.log1p(bpk*(1.0 - t*t)))

    integral_val, err = quad(integrand_pk_lv, 0.0, 1.0, limit=100,
                             epsabs=0.0, epsrel=1.49e-08)
    per_err = 100.0*abs(err/integral_val)
    if per_err > 0.1:
        print("Warning: Product Kappa scipy quad Percent error=",
              per_err, " for T_par, T_perp =", T_par,T_perp, " and kappa_par, kappa_perp =", kappa_par, kappa_perp, flush=True)
    return pref*integral_val

def fried_egg(E, T_par, T_perp, kappa):
    afe = E / T_par
    bfe = E/(kappa*T_perp)
    pref = 1./(T_perp*np.sqrt(T_par))
    def integrand_fe_lv(t):
        return np.exp(- afe*t*t -(1.0 + kappa)*np.log1p(bfe*(1.0 - t*t)) )
    integral_val, err = quad(integrand_fe_lv, 0.0, 1.0, limit=100,
                             epsabs=0.0, epsrel=1.49e-08)
    per_err = 100.0*abs(err/integral_val)
    if per_err > 0.1:
        print("Warning: Fried-Egg scipy quad Percent error=",
              per_err, " for T_par, T_perp =", T_par, T_perp, " and kappa =", kappa, flush=True)
    return pref*integral_val

# ---------- (softplus, q_to_f, f_to_q, logistic utilities) ----------
def softplus(x):
    return np.logaddexp(0.0, x)

def q_to_f(q):
    q = np.asarray(q, dtype=np.float64)
    d = softplus(q)
    b = np.concatenate(([0.0], np.cumsum(d)))
    logits = -b
    exp_logits = np.exp(logits - logits.max())
    return exp_logits/exp_logits.sum()

def f_to_q(f, eps=1e-200):
    f = np.asarray(f, dtype=np.float64)
    d = np.log(np.clip(f[:-1]/f[1:], eps, None))
    return np.log(np.clip(np.exp(d) - 1.0, eps, None))

def logistic(x):
    return 1.0/(1.0 + np.exp(-x))

def logistic_inv(u):
    return np.log(u/(1.-u))

# ---------- (transform_params_to_T, five_maxwellians, etc.) ----------
def inverse_stick_breaking(T, T_min, T_max):
    log_T       = np.log(T)
    log_T_min   = np.log(T_min)
    log_T_max   = np.log(T_max)
    alpha       = (log_T - log_T_min)/(log_T_max - log_T_min)
    u           = np.zeros(5)
    u[0]        = alpha[0]
    for i in range(1, 5):
        u[i]    = (alpha[i] - alpha[i-1])/(1.0 - alpha[i-1])
    return logistic_inv(u)

def transform_params_to_T(x, T_min, T_max):
    p     = x[4:]
    u     = logistic(p)
    alpha = np.zeros(5)
    alpha[0] = u[0]
    for i in range(1, 5):
        alpha[i] = alpha[i-1] + u[i]*(1 - alpha[i-1])
    log_T = np.log(T_min) + alpha*(np.log(T_max) - np.log(T_min))
    return np.exp(log_T)

def five_maxwellians(E, params, T_min, T_max):
    q   = params[:4]
    fi   = q_to_f(q)
    Ti   = transform_params_to_T(params, T_min, T_max)
    E   = np.atleast_1d(E)
    pref = fi/(Ti**1.5) # f_{i}/(T_{i}^{3/2})
    log_terms = -E[:, None]/Ti[None, :]
    np.clip(log_terms, -700, None, out=log_terms)
    return (pref*np.exp(log_terms)).sum(axis=1)

LN10 = np.log(1.10)
def resid_vec(params, E, data, T_min, T_max):
    model = five_maxwellians(E, params, T_min, T_max)
    EPS   = 1e-200
    r = (np.log(np.maximum(model, EPS)) -
         np.log(np.maximum(data , EPS)))/LN10
    return r

def fit_five_maxwellians(E, data, T_input, T_min, T_max):
    lower = np.full(9, -10.0)
    upper = np.full(9,  10.0)
    T_maxxx = 50.*T_input if (50.*T_input) < 500. else 500.
    T_init  = np.logspace(np.log10(T_input), np.log10(T_maxxx), 5)
    p_init  = inverse_stick_breaking(T_init, T_min, T_max)
    f_init  = np.array([0.4, 0.25, 0.18, 0.12, 0.05])
    q_init  = f_to_q(f_init)
    x0      = np.concatenate([q_init, p_init])
    result  = least_squares(resid_vec, x0,
                            args=(E, data, T_min, T_max),
                            bounds=(lower, upper),
                            method='trf', jac='3-point',
                            x_scale='jac', loss='linear',
                            f_scale=1.0, ftol=1e-15,
                            xtol=1e-15, gtol=1e-15,
                            max_nfev=100000, verbose=0)
    return result.x, result.cost

def extract_f_T(params, T_min, T_max):
    q = params[:4]
    return q_to_f(q), transform_params_to_T(params, T_min, T_max)

def compute_percent_errors(E, data, params, T_min, T_max):
    model = five_maxwellians(E, params, T_min, T_max)
    eps = 1e-200
    return 100.0*np.abs(model - data)/np.maximum(data, eps)

# -------------------------------------------------------------------
# NEW helper executed in worker process (1 κ‑value task)
# -------------------------------------------------------------------
def _process_kappa(j, E_vals, T_val_par, T_val_perp, k1, k2, T_min_fit, T_max_fit):
    data_standard  = np.array([standard_kappa(E, T_val_par, T_val_perp, k1) for E in E_vals])
    #data_product   = np.array([product_kappa(E,  T_val_par, T_val_perp, k2) for E in E_vals])
    data_fried_egg = np.array([fried_egg(E,      T_val_par, T_val_perp, k2) for E in E_vals])

    # Standard κ
    p_std, _ = fit_five_maxwellians(E_vals, data_standard,
                                    max(T_val_par, T_val_perp), T_min_fit, T_max_fit)
    f_std, T_std = extract_f_T(p_std, T_min_fit, T_max_fit)
    perr_std = compute_percent_errors(E_vals, data_standard,
                                      p_std, T_min_fit, T_max_fit).max()

    # Product κ
    # p_pk , _ = fit_five_maxwellians(E_vals, data_product,
    #                                T_val_par, T_val_perp, T_min_fit, T_max_fit)
    #f_pk, T_pk = extract_f_T(p_pk, T_min_fit, T_max_fit)
    #perr_pk = compute_percent_errors(E_vals, data_product,
    #                                 p_pk, T_min_fit, T_max_fit).max()

    # Fried‑Egg
    p_fe , _ = fit_five_maxwellians(E_vals, data_fried_egg,
                                    max(T_val_par, T_val_perp), T_min_fit, T_max_fit)
    f_fe, T_fe = extract_f_T(p_fe, T_min_fit, T_max_fit)
    perr_fe = compute_percent_errors(E_vals, data_fried_egg,
                                     p_fe, T_min_fit, T_max_fit).max()

    # Return everything needed to fill the big arrays
    return (j, T_std, f_std, T_fe, f_fe,
            perr_std, perr_fe)
    #return (j, T_std, f_std, T_pk, f_pk, T_fe, f_fe,
    #        perr_std, perr_pk, perr_fe)

# -------------------------------------------------------------------
# MAIN SCRIPT
# -------------------------------------------------------------------
if __name__ == "__main__":
    start_time = time.time()

    # ---------------- grid definitions ----------------
    T_min_par, T_max_par = 0.99, 100.
    num_points_T_par = 5
    T_grid_par = np.logspace(np.log10(T_min_par), np.log10(T_max_par), num_points_T_par)
    
    T_min_perp, T_max_perp = 0.99, 100.
    num_points_T_perp = 5
    T_grid_perp = np.logspace(np.log10(T_min_perp), np.log10(T_max_perp), num_points_T_perp)

    kappa_min, kappa_max = 1.51, 300.
    num_points_kappa = 5
    kappa_grid  = np.logspace(np.log10(kappa_min),  np.log10(kappa_max),
                              num_points_kappa)
    kappa_min2  = 1.01
    kappa_grid2 = np.logspace(np.log10(kappa_min2), np.log10(kappa_max),
                              num_points_kappa)

    E_min, E_max, num_points_E = 0.01, 500.0, 150

    # ---------------- result arrays ----------------
    T_i_max4_sk = np.full((5, num_points_T_par, num_points_T_perp, num_points_kappa), np.nan)
    f_i_max4_sk = np.full_like(T_i_max4_sk, np.nan)
    #T_i_max4_pk = np.full_like(T_i_max4_sk, np.nan)
    #f_i_max4_pk = np.full_like(T_i_max4_sk, np.nan)
    T_i_max4_fe = np.full_like(T_i_max4_sk, np.nan)
    f_i_max4_fe = np.full_like(T_i_max4_sk, np.nan)

    max_rel_err_sk = np.full((num_points_T_par, num_points_T_perp, num_points_kappa), np.nan)
    # max_rel_err_pk = np.full_like(max_rel_err_sk, np.nan)
    max_rel_err_fe = np.full_like(max_rel_err_sk, np.nan)

    
    # before
    # pool = mp.Pool(processes=mp.cpu_count())   # grabs 24 logical threads

    # after (16 workers = 8 P‑cores + 8 E‑cores)
    pool = mp.Pool(processes=16)

    # --------------- main double loop ---------------
    for i, T_val_par in enumerate(T_grid_par):
        for k, T_val_perp in enumerate(T_grid_perp):
            print(" i =", i, f"of {num_points_T_par - 1}, k = ", k,f"of {num_points_T_perp - 1}", flush=True)
            E_vals_i = np.logspace(np.log10(E_min),
                                   np.log10(min(E_max,max( 100.*T_val_par, 100.*T_val_perp))),
                                   num_points_E)
    
            # prepare argument tuples for every kappa j
            arg_list = [(j, E_vals_i, T_val_par, T_val_perp,
                         kappa_grid[j], kappa_grid2[j],
                         0.01, 750.0)
                        for j in range(num_points_kappa)]
    
            # parallel execution over all j’s
            results = pool.starmap(_process_kappa, arg_list)
    
            # write results back into the big arrays
            #for (j, T_std, f_std, T_pk, f_pk, T_fe, f_fe,
            #     perr_std, perr_pk, perr_fe) in results:
            for (j, T_std, f_std, T_fe, f_fe,
                 perr_std, perr_fe) in results:
    
                T_i_max4_sk[:, i, k, j] = T_std
                f_i_max4_sk[:, i, k, j] = f_std
    
                #T_i_max4_pk[:, i, k, j] = T_pk
                #f_i_max4_pk[:, i, k, j] = f_pk
    
                T_i_max4_fe[:, i, k, j] = T_fe
                f_i_max4_fe[:, i, k, j] = f_fe
    
                max_rel_err_sk[i, k, j] = perr_std
                #max_rel_err_pk[i, k, j] = perr_pk
                max_rel_err_fe[i, k, j] = perr_fe

    # tidy up pool
    pool.close()
    pool.join()

    # ---------------- save ----------------
    np.savez(
        "lsquares_fixed_order_fi_and_Ti_local_solver_5x5x5_Tpar_vs_Tperp_vs_kappa_indiv_over_Elog150pnts_0.01-min_of_500_or_100xTpar_or_100xTperp_sk_and_fe_only.npz",
        Tpar=T_grid_par,
        Tperp=T_grid_perp,
        kappa=kappa_grid,
        kappa2=kappa_grid2,
        max_rel_err_sk=max_rel_err_sk,
        #max_rel_err_pk=max_rel_err_pk,
        max_rel_err_fe=max_rel_err_fe,
        T_i_max4_sk=T_i_max4_sk,
        #T_i_max4_pk=T_i_max4_pk,
        T_i_max4_fe=T_i_max4_fe,
        f_i_max4_sk=f_i_max4_sk,
        #f_i_max4_pk=f_i_max4_pk,
        f_i_max4_fe=f_i_max4_fe
    )

    print('nanmax of max rel err for  sk, pk, fe =',
          np.nanmax(max_rel_err_sk),
          #np.nanmax(max_rel_err_pk),
          np.nanmax(max_rel_err_fe), flush=True)
    print("--- %s seconds ---" % (time.time() - start_time), flush=True)
