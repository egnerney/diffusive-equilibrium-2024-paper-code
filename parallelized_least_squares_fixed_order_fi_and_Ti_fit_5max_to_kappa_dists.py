# -*- coding: utf-8 -*-
"""
Created on Sun Apr 20 05:56:56 2025

@author: Owner
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Apr 17 23:35:04 2025
@author: Owner
Parallelised inner‑kappa loop with multiprocessing.Pool (MEDIUM impact)
"""
import numpy as np
# import matplotlib.pyplot as plt   # still commented out
from scipy.special import gammaln
from scipy.integrate import quad
from scipy.optimize import least_squares
from scipy import stats
import multiprocessing as mp          # <<< NEW
import time

# -------------------------------------------------------------------
# (all the distribution, utility, fitting routines are UNCHANGED)
# -------------------------------------------------------------------
def standard_kappa(E, T, kappa):
    result = np.exp(gammaln(kappa) - gammaln(kappa - 0.5)
                    - 0.5*np.log(kappa) - (1.0 + kappa)*np.log1p(E /(kappa*T)))
    return result

def product_kappa(E, T, kappa):
    a2   = E / (T * kappa)
    pref = np.exp(gammaln(kappa + 1.0) - gammaln(kappa + 0.5)
                  - 0.5*np.log(kappa))

    def integrand_pk_lv(t):
        return np.exp(-(1.0 + kappa)
                      * (np.log1p(a2*(1.0 - t*t)) + np.log1p(a2*t*t)))

    integral_val, err = quad(integrand_pk_lv, 0.0, 1.0, limit=100,
                             epsabs=0.0, epsrel=1.49e-08)
    per_err = 100.0*abs(err/integral_val)
    if per_err > 0.1:
        print("Warning: Product Kappa scipy quad Percent error=",
              per_err, " for T =", T, " and kappa =", kappa, flush=True)
    return pref*integral_val

def fried_egg(E, T, kappa):
    a = E / T
    def integrand_fe_lv(t):
        return np.exp(-(1.0 + kappa)*np.log1p((a/kappa)*(1.0 - t*t)) - a*t*t)
    integral_val, err = quad(integrand_fe_lv, 0.0, 1.0, limit=100,
                             epsabs=0.0, epsrel=1.49e-08)
    per_err = 100.0*abs(err/integral_val)
    if per_err > 0.1:
        print("Warning: Fried-Egg scipy quad Percent error=",
              per_err, " for T =", T, " and kappa =", kappa, flush=True)
    return integral_val

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
    f   = q_to_f(q)
    T   = transform_params_to_T(params, T_min, T_max)
    E   = np.atleast_1d(E)
    log_terms = -E[:, None]/T[None, :]
    np.clip(log_terms, -700, None, out=log_terms)
    return (f*np.exp(log_terms)).sum(axis=1)

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
def _process_kappa(j, E_vals, T_val, k1, k2, T_min_fit, T_max_fit):
    data_standard  = np.array([standard_kappa(E, T_val, k1) for E in E_vals])
    data_product   = np.array([product_kappa(E,  T_val, k2) for E in E_vals])
    data_fried_egg = np.array([fried_egg(E,      T_val, k2) for E in E_vals])

    # Standard κ
    p_std, _ = fit_five_maxwellians(E_vals, data_standard,
                                    T_val, T_min_fit, T_max_fit)
    f_std, T_std = extract_f_T(p_std, T_min_fit, T_max_fit)
    perr_std = compute_percent_errors(E_vals, data_standard,
                                      p_std, T_min_fit, T_max_fit).max()

    # Product κ
    p_pk , _ = fit_five_maxwellians(E_vals, data_product,
                                    T_val, T_min_fit, T_max_fit)
    f_pk, T_pk = extract_f_T(p_pk, T_min_fit, T_max_fit)
    perr_pk = compute_percent_errors(E_vals, data_product,
                                     p_pk, T_min_fit, T_max_fit).max()

    # Fried‑Egg
    p_fe , _ = fit_five_maxwellians(E_vals, data_fried_egg,
                                    T_val, T_min_fit, T_max_fit)
    f_fe, T_fe = extract_f_T(p_fe, T_min_fit, T_max_fit)
    perr_fe = compute_percent_errors(E_vals, data_fried_egg,
                                     p_fe, T_min_fit, T_max_fit).max()

    # Return everything needed to fill the big arrays
    return (j, T_std, f_std, T_pk, f_pk, T_fe, f_fe,
            perr_std, perr_pk, perr_fe)

# -------------------------------------------------------------------
# MAIN SCRIPT
# -------------------------------------------------------------------
if __name__ == "__main__":
    start_time = time.time()

    # ---------------- grid definitions ----------------
    T_min, T_max = 0.99, 100.
    num_points_T = 500
    T_grid = np.logspace(np.log10(T_min), np.log10(T_max), num_points_T)

    kappa_min, kappa_max = 1.51, 300.
    num_points_kappa = 500
    kappa_grid  = np.logspace(np.log10(kappa_min),  np.log10(kappa_max),
                              num_points_kappa)
    kappa_min2  = 1.01
    kappa_grid2 = np.logspace(np.log10(kappa_min2), np.log10(kappa_max),
                              num_points_kappa)

    E_min, E_max, num_points_E = 0.01, 500.0, 150

    # ---------------- result arrays ----------------
    T_i_max4_sk = np.full((5, num_points_T, num_points_kappa), np.nan)
    f_i_max4_sk = np.full_like(T_i_max4_sk, np.nan)
    T_i_max4_pk = np.full_like(T_i_max4_sk, np.nan)
    f_i_max4_pk = np.full_like(T_i_max4_sk, np.nan)
    T_i_max4_fe = np.full_like(T_i_max4_sk, np.nan)
    f_i_max4_fe = np.full_like(T_i_max4_sk, np.nan)

    max_rel_err_sk = np.full((num_points_T, num_points_kappa), np.nan)
    max_rel_err_pk = np.full_like(max_rel_err_sk, np.nan)
    max_rel_err_fe = np.full_like(max_rel_err_sk, np.nan)

    # --------------- multiprocessing pool -----------
    pool = mp.Pool(processes=mp.cpu_count())

    # --------------- main double loop ---------------
    for i, T_val in enumerate(T_grid):
        print(" i =", i, f"of {num_points_T - 1}", flush=True)
        E_vals_i = np.logspace(np.log10(E_min),
                               np.log10(min(E_max, 100.*T_val)),
                               num_points_E)

        # prepare argument tuples for every kappa j
        arg_list = [(j, E_vals_i, T_val,
                     kappa_grid[j], kappa_grid2[j],
                     0.01, 750.0)
                    for j in range(num_points_kappa)]

        # parallel execution over all j’s
        results = pool.starmap(_process_kappa, arg_list)

        # write results back into the big arrays
        for (j, T_std, f_std, T_pk, f_pk, T_fe, f_fe,
             perr_std, perr_pk, perr_fe) in results:

            T_i_max4_sk[:, i, j] = T_std
            f_i_max4_sk[:, i, j] = f_std

            T_i_max4_pk[:, i, j] = T_pk
            f_i_max4_pk[:, i, j] = f_pk

            T_i_max4_fe[:, i, j] = T_fe
            f_i_max4_fe[:, i, j] = f_fe

            max_rel_err_sk[i, j] = perr_std
            max_rel_err_pk[i, j] = perr_pk
            max_rel_err_fe[i, j] = perr_fe

    # tidy up pool
    pool.close()
    pool.join()

    # ---------------- save ----------------
    np.savez(
        "lsquares_fixed_order_fi_and_Ti_local_solver_500x500_T_vs_kappa_indiv_over_Elog150pnts_0.01-min_of_500_or_100xT.npz",
        T=T_grid,
        kappa=kappa_grid,
        kappa2=kappa_grid2,
        E_vals=E_vals_i,
        max_rel_err_sk=max_rel_err_sk,
        max_rel_err_pk=max_rel_err_pk,
        max_rel_err_fe=max_rel_err_fe,
        T_i_max4_sk=T_i_max4_sk,
        T_i_max4_pk=T_i_max4_pk,
        T_i_max4_fe=T_i_max4_fe,
        f_i_max4_sk=f_i_max4_sk,
        f_i_max4_pk=f_i_max4_pk,
        f_i_max4_fe=f_i_max4_fe
    )

    print('nanmax of max rel err for  sk, pk, fe =',
          np.nanmax(max_rel_err_sk),
          np.nanmax(max_rel_err_pk),
          np.nanmax(max_rel_err_fe), flush=True)
    print("--- %s seconds ---" % (time.time() - start_time), flush=True)
