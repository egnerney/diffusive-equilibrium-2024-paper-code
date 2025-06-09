# -*- coding: utf-8 -*-
"""
Created on Thu Apr 17 23:35:04 2025

@author: Owner
"""




import numpy as np
#import matplotlib.pyplot as plt

from scipy.special import  gammaln # gamma
from scipy.integrate import quad
#from scipy.optimize import minimize 
from scipy.optimize import least_squares
#from scipy.optimize import differential_evolution
from scipy import stats


import time






# ------------------------------
# 1. Define the "true" distributions
# ------------------------------
def standard_kappa(E, T, kappa):
    """
    Standard kappa distribution:
      f(E) = (Gamma(kappa) / [sqrt(kappa) * Gamma(kappa - 1/2)])
             * (1 + E/(kappa*T))^(-(1 + kappa))
    """
    
    result = np.exp(gammaln(kappa) - gammaln(kappa - 0.5) - 0.5*np.log(kappa) - (1.0 + kappa)*np.log1p(E /(kappa*T)) )
    return result

def product_kappa(E, T, kappa):
    """
    Product kappa distribution:
      f(E) = [Gamma(kappa+1)/( sqrt(kappa)*Gamma(kappa+0.5))]
             * Integral_{0 to 1} (1 + (E/(kappa*T))*(1 - t^2))^(-(1 + kappa))
                                 * (1 + (E/(kappa*T))*t^2)^(-(1 + kappa)) dt
    """

    a2 = E / (T * kappa)
    pref = np.exp(gammaln(kappa + 1.0) - gammaln(kappa + 0.5) - 0.5*np.log(kappa))
    def integrand_pk_lv(t):
        return np.exp( - (1.0 + kappa)*(np.log1p(a2*(1.0 - t*t)) + np.log1p( a2*t*t)))
                       

    integral_val, err = quad(integrand_pk_lv, 0.0, 1.0, limit=100, epsabs=0.0, epsrel=1.49e-08)
    per_err = 100.0 * abs(err / integral_val)
    if per_err > 0.1:
        print("Warning: Product Kappa scipy quad Percent error=", per_err,
              " for T =", T, " and kappa =", kappa, flush=True)
               
    return pref*integral_val

def fried_egg(E, T, kappa):
    """
    Fried-Egg distribution:
      f(E) = Integral_{0 to 1} [
                (1 + (E/(kappa*T))*(1 - t**2))^(-(1 + kappa)) * exp(-E*t^2 / T)
             ] dt
    """

    a = E / T

    def integrand_fe_lv(t):
        return np.exp(-(1.0 + kappa)*np.log1p((a/kappa)*(1.0 - t*t)) -a*t*t)
               
    integral_val, err = quad(integrand_fe_lv, 0.0, 1.0, limit=100, epsabs=0.0, epsrel=1.49e-08)
    per_err = 100.0 * abs(err / integral_val)
    if per_err > 0.1:
        print("Warning: Fried-Egg scipy quad Percent error=", per_err,
              " for T =", T, " and kappa =", kappa, flush=True)
     
    return integral_val


# ---- utility --------------------------------------------------------------
def softplus(x):
    # numerically stable: log(1+e^x)
    return np.logaddexp(0.0, x)

# ---- main forward map -----------------------------------------------------
def q_to_f(q):
    """
    Parameters
    ----------
    q : array_like, shape (4,)
        Unconstrained optimisation parameters.

    Returns
    -------
    f : ndarray, shape (5,)
        Strictly descending weights that sum to 1.
    """
    q = np.asarray(q, dtype=np.float64)
    if q.shape != (4,):
        raise ValueError("q must have shape (4,)")

    # 1) positive gaps  d_i = softplus(q_i)
    d = softplus(q)

    # 2) cumulative barriers  b_0 = 0,  b_i = Σ_{k<i} d_k   (i = 1..4)
    b = np.concatenate(([0.0], np.cumsum(d)))           # shape (5,)

    # 3) descending logits  a_i = -b_i   ⇒  a_0 = 0 > a_1 > … > a_4
    logits = -b

    # 4) soft‑max
    exp_logits = np.exp(logits - logits.max())           # subtract max for stability
    f = exp_logits / exp_logits.sum()

    return f

def f_to_q(f, eps=1e-200):
    """
    Parameters
    ----------
    f : array_like, shape (5,)
        Strictly descending positive probabilities that sum to one.

    Returns
    -------
    q : ndarray, shape (4,)
        Parameters that reproduce `f` under q_to_f.
    """
    f = np.asarray(f, dtype=np.float64)
    if f.shape != (5,):
        raise ValueError("f must have shape (5,)")

    if not (np.all(f[:-1] > f[1:]) and np.all(f > 0) and abs(f.sum() - 1.0) < 1e-10):
        raise ValueError("f must be positive, strictly descending and sum to 1")

    #     d_i = -ln(f_{i+1}/f_i) = ln(f_i/f_{i+1})   (positive)
    #d = np.log(f[:-1] / f[1:] + eps)
    d = np.log(np.clip(f[:-1] / f[1:], eps, None))

    #     q_i = softplus^{-1}(d_i)  =  ln(e^{d_i} - 1)
    #q = np.log(np.exp(d) - 1.0 + eps)
    q = np.log(np.clip(np.exp(d) - 1.0, eps, None))

    return q




def logistic(x):
    return 1.0 / (1.0 + np.exp(-x))

def logistic_inv(u):
    return np.log(u/(1. - u))

def inverse_stick_breaking(T, T_min, T_max):
    """
    T: array of length 5 (T0..T4), each in [T_min, T_max], strictly ascending
    T_min, T_max: scalars
    Returns: p array of length 5
    """
    # 1) b_i = ln(T_i)
    log_T = np.log(T)

    # 2) alpha_i in [0,1]
    log_T_min = np.log(T_min)
    log_T_max = np.log(T_max)

    alpha = (log_T - log_T_min) / (log_T_max - log_T_min)

    # 3) invert recursion alpha -> u
    u = np.zeros(5)
    u[0] = alpha[0]
    for i in range(1, 5):
        numerator = alpha[i] - alpha[i-1]
        denom     = 1.0 - alpha[i-1]
        u[i]      = numerator / denom

    # 4) p_i = logit(u_i)
    p = logistic_inv(u)
    return p

def transform_params_to_T(x, T_min, T_max):
    """
    x: array of length >=9
       x[:4] are a_i, x[4:] are p_i
    T_min, T_max: e.g. (T_input/2, 500)

    Returns: T array of length 5 guaranteed in [T_min, T_max] and sorted
    """
    # Extract the p_i
    p = x[4:]  # length 5
    u = logistic(p)  # each in (0,1)

    alpha = np.zeros(5)
    alpha[0] = u[0]
    for i in range(1,5):
        alpha[i] = alpha[i-1] + u[i]*(1 - alpha[i-1])

    # Now map alpha_i to b_i in [ln(T_min), ln(T_max)]
    log_T_min = np.log(T_min)
    log_T_max = np.log(T_max)
    b = log_T_min + alpha * (log_T_max - log_T_min)  # length 5

    # return T = exp(b)
    return np.exp(b)
# ------------------------------
# 2. Five-Maxwellian model
# ------------------------------
def five_maxwellians(E, params, T_min, T_max):
    """
    E : 1‑D array of energies
    params : 9‑vector  [q1..q4, p0..p4]
    """
    q  = params[:4]
    #p  = params[4:]

    f  = q_to_f(q)                          # shape (5,)
    T  = transform_params_to_T(params, T_min, T_max)   # shape (5,)

    E   = np.atleast_1d(E)                  # shape (m,)
    # log of each Maxwellian term
    log_terms = -E[:, None] / T[None, :]    # (m,5)

    # keep every exponent ≥ –700 to avoid underflow to zero
    np.clip(log_terms, -700, None, out=log_terms)

    # matrix form:   model_j = Σ_i  f_i * exp(log_terms_ji)
    model = (f * np.exp(log_terms)).sum(axis=1)
    return model

# ------------------------------
# 3. residuals vector function
# ------------------------------
"""
def residuals_vec(params, E, data, T_min, T_max):
    EPS  = 1e-120                        # numerical floor (float64 safe)
    model = five_maxwellians(E, params, T_min, T_max)

    # 10 % error bar, but never smaller than EPS
    err   = np.maximum(0.1 * data, EPS)

    return (model - data) / err         # ±1  <=>  ±10 % linear error
"""

LN10 = np.log(1.10)        # ln(1 + 10 %)  ≃ 0.09531

def resid_vec(params, E, data, T_min, T_max):
    model = five_maxwellians(E, params, T_min, T_max)
    EPS   = 1e-200                        # ≥ exp(‑700)
    # log residual; ±1  ↔  ±10 % relative error
    r = (np.log(np.maximum(model, EPS)) -
         np.log(np.maximum(data , EPS))) / LN10
    
    return r

# ------------------------------
# 4A. Fit routine: local solver (replacing differential evolution)
# ------------------------------
def fit_five_maxwellians(E, data, T_input, T_min, T_max):
    """
    """
    
    #bounds = [(-10.,10.)]*4 + [(-10.,10.)]*5
    #lower = np.array([-10.]*4 + [-10.]*5)
    #upper = np.array([ 10.]*4 + [ 10.]*5)
    lower = np.full(9, -10.0)
    upper = np.full(9,  10.0)

    #T_init = np.logspace(np.log10(T_input), np.log10(500), 5)
    if ((50. * T_input) < 500.):
        T_maxxx = 50. * T_input
    else:
        T_maxxx = 500.
    #T_init = np.logspace(np.log10(T_input/3), np.log10(3*T_input), 5)
    T_init = np.logspace(np.log10(T_input), np.log10(T_maxxx), 5)
    
    p_init = inverse_stick_breaking(T_init, T_min, T_max)
    
    # f_init = np.array([0.884, 0.1  , 0.01 , 0.005, 0.001])
    f_init = np.array([0.4, 0.25, 0.18, 0.12, 0.05])
    q_init = f_to_q(f_init) 
    x0 = np.concatenate([q_init, p_init])
    # Local solve
    # ----- run least_squares ------------------------------------------------
    result = least_squares(
        fun      = resid_vec,
        x0       = x0,
        args     = (E, data, T_min, T_max),   # True → log residuals
        bounds   = (lower, upper),
        method   = 'trf',
        jac      = '3-point',    # keep finite‑diff for now; exact Jacobian optional
        x_scale  = 'jac',        # scale variables automatically
        loss     = 'linear',   
        f_scale  = 1.0,         # 0.10 ⇒ residuals |r|>0.10 start to be down‑weighted
        ftol     = 1e-15,
        xtol     = 1e-15,
        gtol     = 1e-15,
        max_nfev = 100000,
        verbose  = 0            # 0 = silent, 1 = summary, 2 = per‑iter
    )

    #print("Converged?", result.success, " |  cost =", result.cost,
    #      " |  max|r| =", np.max(np.abs(result.fun)))
    return result.x, result.cost




def extract_f_T(params, T_min, T_max):
    q = params[:4]
    f = q_to_f(q)  
    T = transform_params_to_T(params, T_min, T_max)
    return f, T

def compute_percent_errors(E, data, params, T_min, T_max):
    model = five_maxwellians(E, params, T_min, T_max)
    eps = 1e-200
    pct_err = 100.0 * np.abs(model - data) / np.maximum(data, eps)
    return pct_err

def summary_stats(name, pct_err):
    mn = np.min(pct_err)
    mx = np.max(pct_err)
    mean_ = np.mean(pct_err)
    std_ = np.std(pct_err)
    q25 = np.percentile(pct_err, 25)
    q50 = np.median(pct_err)
    q75 = np.percentile(pct_err, 75)
    mode_ = stats.mode(np.round(pct_err, decimals=3), keepdims=True)[0][0]
    
    print(f"--- {name} Stats on % absolute relative error ---",flush=True)
    print(f"  Min:  {mn:10.4g}",flush=True)
    print(f"  Max:  {mx:10.4g}",flush=True)
    print(f"  Mean: {mean_:10.4g}",flush=True)
    print(f"  Std:  {std_:10.4g}",flush=True)
    print(f"  25%:  {q25:10.4g}",flush=True)
    print(f"  50%:  {q50:10.4g}  (median)",flush=True)
    print(f"  75%:  {q75:10.4g}",flush=True)
    print(f"  Mode(approx): {mode_:10.4g}\n",flush=True)







if __name__=="__main__":
    start_time = time.time()
    # ------------------------------
    # Example usage or main loop
    # ------------------------------
    T_min = 0.99
    T_max = 100.
    num_points_T = 200
    T = np.logspace(np.log10(T_min), np.log10(T_max), num_points_T)
    
    
    
    
    kappa_min = 1.51
    kappa_max = 300.
    num_points_kappa = 200
    kappa = np.logspace(np.log10(kappa_min), np.log10(kappa_max), num_points_kappa)
    
    kappa_min2 = 1.01
    kappa2 = np.logspace(np.log10(kappa_min2), np.log10(kappa_max), num_points_kappa)

    E_min = 0.01
    E_max = 500.0
    num_points_E = 150

    T_i_max4_sk = np.full((5,len(T),len(kappa)),np.nan)
    f_i_max4_sk = np.full((5,len(T),len(kappa)),np.nan)

    T_i_max4_pk = np.full((5,len(T),len(kappa)),np.nan)
    f_i_max4_pk = np.full((5,len(T),len(kappa)),np.nan)

    T_i_max4_fe = np.full((5,len(T),len(kappa)),np.nan)
    f_i_max4_fe = np.full((5,len(T),len(kappa)),np.nan)

    max_rel_err_sk = np.full((len(T),len(kappa)),np.nan)
    max_rel_err_pk = np.full((len(T),len(kappa)),np.nan)
    max_rel_err_fe = np.full((len(T),len(kappa)),np.nan)

    for i in range(len(T)):
        E_vals =  np.logspace(np.log10(E_min), np.log10(min(E_max,100.*T[i])), num_points_E)
        print(" i = ",i, f" of {num_points_T - 1}",flush=True)
        for j in range(len(kappa)):
            #print(" i+1 = ",i+1, f" of {num_points_T}, j+1 = ",j+1, f" of {num_points_kappa}")

            data_standard  = np.array([standard_kappa(E, T[i], kappa[j]) for E in E_vals])
            data_product   = np.array([product_kappa(E,  T[i], kappa2[j]) for E in E_vals])
            data_fried_egg = np.array([fried_egg(E,      T[i], kappa2[j]) for E in E_vals])
            """
            if np.min(data_standard) < 0.0:
                print("Warning Standard Kappa < 0, min =", np.min(data_standard),flush=True)
            if np.min(data_product) < 0.0:
                print("Warning Product Kappa < 0, min =", np.min(data_product),flush=True)
            if np.min(data_fried_egg) < 0.0:
                print("Warning Fried Egg < 0, min =", np.min(data_fried_egg),flush=True)
            """
            T_min, T_max = 0.01, 750. # make if lower T_max below 401 make sure lower T_init[4] in fit_five_maxwellians above
            # Fit standard kappa
            params_std, cost_std = fit_five_maxwellians(E_vals, data_standard, T[i], T_min, T_max)
            f_std, T_std = extract_f_T(params_std, T_min, T_max)
            T_i_max4_sk[:, i, j] = T_std
            f_i_max4_sk[:, i, j] = f_std
    
            # Fit product kappa
            params_prod, cost_prod = fit_five_maxwellians(E_vals, data_product, T[i], T_min, T_max)
            f_prod, T_prod = extract_f_T(params_prod, T_min, T_max)
            T_i_max4_pk[:, i, j] = T_prod
            f_i_max4_pk[:, i, j] = f_prod
    
            # Fit fried-egg
            params_fried, cost_fried = fit_five_maxwellians(E_vals, data_fried_egg, T[i], T_min, T_max)
            f_fried, T_fried = extract_f_T(params_fried, T_min, T_max)
            T_i_max4_fe[:, i, j] = T_fried
            f_i_max4_fe[:, i, j] = f_fried
    
            # % error
            pct_err_std   = compute_percent_errors(E_vals, data_standard,  params_std, T_min, T_max)
            pct_err_prod  = compute_percent_errors(E_vals, data_product,   params_prod, T_min, T_max)
            pct_err_fried = compute_percent_errors(E_vals, data_fried_egg, params_fried, T_min, T_max)
    
            max_rel_err_sk[i, j]  = np.max(pct_err_std)
            max_rel_err_pk[i, j]  = np.max(pct_err_prod)
            max_rel_err_fe[i, j]  = np.max(pct_err_fried)
    
            #  if >20% error give warning
            """
            if max_rel_err_sk[i, j] > 20.0:
                print("Warning Max Percent Error Standard-Kappa > 10%:", max_rel_err_sk[i, j],flush=True)

    
            if max_rel_err_pk[i, j] > 20.0:
                print("Warning Max Percent Error Product-Kappa > 10%:", max_rel_err_pk[i, j],flush=True)

    
            if max_rel_err_fe[i, j] > 20.0:
                print("Warning Max Percent Error Fried-Egg > 10%:", max_rel_err_fe[i, j],flush=True)
            """
            """
            for k in range(4):
                if (T_i_max4_sk[k, i, j] > T_i_max4_sk[k + 1, i, j] + 1e-7):
                    print(f"Standard kappa Temp Constraints not Met for T {k}",flush=True)
                if (T_i_max4_pk[k, i, j] > T_i_max4_pk[k + 1, i, j] + 1e-7):
                    print(f"Product kappa Temp Constraints not Met for T {k}",flush=True)
                if (T_i_max4_fe[k, i, j] > T_i_max4_fe[k + 1, i, j] + 1e-7):
                    print(f"Fried-Egg Temp Constraints not Met for T {k}",flush=True)
                if (f_i_max4_sk[k, i, j] < f_i_max4_sk[k + 1, i, j] - 1e-7):
                    print(f"Standard kappa fi Constraints not Met for f {k}",flush=True)
                if (f_i_max4_pk[k, i, j] < f_i_max4_pk[k + 1, i, j] - 1e-7):
                    print(f"Product kappa fi Constraints not Met for f {k}",flush=True)
                if (f_i_max4_fe[k, i, j] < f_i_max4_fe[k + 1, i, j] - 1e-7):
                    print(f"Fried-Egg fi Constraints not Met for f {k}",flush=True)
            """
    # Save
    np.savez(
        "lsquares_fixed_order_fi_and_Ti_local_solver_100x100_T_vs_kappa_indiv_over_Elog150pnts_0.01-min_of_500_or_100xT.npz",
        T=T,
        kappa=kappa,  
        kappa2=kappa2,  
        E_vals=E_vals,
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
    
    print('nanmax of max rel err for  sk, pk,and  fe = ',np.nanamax(max_rel_err_sk),np.nanamax(max_rel_err_pk),np.nanamax(max_rel_err_fe),flush=True)
    print("--- %s seconds ---" % (time.time() - start_time),flush=True)
    """
    # Plotting (unchanged)
    for j in range(len(kappa)):
        # ------------------------------------------------------------
        # 1) Figure for T_i and f_i in 3 rows x 2 columns
        # ------------------------------------------------------------
        fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(12, 10), sharex=True)
    
        # ============ Row 0: Standard Kappa ============
        axes[0,0].plot(T, T_i_max4_sk[0, :, j], label=r'T$_1$')
        axes[0,0].plot(T, T_i_max4_sk[1, :, j], label=r'T$_2$')
        axes[0,0].plot(T, T_i_max4_sk[2, :, j], label=r'T$_3$')
        axes[0,0].plot(T, T_i_max4_sk[3, :, j], label=r'T$_4$')
        axes[0,0].plot(T, T_i_max4_sk[4, :, j], label=r'T$_5$')
        axes[0,0].set_yscale('log')
        axes[0,0].set_ylim(0.01, 750)
        axes[0,0].set_ylabel(r"T$_i$ (eV)")
        axes[0,0].set_title(f"Standard Kappa T$_i$, kappa={kappa[j]:.2f}")
        axes[0,0].legend()
    
        axes[0,1].plot(T, f_i_max4_sk[0, :, j], label=r'f$_1$')
        axes[0,1].plot(T, f_i_max4_sk[1, :, j], label=r'f$_2$')
        axes[0,1].plot(T, f_i_max4_sk[2, :, j], label=r'f$_3$')
        axes[0,1].plot(T, f_i_max4_sk[3, :, j], label=r'f$_4$')
        axes[0,1].plot(T, f_i_max4_sk[4, :, j], label=r'f$_5$')
        axes[0,1].set_yscale('log')
        axes[0,1].set_ylim(0.001, 1)
        axes[0,1].set_ylabel(r"f$_i$")
        axes[0,1].set_title(f"Standard Kappa f$_i$, kappa={kappa[j]:.2f}")
        axes[0,1].legend()
    
        # ============ Row 1: Product Kappa ============
        axes[1,0].plot(T, T_i_max4_pk[0, :, j], label=r'T$_1$')
        axes[1,0].plot(T, T_i_max4_pk[1, :, j], label=r'T$_2$')
        axes[1,0].plot(T, T_i_max4_pk[2, :, j], label=r'T$_3$')
        axes[1,0].plot(T, T_i_max4_pk[3, :, j], label=r'T$_4$')
        axes[1,0].plot(T, T_i_max4_pk[4, :, j], label=r'T$_5$')
        axes[1,0].set_yscale('log')
        axes[1,0].set_ylim(0.01, 750)
        axes[1,0].set_ylabel(r"T$_i$ (eV)")
        axes[1,0].set_title(f"Product Kappa T$_i$, kappa={kappa2[j]:.2f}")
        axes[1,0].legend()
    
        axes[1,1].plot(T, f_i_max4_pk[0, :, j], label=r'f$_1$')
        axes[1,1].plot(T, f_i_max4_pk[1, :, j], label=r'f$_2$')
        axes[1,1].plot(T, f_i_max4_pk[2, :, j], label=r'f$_3$')
        axes[1,1].plot(T, f_i_max4_pk[3, :, j], label=r'f$_4$')
        axes[1,1].plot(T, f_i_max4_pk[4, :, j], label=r'f$_5$')
        axes[1,1].set_yscale('log')
        axes[1,1].set_ylim(0.001, 1)
        axes[1,1].set_ylabel(r"f$_i$")
        axes[1,1].set_title(f"Product Kappa f$_i$, kappa={kappa2[j]:.2f}")
        axes[1,1].legend()
    
        # ============ Row 2: Fried-Egg ============
        axes[2,0].plot(T, T_i_max4_fe[0, :, j], label=r'T$_1$')
        axes[2,0].plot(T, T_i_max4_fe[1, :, j], label=r'T$_2$')
        axes[2,0].plot(T, T_i_max4_fe[2, :, j], label=r'T$_3$')
        axes[2,0].plot(T, T_i_max4_fe[3, :, j], label=r'T$_4$')
        axes[2,0].plot(T, T_i_max4_fe[4, :, j], label=r'T$_5$')
        axes[2,0].set_yscale('log')
        axes[2,0].set_ylim(0.01, 750)
        axes[2,0].set_ylabel(r"T$_i$ (eV)")
        axes[2,0].set_xlabel(r"T (eV)")
        axes[2,0].set_title(f"Fried-Egg T$_i$, kappa={kappa2[j]:.2f}")
        axes[2,0].legend()
    
        axes[2,1].plot(T, f_i_max4_fe[0, :, j], label=r'f$_1$')
        axes[2,1].plot(T, f_i_max4_fe[1, :, j], label=r'f$_2$')
        axes[2,1].plot(T, f_i_max4_fe[2, :, j], label=r'f$_3$')
        axes[2,1].plot(T, f_i_max4_fe[3, :, j], label=r'f$_4$')
        axes[2,1].plot(T, f_i_max4_fe[4, :, j], label=r'f$_5$')
        axes[2,1].set_yscale('log')
        axes[2,1].set_ylim(0.001, 1)
        axes[2,1].set_ylabel(r"f$_i$")
        axes[2,1].set_xlabel(r"T (eV)")
        axes[2,1].set_title(f"Fried-Egg f$_i$, kappa={kappa2[j]:.2f}")
        axes[2,1].legend()
    
        plt.tight_layout()
        plt.savefig(f"lsquares_fixed_order_fi_and_Ti_local_solver_10x10_T_kappa_Ti_fi_subplots_5MaxFits_kappa={kappa[j]:.2f}_kappa2={kappa2[j]:.2f}_Elog150pnts_0.01-min_of_500_or_100xT.pdf", bbox_inches='tight')
        plt.show()
    
        # ------------------------------------------------------------
        # 2) Figure for max_abs_percent_error in 3 rows × 1 column
        # ------------------------------------------------------------
        fig2, axes2 = plt.subplots(nrows=3, ncols=1, figsize=(8, 10), sharex=True)
    
        # Standard Kappa
        axes2[0].plot(T, max_rel_err_sk[:, j], label='Standard-Kappa')
        axes2[0].set_yscale('log')
        axes2[0].set_ylim(0.1, 20)
        axes2[0].set_ylabel("Max % Rel. Error")
        axes2[0].set_title(f"Max % Error vs. T, Standard Kappa, kappa={kappa[j]:.2f}")
        axes2[0].legend()
    
        # Product Kappa
        axes2[1].plot(T, max_rel_err_pk[:, j], label='Product-Kappa')
        axes2[1].set_yscale('log')
        axes2[1].set_ylim(0.1, 20)
        axes2[1].set_ylabel("Max % Rel. Error")
        axes2[1].set_title(f"Max % Error vs. T, Product Kappa, kappa={kappa2[j]:.2f}")
        axes2[1].legend()
    
        # Fried-Egg
        axes2[2].plot(T, max_rel_err_fe[:, j], label='Fried-Egg')
        axes2[2].set_yscale('log')
        axes2[2].set_ylim(0.1, 20)
        axes2[2].set_xlabel(r"T (eV)")
        axes2[2].set_ylabel("Max % Rel. Error")
        axes2[2].set_title(f"Max % Error vs. T, Fried-Egg, kappa={kappa2[j]:.2f}")
        axes2[2].legend()
    
        plt.tight_layout()
        plt.savefig(f"lsquares_fixed_order_fi_and_Ti_local_solver_10x10max_rel_error_subplots_5MaxFits_kappa={kappa[j]:.2f}_kappa2={kappa2[j]:.2f}_Elog150pnts_0.01-min_of_500_or_100xT.pdf", bbox_inches='tight')
        plt.show()

    """