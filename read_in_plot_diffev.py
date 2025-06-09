# -*- coding: utf-8 -*-
"""
Created on Wed Apr 16 16:30:01 2025

@author: Owner
"""


import numpy as np
import matplotlib.pyplot as plt

from scipy.special import gamma
from scipy.integrate import quad
from scipy.optimize import minimize, LinearConstraint, differential_evolution
from scipy import stats

import pandas as pd
import mpmath as mp


# ------------------------------
# 1. Define the "true" distributions
# ------------------------------
def standard_kappa(E, T, kappa):
    """
    Standard kappa distribution:
      f(E) = (Gamma(kappa) / [sqrt(kappa) * Gamma(kappa - 1/2)])
             * (1 + E/(kappa*T))^(-(1 + kappa))
    """
    nnn=0.
    facc=(E**nnn)
    if kappa < 50.:
        prefactor = facc*gamma(kappa) / (np.sqrt(kappa) * gamma(kappa - 0.5))
        factorr = (1.0 + E/(kappa*T))**(-(1.0 + kappa))
        result = prefactor * factorr
    else:
        mp.mp.dps = 30
        prefactor = mp.mpf(facc)*mp.gamma(kappa) / (mp.sqrt(kappa) * mp.gamma(kappa - 0.5))
        factorr = mp.mpf((1.0 + E/(kappa*T))**(-(1.0 + kappa)))
        result = float(prefactor * factorr)
    
    return result

def product_kappa(E, T, kappa):
    """
    Product kappa distribution:
      f(E) = [Gamma(kappa+1)/( sqrt(kappa)*Gamma(kappa+0.5))]
             * Integral_{0 to 1} (1 + (E/(kappa*T))*(1 - t^2))^(-(1 + kappa))
                                 * (1 + (E/(kappa*T))*t^2)^(-(1 + kappa)) dt
    """
    nnn=0.
    facc=(E**nnn)
    a2 = E / (T * kappa)
    prefactor = facc*mp.gamma(kappa + 1.0) / (mp.sqrt(kappa) * mp.gamma(kappa + 0.5))
    
    def integrand_pk(t):
        return ((1.0 + a2*(1.0 - t**2))**(-(1.0 + kappa))) * \
               ((1.0 + a2*(t**2))**(-(1.0 + kappa)))
                       
    integral_val, err = quad(integrand_pk, 0.0, 1.0, limit=100)
    per_err = 100.0 * abs(err / integral_val)
    if float(per_err) > 0.1:
        mp.mp.dps = 50
        integral_val, err = mp.quad(integrand_pk, [0.0, 1.0], error=True)
        per_err = 100.0 * abs(err / integral_val)
        if float(per_err) > 0.1:
            mp.mp.dps = 200
            integral_val, err = mp.quad(integrand_pk, [0.0, 1.0], error=True, maxdegree=10)
            per_err = 100.0 * abs(err / integral_val)
            if float(per_err) > 0.1:
                print("Warning: Product Kappa Percent error2 for maxdegree=100 =", per_err,
                      " for T =", T, " and kappa =", kappa, flush=True)
    return float(prefactor * integral_val)

def fried_egg(E, T, kappa):
    """
    Fried-Egg distribution:
      f(E) = Integral_{0 to 1} [
                (1 + (E/(kappa*T))*(1 - t**2))^(-(1 + kappa)) * exp(-E*t^2 / T)
             ] dt
    """
    
    nnn=0.
    facc=(E**nnn)
    def integrand_fe(t):
        return ((1.0 + (E/(kappa*T))*(1.0 - t**2))**(-(1.0 + kappa))) * \
                mp.e**(-E*(t**2)/T)
    
    integral_val, err = quad(integrand_fe, 0.0, 1.0, limit=100)
    per_err = 100.0 * abs(err / integral_val)
    if float(per_err) > 0.1:
        mp.mp.dps = 50
        integral_val, err = mp.quad(integrand_fe, [0.0, 1.0], error=True)
        per_err = 100.0 * abs(err / integral_val)
        if float(per_err) > 0.1:
            mp.mp.dps = 200
            integral_val, err = mp.quad(integrand_fe, [0.0, 1.0], error=True, maxdegree=10)
            per_err = 100.0 * abs(err / integral_val)
            if float(per_err) > 0.1:
                print("Warning: Fried-Egg Percent error2 for maxdegree=100 =", per_err,
                      " for T =", T, " and kappa =", kappa, flush=True)
    return float(facc*integral_val)




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

    # Option (A): just sort them
    # alpha = np.sort(u)  # ensures alpha_0 <= alpha_1 <= ... <= alpha_4

    # Option (B): "stick-breaking" approach, fully differentiable
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
   params: 
     - first 4 entries are a_1..a_4 (with a_0 pinned = 0)
     - next 5 entries are p_0..p_4 for stick-breaking
   ...
   """
    a = [0., params[0], params[1], params[2], params[3]]
    exp_a = np.exp(a)
    sum_exp_a = np.sum(exp_a)
    f = exp_a / sum_exp_a 
    T = transform_params_to_T(params, T_min, T_max)
    nnn=0.
    facc=(E**nnn)
    return (facc*(
        f[0]*np.exp(-np.clip(E/T[0], a_min=0., a_max=700.)) +
        f[1]*np.exp(-np.clip(E/T[1], a_min=0., a_max=700.)) +
        f[2]*np.exp(-np.clip(E/T[2], a_min=0., a_max=700.)) +
        f[3]*np.exp(-np.clip(E/T[3], a_min=0., a_max=700.)) +
        f[4]*np.exp(-np.clip(E/T[4], a_min=0., a_max=700.)))
    )


# ------------------------------
# 3. Cost function
# ------------------------------
def cost_function_de_global(params, E, data, T_min, T_max):
    """
    Modified cost function that:
      1) Computes the standard sum of squared log-residuals.
      2) Computes max % absolute relative error across energies.
      3) Adds a penalty on max error and on T ordering.
    """
    EPS = 1e-100
    model_vals = five_maxwellians(E, params, T_min, T_max)
    model_vals_safe = np.maximum(model_vals, EPS)
    data_safe = np.maximum(data, EPS)
    
    # --- base cost (sum of squares of log-differences) ---
    residuals = np.log(model_vals_safe) - np.log(data_safe)
    base_cost = np.sum(residuals**2.)
    
    # --- compute % errors for penalty check ---
    pct_err_arr = 100.0 * np.abs(model_vals_safe - data_safe)/np.maximum(data_safe, EPS)
    max_err     = np.max(pct_err_arr)      
    
    # define penalty 
    alpha = 5e-4
    threshold   = 10.0
    penalty = alpha* np.sum( np.maximum( pct_err_arr - threshold, 0.0)**2.)
        
    return base_cost + penalty


# ------------------------------
# 4A. Fit routine: local solver (replaced here with differential evolution)
# ------------------------------
def fit_five_maxwellians(E, data, T_input, T_min, T_max):
    """
    Fit sum of 5 Maxwellians using a local solver (L-BFGS-B).
    -30 <= a_i <= 30
    ln(0.01) <= b_i <= ln(500)
    """
    
    # Build the bounds
    bounds = [(-10.,10.)]*4 + [(-10.,10.)]*5
    
    # Example initial guess (not used by differential_evolution but left for reference):
    T_init = np.logspace(np.log10(T_input), np.log10(500), 5)
    p_init = inverse_stick_breaking(T_init, T_min, T_max)
    a_init = [0.]*4
    x0 = list(a_init) + list(p_init)
    
    # Global solve using differential evolution
    result = differential_evolution(
        cost_function_de_global,
        bounds=bounds,
        args=(E, data, T_min, T_max),
        maxiter=20000,
        strategy='rand1bin',#strategy='best1bin',
        polish=True,
        tol=1e-5,
        mutation=(0.5, 1),
        recombination=0.7,
        disp=False,
        workers=-1,
        updating='deferred',
        popsize=50,
        init='sobol'#'latinhypercube'
    )
    print("Converged?", result.success)
    print("Message:", result.message)
    return result.x, result.fun



def extract_f_T(params, T_min, T_max):
    a = [0., params[0], params[1], params[2], params[3]]
    exp_a = np.exp(a)
    sum_exp_a = np.sum(exp_a)
    f = exp_a / sum_exp_a
    T = transform_params_to_T(params, T_min, T_max)
    return f, T

def compute_percent_errors(E, data, params, T_min, T_max):
    model = five_maxwellians(E, params, T_min, T_max)
    eps = 1e-30
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
    
    print(f"--- {name} Stats on % absolute relative error ---")
    print(f"  Min:  {mn:10.4g}")
    print(f"  Max:  {mx:10.4g}")
    print(f"  Mean: {mean_:10.4g}")
    print(f"  Std:  {std_:10.4g}")
    print(f"  25%:  {q25:10.4g}")
    print(f"  50%:  {q50:10.4g}  (median)")
    print(f"  75%:  {q75:10.4g}")
    print(f"  Mode(approx): {mode_:10.4g}\n")


if __name__=="__main__":
    # Load the .npz archive
    data = np.load("diffev_10x10_T_vs_kappa_indiv_dif_ev_fully_contrained_fits_of_Evals.npz")
    
    # Extract each array into a variable of the same name
    T               = data["T"]
    kappa           = data["kappa"]
    E_vals          = data["E_vals"]
    max_rel_err_sk  = data["max_rel_err_sk"]
    max_rel_err_pk  = data["max_rel_err_pk"]
    max_rel_err_fe  = data["max_rel_err_fe"]
    T_i_max4_sk     = data["T_i_max4_sk"]
    T_i_max4_pk     = data["T_i_max4_pk"]
    T_i_max4_fe     = data["T_i_max4_fe"]
    f_i_max4_sk     = data["f_i_max4_sk"]
    f_i_max4_pk     = data["f_i_max4_pk"]
    f_i_max4_fe     = data["f_i_max4_fe"]
    
    # (Optional) close the archive when done
    data.close()

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
        axes[0,0].set_ylim(0.1, 500)
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
        axes[1,0].set_ylim(0.01, 500)
        axes[1,0].set_ylabel(r"T$_i$ (eV)")
        axes[1,0].set_title(f"Product Kappa T$_i$, kappa={kappa[j]:.2f}")
        axes[1,0].legend()
    
        axes[1,1].plot(T, f_i_max4_pk[0, :, j], label=r'f$_1$')
        axes[1,1].plot(T, f_i_max4_pk[1, :, j], label=r'f$_2$')
        axes[1,1].plot(T, f_i_max4_pk[2, :, j], label=r'f$_3$')
        axes[1,1].plot(T, f_i_max4_pk[3, :, j], label=r'f$_4$')
        axes[1,1].plot(T, f_i_max4_pk[4, :, j], label=r'f$_5$')
        axes[1,1].set_yscale('log')
        axes[1,1].set_ylim(0.001, 1)
        axes[1,1].set_ylabel(r"f$_i$")
        axes[1,1].set_title(f"Product Kappa f$_i$, kappa={kappa[j]:.2f}")
        axes[1,1].legend()
    
        # ============ Row 2: Fried-Egg ============
        axes[2,0].plot(T, T_i_max4_fe[0, :, j], label=r'T$_1$')
        axes[2,0].plot(T, T_i_max4_fe[1, :, j], label=r'T$_2$')
        axes[2,0].plot(T, T_i_max4_fe[2, :, j], label=r'T$_3$')
        axes[2,0].plot(T, T_i_max4_fe[3, :, j], label=r'T$_4$')
        axes[2,0].plot(T, T_i_max4_fe[4, :, j], label=r'T$_5$')
        axes[2,0].set_yscale('log')
        axes[2,0].set_ylim(0.01, 500)
        axes[2,0].set_ylabel(r"T$_i$ (eV)")
        axes[2,0].set_xlabel(r"T (eV)")
        axes[2,0].set_title(f"Fried-Egg T$_i$, kappa={kappa[j]:.2f}")
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
        axes[2,1].set_title(f"Fried-Egg f$_i$, kappa={kappa[j]:.2f}")
        axes[2,1].legend()
    
        plt.tight_layout()
        plt.savefig(f"diffev_10x10_T_kappa_Ti_fi_subplots_5MaxFits_kappa={kappa[j]:.2f}.pdf", bbox_inches='tight')
        plt.show()
    
        # ------------------------------------------------------------
        # 2) Figure for max_abs_percent_error in 3 rows × 1 column
        # ------------------------------------------------------------
        fig2, axes2 = plt.subplots(nrows=3, ncols=1, figsize=(8, 10), sharex=True)
    
        # Standard Kappa
        axes2[0].plot(T, max_rel_err_sk[:, j], label='Standard-Kappa')
        axes2[0].set_yscale('log')
        axes2[0].set_ylim(0.1, 50)
        axes2[0].set_ylabel("Max % Rel. Error")
        axes2[0].set_title(f"Max % Error vs. T, Standard Kappa, kappa={kappa[j]:.2f}")
        axes2[0].legend()
    
        # Product Kappa
        axes2[1].plot(T, max_rel_err_pk[:, j], label='Product-Kappa')
        axes2[1].set_yscale('log')
        axes2[1].set_ylim(0.1, 50)
        axes2[1].set_ylabel("Max % Rel. Error")
        axes2[1].set_title(f"Max % Error vs. T, Product Kappa, kappa={kappa[j]:.2f}")
        axes2[1].legend()
    
        # Fried-Egg
        axes2[2].plot(T, max_rel_err_fe[:, j], label='Fried-Egg')
        axes2[2].set_yscale('log')
        axes2[2].set_ylim(0.1, 50)
        axes2[2].set_xlabel(r"T (eV)")
        axes2[2].set_ylabel("Max % Rel. Error")
        axes2[2].set_title(f"Max % Error vs. T, Fried-Egg, kappa={kappa[j]:.2f}")
        axes2[2].legend()
    
        plt.tight_layout()
        plt.savefig(f"diffev_10x10max_rel_error_subplots_5MaxFits_kappa={kappa[j]:.2f}.pdf", bbox_inches='tight')
        plt.show()