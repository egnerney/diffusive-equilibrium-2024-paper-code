# -*- coding: utf-8 -*-
"""
Created on Sat Apr 12 02:19:09 2025

@author: Owner
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Apr 11 07:26:05 2025

@author: Owner
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 11 00:14:29 2025

@author: edne8319
"""

import numpy as np
import matplotlib.pyplot as plt

from scipy.special import gamma
from scipy.integrate import quad
from scipy.optimize import differential_evolution, LinearConstraint


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
    if kappa< 50.:
        prefactor = gamma(kappa) / (np.sqrt(kappa) *gamma(kappa - 0.5))
        factorr = (1.0 + E/(kappa*T))**(-(1.0 + kappa))
        result= prefactor * factorr
    else:
        mp.mp.dps = 30
        prefactor = mp.gamma(kappa) / (mp.sqrt(kappa) * mp.gamma(kappa - 0.5))
        factorr = mp.mpf((1.0 + E/(kappa*T))**(-(1.0 + kappa)))
        result= float(prefactor * factorr)
    
    return result

def product_kappa(E, T, kappa):
    """
    Product kappa distribution:
      f(E) = [Gamma(kappa+1)/( sqrt(kappa)*Gamma(kappa+0.5))]
             * Integral_{0 to 1} (1 + (E/(kappa*T))*(1 - t^2))^(-(1 + kappa))
                                 * (1 + (E/(kappa*T))*t^2)^(-(1 + kappa)) dt
    """
    #mp.mp.dps = 30
    a2 = E / (T * kappa)
    prefactor = mp.gamma(kappa + 1.0) / (mp.sqrt(kappa) * mp.gamma(kappa + 0.5))
    
    def integrand_pk(t):
        return ((1.0 + a2*(1.0 - t**2))**(-(1.0 + kappa))) * \
               ((1.0 + a2*(t**2))**(-(1.0 + kappa)))
                       
    integral_val, err = quad(integrand_pk,0.0, 1.0, limit=100)
    per_err = 100.0 * abs(err / integral_val)
    if float(per_err) > 0.1:
        mp.mp.dps = 50
        #print("Warning: Product Kappa Percent error =", per_err," for T =", T, " and kappa =", kappa, flush=True)
        integral_val, err = mp.quad(integrand_pk, [0.0, 1.0], error=True)
        per_err = 100.0 * abs(err / integral_val)
        if float(per_err) > 0.1:
            mp.mp.dps = 200
            integral_val, err = mp.quad(integrand_pk, [0.0, 1.0], error=True, maxdegree=10)
            per_err = 100.0 * abs(err / integral_val)
            if float(per_err) > 0.1:
                print("Warning: Product Kappa Percent error2 for maxdegree=100 =", per_err," for T =", T, " and kappa =", kappa, flush=True)         
    return float(prefactor * integral_val)

def fried_egg(E, T, kappa):
    """
    Fried-Egg distribution:
      f(E) = Integral_{0 to 1} [
                (1 + (E/(kappa*T))*(1 - t**2))^(-(1 + kappa)) * exp(-E*t^2 / T)
             ] dt
    """
    #mp.mp.dps = 30
    def integrand_fe(t):
        return ((1.0 + (E/(kappa*T))*(1.0 - t**2))**(-(1.0 + kappa))) * \
                mp.e**(-E*(t**2)/T)
    
    integral_val, err = quad(integrand_fe, 0.0, 1.0, limit=100)
    per_err = 100.0 * abs(err / integral_val)
    if float(per_err) > 0.1:
        mp.mp.dps = 50
        #print("Warning: Fried-Egg Percent error =", per_err," for T =", T, " and kappa =", kappa, flush=True)
        integral_val, err = mp.quad(integrand_fe, [0.0, 1.0], error=True)
        per_err = 100.0 * abs(err / integral_val)
        if float(per_err) > 0.1:
            mp.mp.dps = 200
            #print("Warning: Fried-Egg Percent error =", per_err," for T =", T, " and kappa =", kappa, flush=True)
            integral_val, err = mp.quad(integrand_fe, [0.0, 1.0], error=True, maxdegree=10)
            per_err = 100.0 * abs(err / integral_val)
            if float(per_err) > 0.1:
                print("Warning: Fried-Egg Percent error2 for maxdegree=100 =", per_err," for T =", T, " and kappa =", kappa, flush=True)
    return float(integral_val)


# ------------------------------
# 2. Five-Maxwellian model
# ------------------------------
def five_maxwellians(E, params):
    """
    params: [a_1..a_5, b_1..b_5]
    f_i = exp(a_i)/sum_j exp(a_j)
    T_i = exp(b_i)
    sum of Maxwellians
    """
    a = params[:5]
    b = params[5:]
    
    exp_a = np.exp(a)
    sum_exp_a = np.sum(exp_a)
    f = exp_a / sum_exp_a
    T = np.exp(b)
    
    return (
        f[0]*np.exp(-E/T[0]) +
        f[1]*np.exp(-E/T[1]) +
        f[2]*np.exp(-E/T[2]) +
        f[3]*np.exp(-E/T[3]) +
        f[4]*np.exp(-E/T[4])
    )


# ------------------------------
# 3. Cost function (top-level!)
# ------------------------------

def cost_function_de_global(params, E, data):
    """
    Modified cost function that:
      1) Computes the standard sum of squared log-residuals.
      2) Computes max % absolute relative error across energies.
      3) Adds a smoothly-varying penalty if max error > 10%.
         The penalty's influence is tapered at the very beginning
         if the median % error is too large (e.g., >200%), so as not
         to dominate the initial search phase.
    """
    EPS = 1e-30
    model_vals = five_maxwellians(E, params)
    model_vals_safe = np.maximum(model_vals, EPS)
    data_safe       = np.maximum(data,       EPS)
    
    # --- base cost (sum of squares of log-differences) ---
    residuals = np.log(model_vals_safe) - np.log(data_safe)
    base_cost = np.sum(residuals**2)
    
    # --- compute % errors for penalty check ---
    pct_err_arr = 100.0 * np.abs(model_vals_safe - data_safe)/np.maximum(data_safe, EPS)
    max_err     = np.max(pct_err_arr)      # maximum % error across energy bins
    median_err  = np.median(pct_err_arr)   # median % error across energy bins
    
    # define penalty if max_err > threshold
    threshold_err = 5.0    # 10% threshold
    penalty = 0.0
    
    """
    if max_err > threshold_err:
        # If the median error is extremely large (e.g., > 200%),
        # then we reduce or eliminate the penalty to allow DE to
        # make big corrections first without being swamped by penalty.
        median_cutoff = 200.0
        if median_err < median_cutoff:
            # We'll let the penalty scale up smoothly as median_err
            # drops below median_cutoff. Example approach:
            #scale_factor = (median_cutoff - median_err) / median_cutoff
            # A small alpha factor so penalty doesn't completely dominate
            alpha = 0.2
            
            # Quadratic or exponential growth from the threshold:
            # e.g. (max_err - threshold_err)**2
            # penalty = alpha * scale_factor * (max_err - threshold_err)**2
            alpha_scaled = alpha * np.clip( (median_cutoff - median_err)/median_cutoff, 0.0, 0.5 )
            penalty = alpha_scaled * (max_err - threshold_err)**2
    """
    alpha = 0.5
    
    # Quadratic or exponential growth from the threshold:
    # e.g. (max_err - threshold_err)**2
    # penalty = alpha * scale_factor * (max_err - threshold_err)**2
    alpha_scaled = alpha #* np.clip( (median_cutoff - median_err)/median_cutoff, 0.0, 0.5 )
    penalty = alpha_scaled *( max_err**2.) #(max_err - threshold_err)**2            
            
    b = params[5:]
    T = np.exp(b)
    penalty_temp = 0.0
    alpha2 = 1. # 0.2, 0.5, 1.0, 10., 100.0
    for i in range(4):
        if T[i] > T[i+1]:
            diff = T[i] - T[i+1]
            penalty_temp += alpha2 * diff**2   # alpha>0
    
    return base_cost + penalty + penalty_temp


# ------------------------------
# 4. Fit routine with differential_evolution
# ------------------------------
def fit_five_maxwellians(E, data):
    """
    Fit sum of 5 Maxwellians using differential_evolution.
    -30 <= a_i <= 30
    ln(0.01) <= b_i <= ln(500)
    """
    # We'll build a single A matrix of shape (8,10):
    # rows 0..3 => b_i - b_{i+1} <= 0
    # rows 4..7 => a_i - a_{i+1} >= 0
    """
    A = np.zeros((4, 10))

    # Indices: a1->p[0], a2->p[1], a3->p[2], a4->p[3], a5->p[4]
    #          b1->p[5], b2->p[6], b3->p[7], b4->p[8], b5->p[9].

    # --- Enforce T1 <= T2 <= T3 <= T4 <= T5 ---
    # b1 - b2 <= 0
    A[0,5] =  1.0   # b1
    A[0,6] = -1.0   # b2

    # b2 - b3 <= 0
    A[1,6] =  1.0
    A[1,7] = -1.0

    # b3 - b4 <= 0
    A[2,7] =  1.0
    A[2,8] = -1.0

    # b4 - b5 <= 0
    A[3,8] =  1.0
    A[3,9] = -1.0

    # --- Enforce f1 >= f2 >= f3 >= f4 >= f5 ---
    # means a1 >= a2, a2 >= a3, etc.
    

    # a1 - a2 >= 0
    A[4,0] =  1.0
    A[4,1] = -1.0

    # a2 - a3 >= 0
    A[5,1] =  1.0
    A[5,2] = -1.0

    # a3 - a4 >= 0
    A[6,2] =  1.0
    A[6,3] = -1.0

    # a4 - a5 >= 0
    A[7,3] =  1.0
    A[7,4] = -1.0
    """
    # Now define lb, ub
    #lb = np.full(4, -np.inf)
    #ub = np.full(4, +np.inf)

    # b_i - b_{i+1} <= 0 => row 0..3 =>  -inf <= A[i,:].dot(p) <= 0
    #ub[0:4] = 0.0

    # a_i - a_{i+1} >= 0 => row 4..7 =>  0 <= A[i,:].dot(p) <= +inf
    #lb[4:8] = 0.0

    #linear_constraint = LinearConstraint(A, lb, ub)

    # Build the bounds
    bounds = [(-30, 30)]*5 + [(np.log(0.1), np.log(500.0))]*5
    
    # We pass the top-level cost function + args
    result = differential_evolution(
        cost_function_de_global,
        bounds=bounds,
        args=(E, data),         # <--- pass E, data here
        updating='deferred',
        workers=-1,        
        polish=True,
        maxiter= 20000,
        popsize=15, #constraints=linear_constraint,
        init='latinhypercube',
        strategy='best1bin'
    )
    
    return result.x, result.fun

# ------------------------------
# 4. Fit routine with differential_evolution
# ------------------------------
def fit_five_maxwellians2(E, data):
    """
    Fit sum of 5 Maxwellians using differential_evolution.
    -30 <= a_i <= 30
    ln(0.01) <= b_i <= ln(500)
    """
    
    """
    # We'll build a single A matrix of shape (8,10):
    # rows 0..3 => b_i - b_{i+1} <= 0
    # rows 4..7 => a_i - a_{i+1} >= 0

    A = np.zeros((4, 10))

    # Indices: a1->p[0], a2->p[1], a3->p[2], a4->p[3], a5->p[4]
    #          b1->p[5], b2->p[6], b3->p[7], b4->p[8], b5->p[9].

    # --- Enforce T1 <= T2 <= T3 <= T4 <= T5 ---
    # b1 - b2 <= 0
    A[0,5] =  1.0   # b1
    A[0,6] = -1.0   # b2

    # b2 - b3 <= 0
    A[1,6] =  1.0
    A[1,7] = -1.0

    # b3 - b4 <= 0
    A[2,7] =  1.0
    A[2,8] = -1.0

    # b4 - b5 <= 0
    A[3,8] =  1.0
    A[3,9] = -1.0

    # --- Enforce f1 >= f2 >= f3 >= f4 >= f5 ---
    # means a1 >= a2, a2 >= a3, etc.

    # a1 - a2 >= 0
    A[4,0] =  1.0
    A[4,1] = -1.0

    # a2 - a3 >= 0
    A[5,1] =  1.0
    A[5,2] = -1.0

    # a3 - a4 >= 0
    A[6,2] =  1.0
    A[6,3] = -1.0

    # a4 - a5 >= 0
    A[7,3] =  1.0
    A[7,4] = -1.0
    """
    # Now define lb, ub
    #lb = np.full(4, -np.inf)
    # ub = np.full(4, +np.inf)

    # b_i - b_{i+1} <= 0 => row 0..3 =>  -inf <= A[i,:].dot(p) <= 0
    # ub[0:4] = 0.0

    # a_i - a_{i+1} >= 0 => row 4..7 =>  0 <= A[i,:].dot(p) <= +inf
    #lb[4:8] = 0.0

    # linear_constraint = LinearConstraint(A, lb, ub)
    
    # Build the bounds
    bounds = [(-30, 30)]*5 + [(np.log(0.1), np.log(500.0))]*5
    
    # We pass the top-level cost function + args
    result = differential_evolution(
        cost_function_de_global,
        bounds=bounds,
        args=(E, data),         # <--- pass E, data here
        updating='deferred',
        workers=-1, 
        polish=True,
        maxiter= 20000,
        popsize=50, #constraints=linear_constraint,
        init='latinhypercube',
        strategy='best1bin'
    )
    
    return result.x, result.fun


def extract_f_T(params):
    a = params[:5]
    b = params[5:]
    exp_a = np.exp(a)
    sum_exp_a = np.sum(exp_a)
    f = exp_a / sum_exp_a
    T = np.exp(b)
    return f, T


def compute_percent_errors(E, data, params):
    model = five_maxwellians(E, params)
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
    # ------------------------------
    # Example usage or main loop
    # ------------------------------
    T_min = 1.0
    T_max = 100.
    num_points_T = 5
    T = np.logspace(np.log10(T_min), np.log10(T_max), num_points_T)
    
    kappa_min = 1.51
    kappa_max = 100.0
    num_points_kappa = 5
    kappa = np.logspace(np.log10(kappa_min), np.log10(kappa_max), num_points_kappa)

    E_min = 0.01
    E_max = 500.0
    num_points = 100
    E_vals = np.logspace(np.log10(E_min), np.log10(E_max), num_points)

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
        for j in range(len(kappa)):
            print(" i = ",i, " of 9, j = ",j, " of 9")

            data_standard  = np.array([standard_kappa(E, T[i], kappa[j]) for E in E_vals])
            data_product   = np.array([product_kappa(E,  T[i], kappa[j]) for E in E_vals])
            data_fried_egg = np.array([fried_egg(E,      T[i], kappa[j]) for E in E_vals])
    
            if np.min(data_standard) < 0.0:
                print("Warning Standard Kappa < 0, min =", np.min(data_standard))
            if np.min(data_product) < 0.0:
                print("Warning Product Kappa < 0, min =", np.min(data_product))
            if np.min(data_fried_egg) < 0.0:
                print("Warning Fried Egg < 0, min =", np.min(data_fried_egg))
    
            # Fit standard kappa
            params_std, cost_std = fit_five_maxwellians(E_vals, data_standard)
            f_std, T_std = extract_f_T(params_std)
            T_i_max4_sk[:, i, j] = T_std
            f_i_max4_sk[:, i, j] = f_std
    
            # Fit product kappa
            params_prod, cost_prod = fit_five_maxwellians(E_vals, data_product)
            f_prod, T_prod = extract_f_T(params_prod)
            T_i_max4_pk[:, i, j] = T_prod
            f_i_max4_pk[:, i, j] = f_prod
    
            # Fit fried-egg
            params_fried, cost_fried = fit_five_maxwellians(E_vals, data_fried_egg)
            f_fried, T_fried = extract_f_T(params_fried)
            T_i_max4_fe[:, i, j] = T_fried
            f_i_max4_fe[:, i, j] = f_fried
    
    
            # % error
            pct_err_std   = compute_percent_errors(E_vals, data_standard,  params_std)
            pct_err_prod  = compute_percent_errors(E_vals, data_product,   params_prod)
            pct_err_fried = compute_percent_errors(E_vals, data_fried_egg, params_fried)
    
            max_rel_err_sk[i, j]  = np.max(pct_err_std)
            max_rel_err_pk[i, j]  = np.max(pct_err_prod)
            max_rel_err_fe[i, j]  = np.max(pct_err_fried)
            

    
            if max_rel_err_sk[i, j] > 10.0:
                params_std, cost_std = fit_five_maxwellians2(E_vals, data_standard)
                f_std, T_std = extract_f_T(params_std)
                T_i_max4_sk[:, i, j] = T_std
                f_i_max4_sk[:, i, j] = f_std
                pct_err_std   = compute_percent_errors(E_vals, data_standard,  params_std)
                max_rel_err_sk[i, j]  = np.max(pct_err_std)
                if max_rel_err_sk[ i, j] > 10.0:
                    print("Warning Max Percent Error Standard-Kappa > 10%:", max_rel_err_sk[i, j])
            if max_rel_err_pk[i, j] > 10.0:
                # Fit product kappa
                params_prod, cost_prod = fit_five_maxwellians2(E_vals, data_product)
                f_prod, T_prod = extract_f_T(params_prod)
                T_i_max4_pk[:, i, j] = T_prod
                f_i_max4_pk[:, i, j] = f_prod
                pct_err_prod  = compute_percent_errors(E_vals, data_product,   params_prod)
                max_rel_err_pk[i, j]  = np.max(pct_err_prod)
                if max_rel_err_pk[i, j] > 10.0:
                    print("Warning Max Percent Error Product-Kappa > 10%:", max_rel_err_pk[i, j])
            if max_rel_err_fe[i, j] > 10.0:
                params_fried, cost_fried = fit_five_maxwellians2(E_vals, data_fried_egg)
                f_fried, T_fried = extract_f_T(params_fried)
                T_i_max4_fe[:, i, j] = T_fried
                f_i_max4_fe[:, i, j] = f_fried
                pct_err_fried = compute_percent_errors(E_vals, data_fried_egg, params_fried)
                max_rel_err_fe[i, j]  = np.max(pct_err_fried)
                if max_rel_err_fe[i, j] >10.0:
                    print("Warning Max Percent Error Fried-Egg > 10%:", max_rel_err_fe[i, j])
            
            for k in range(4):
                if T_i_max4_sk[k, i, j] > T_i_max4_sk[k + 1, i, j]:
                    print(f"Standard kappa Temp Constraints not Met for T {k}")
                if T_i_max4_pk[k, i, j] > T_i_max4_pk[k + 1, i, j]:
                    print(f"Product kappa Temp Constraints not Met for T {k}")
                if T_i_max4_fe[k, i, j] > T_i_max4_fe[k + 1, i, j]:
                    print(f"Fried-Egg Temp Constraints not Met for T {k}")


    # Save

    
    # 6) Save everything in a single .npz
    np.savez(
        "5x5_1-100_1.51-100_T_kappa_indiv_dif_ev_fully_contrained_fits_of_Evals.npz",
        T=T,
        kappa=kappa,  
        E_vals= E_vals,
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


    
    
    # Suppose we have arrays:
    # T: shape (nT,) of temperature values
    # kappa: shape (nkappa,) of kappa values
    # max_rel_err_sk, max_rel_err_pk, max_rel_err_fe: shape (nT, nkappa)
    # T_i_max4_sk, f_i_max4_sk, etc.: shape (5, nT, nkappa)
    #
    # We'll loop over j in range(len(kappa)) and produce subplots.
    
    for j in range(len(kappa)):
        # ------------------------------------------------------------
        # 1) Figure for T_i and f_i in 3 rows x 2 columns
        # ------------------------------------------------------------
        fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(12, 10), sharex=True)
    
        # ============ Row 0: Standard Kappa ============
        # LEFT subplot: T_i (5 curves)
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
    
        # RIGHT subplot: f_i (5 curves)
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
        # LEFT subplot: T_i
        axes[1,0].plot(T, T_i_max4_pk[0, :, j], label=r'T$_1$')
        axes[1,0].plot(T, T_i_max4_pk[1, :, j], label=r'T$_2$')
        axes[1,0].plot(T, T_i_max4_pk[2, :, j], label=r'T$_3$')
        axes[1,0].plot(T, T_i_max4_pk[3, :, j], label=r'T$_4$')
        axes[1,0].plot(T, T_i_max4_pk[4, :, j], label=r'T$_5$')
        axes[1,0].set_yscale('log')
        axes[1,0].set_ylim(0.1, 500)
        axes[1,0].set_ylabel(r"T$_i$ (eV)")
        axes[1,0].set_title(f"Product Kappa T$_i$, kappa={kappa[j]:.2f}")
        axes[1,0].legend()
    
        # RIGHT subplot: f_i
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
        # LEFT subplot: T_i
        axes[2,0].plot(T, T_i_max4_fe[0, :, j], label=r'T$_1$')
        axes[2,0].plot(T, T_i_max4_fe[1, :, j], label=r'T$_2$')
        axes[2,0].plot(T, T_i_max4_fe[2, :, j], label=r'T$_3$')
        axes[2,0].plot(T, T_i_max4_fe[3, :, j], label=r'T$_4$')
        axes[2,0].plot(T, T_i_max4_fe[4, :, j], label=r'T$_5$')
        axes[2,0].set_yscale('log')
        axes[2,0].set_ylim(0.1, 500)
        axes[2,0].set_ylabel(r"T$_i$ (eV)")
        axes[2,0].set_xlabel(r"T (eV)")
        axes[2,0].set_title(f"Fried-Egg T$_i$, kappa={kappa[j]:.2f}")
        axes[2,0].legend()
    
        # RIGHT subplot: f_i
        axes[2,1].plot(T, f_i_max4_fe[0, :, j], label=r'f$_1$')
        axes[2,1].plot(T, f_i_max4_fe[1, :, j], label=r'f$_2$')
        axes[2,1].plot(T, f_i_max4_fe[2, :, j], label=r'f$_3$')
        axes[2,1].plot(T, f_i_max4_fe[3, :, j], label=r'f$_4$')
        axes[2,1].plot(T, f_i_max4_fe[4, :, j], label=r'f$_5$')
        axes[2,1].set_yscale('log')
        axes[2,1].set_ylim(0.001, 1.0)
        axes[2,1].set_ylabel(r"f$_i$")
        axes[2,1].set_xlabel(r"T (eV)")
        axes[2,1].set_title(f"Fried-Egg f$_i$, kappa={kappa[j]:.2f}")
        axes[2,1].legend()
    
        plt.tight_layout()
        plt.savefig(f"5x5_1-100_1.51-100_T_kappa_Ti_fi_subplots_5MaxFits_kappa={kappa[j]:.2f}.pdf", bbox_inches='tight')
        plt.show()
    
        # ------------------------------------------------------------
        # 2) Figure for max_abs_percent_error in 3 rows × 1 column
        # ------------------------------------------------------------
        fig2, axes2 = plt.subplots(nrows=3, ncols=1, figsize=(8, 10), sharex=True)
    
        # Standard Kappa
        axes2[0].plot(T, max_rel_err_sk[:, j], label='Standard-Kappa')
        axes2[0].set_yscale('log')
        axes2[0].set_ylim(0.1, 100)
        axes2[0].set_ylabel("Max % Rel. Error")
        axes2[0].set_title(f"Max % Error vs. T, Standard Kappa, kappa={kappa[j]:.2f}")
        axes2[0].legend()
    
        # Product Kappa
        axes2[1].plot(T, max_rel_err_pk[:, j], label='Product-Kappa')
        axes2[1].set_yscale('log')
        axes2[1].set_ylim(0.1, 100)
        axes2[1].set_ylabel("Max % Rel. Error")
        axes2[1].set_title(f"Max % Error vs. T, Product Kappa, kappa={kappa[j]:.2f}")
        axes2[1].legend()
    
        # Fried-Egg
        axes2[2].plot(T, max_rel_err_fe[:, j], label='Fried-Egg')
        axes2[2].set_yscale('log')
        axes2[2].set_ylim(0.1, 100)
        axes2[2].set_xlabel(r"T (eV)")
        axes2[2].set_ylabel("Max % Rel. Error")
        axes2[2].set_title(f"Max % Error vs. T, Fried-Egg, kappa={kappa[j]:.2f}")
        axes2[2].legend()
    
        plt.tight_layout()
        plt.savefig(f"5x5_1-100_1.51-100_T_kappa_max_rel_error_subplots_5MaxFits_kappa={kappa[j]:.2f}.pdf", bbox_inches='tight')
        plt.show()

  
   # def rollinsg_median(array, window_size=30):
    #    """Returns the running median for a given 1D array using a specified window size."""
    #    series = pd.Series(array)
     #   return series.rolling(window_size, center=True, min_periods=1).median().values
    """
    npoints=50
    # Compute 30-point running medians
    median_sk = rolling_median(max_rel_err_sk, npoints)
    median_pk = rolling_median(max_rel_err_pk, npoints)
    median_fe = rolling_median(max_rel_err_fe, npoints)
    
    val=0.2
    
    # Create new arrays, keeping only those values within ±some% of the running median
    filtered_sk = np.where(
        (max_rel_err_sk >= val * median_sk) & (max_rel_err_sk <= (1. + val) * median_sk),
        max_rel_err_sk,
        np.nan
    )
    
    filtered_pk = np.where(
        (max_rel_err_pk >= val * median_pk) & (max_rel_err_pk <= (1. + val) * median_pk),
        max_rel_err_pk,
        np.nan
    )
    
    filtered_fe = np.where(
        (max_rel_err_fe >= val * median_fe) & (max_rel_err_fe <= (1. + val) * median_fe),
        max_rel_err_fe,
        np.nan
    )
    
    # Quick summary plot
    plt.figure(figsize=(9,6))
    plt.plot(r, filtered_sk,  label='Standard-Kappa')
    plt.plot(r, filtered_pk,  label='Product-Kappa')
    plt.plot(r, filtered_fe,  label='Fried-Egg')
    plt.yscale('log')
    plt.ylim(0.1, 100)
    plt.xlabel(r"$\rho_{ceq}$")
    plt.ylabel("Max % Rel. Error")
    plt.title("Filtered Max % Relative Error of 5 Maxwellian Fits (0.01 - 500 eV)")
   
    plt.legend()
    plt.tight_layout()
    plt.savefig("filtered_by_median_within_50percent_max_rel_error_of_fits_maxwellian_elec_for_nominal_model_fits_vs_rhoc_10x10_T_vs_kappa.pdf",
                bbox_inches='tight')
    plt.show()
    
    """
