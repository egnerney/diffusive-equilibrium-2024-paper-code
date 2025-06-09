# -*- coding: utf-8 -*-
"""
Created on Sun Apr 13 23:33:28 2025

@author: Owner

Filter & sort (T_i, f_i) for 3 distributions (Standard Kappa, Product Kappa,
Fried-Egg). Then do a second-level filter based on a rolling median of the
(filtered) max_rel_error. Finally, plot the results.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# -------------------------------------------------------------------------
# 0) Read radius array (assuming 4..10 RJ, 601 points)
# -------------------------------------------------------------------------
r = np.linspace(4.0, 10.0, 601)

# -------------------------------------------------------------------------
# 1) Read the fitted 5x601 CSV files for each distribution
# -------------------------------------------------------------------------
f_i_sk = np.loadtxt('f_i_fit_5Max_to_Standard-Kappa_model_elec_for_emiss_calc_5x601.csv', delimiter=',')
T_i_sk = np.loadtxt('T_i_fit_5Max_to_Standard-Kappa_model_elec_for_emiss_calc_5x601.csv', delimiter=',')

f_i_pk = np.loadtxt('f_i_fit_5Max_to_Product-Kappa_model_elec_for_emiss_calc_5x601.csv', delimiter=',')
T_i_pk = np.loadtxt('T_i_fit_5Max_to_Product-Kappa_model_elec_for_emiss_calc_5x601.csv', delimiter=',')

f_i_fe = np.loadtxt('f_i_fit_5Max_to_Fried-Egg_model_elec_for_emiss_calc_5x601.csv', delimiter=',')
T_i_fe = np.loadtxt('T_i_fit_5Max_to_Fried-Egg_model_elec_for_emiss_calc_5x601.csv', delimiter=',')

# -------------------------------------------------------------------------
# 2) Read the max_rel_err arrays for each distribution
# -------------------------------------------------------------------------
max_rel_err_sk = np.loadtxt('max_rel_err_sk.csv', delimiter=',')  # shape (601,)
max_rel_err_pk = np.loadtxt('max_rel_err_pk.csv', delimiter=',')  
max_rel_err_fe = np.loadtxt('max_rel_err_fe.csv', delimiter=',')

# -------------------------------------------------------------------------
# 3) First-level filter: keep only columns with max_rel_err < 10%.
#    Then sort T_i ascending in those columns, reorder f_i.
# -------------------------------------------------------------------------
ERROR_CUTOFF_1 = 10.0  # first-level filter

f_i_sk_sorted1 = np.full_like(f_i_sk, np.nan)
T_i_sk_sorted1 = np.full_like(T_i_sk, np.nan)

f_i_pk_sorted1 = np.full_like(f_i_pk, np.nan)
T_i_pk_sorted1 = np.full_like(T_i_pk, np.nan)

f_i_fe_sorted1 = np.full_like(f_i_fe, np.nan)
T_i_fe_sorted1 = np.full_like(T_i_fe, np.nan)

ncols = f_i_sk.shape[1]

for i in range(ncols):
    # Standard Kappa
    if max_rel_err_sk[i] < ERROR_CUTOFF_1:
        sort_idx = np.argsort(T_i_sk[:, i])
        T_i_sk_sorted1[:, i] = T_i_sk[sort_idx, i]
        f_i_sk_sorted1[:, i] = f_i_sk[sort_idx, i]

    # Product Kappa
    if max_rel_err_pk[i] < ERROR_CUTOFF_1:
        sort_idx = np.argsort(T_i_pk[:, i])
        T_i_pk_sorted1[:, i] = T_i_pk[sort_idx, i]
        f_i_pk_sorted1[:, i] = f_i_pk[sort_idx, i]

    # Fried-Egg
    if max_rel_err_fe[i] < ERROR_CUTOFF_1:
        sort_idx = np.argsort(T_i_fe[:, i])
        T_i_fe_sorted1[:, i] = T_i_fe[sort_idx, i]
        f_i_fe_sorted1[:, i] = f_i_fe[sort_idx, i]

# We'll also create a "filtered" max_rel_err that is NaN for columns >= 10%
filtered_max_rel_err_sk = np.where(max_rel_err_sk < ERROR_CUTOFF_1, max_rel_err_sk, np.nan)
filtered_max_rel_err_pk = np.where(max_rel_err_pk < ERROR_CUTOFF_1, max_rel_err_pk, np.nan)
filtered_max_rel_err_fe = np.where(max_rel_err_fe < ERROR_CUTOFF_1, max_rel_err_fe, np.nan)

# -------------------------------------------------------------------------
# 4) Second-level filter: 
#    - Compute rolling median (window=30) on the "filtered" max_err arrays
#    - Exclude columns more than 30% above that running median
# -------------------------------------------------------------------------
def rolling_median(arr, window_size=30):
    """
    Returns the running median for a given 1D array using a specified window size.
    We use pandas for convenience. We'll do center=True and min_periods=1 to
    get a full-length array of medians.
    """
    s = pd.Series(arr)
    return s.rolling(window_size, center=True, min_periods=1).median().values

# We'll compute separate rolling medians for each distribution
rm_sk = rolling_median(filtered_max_rel_err_sk, 30)
rm_pk = rolling_median(filtered_max_rel_err_pk, 30)
rm_fe = rolling_median(filtered_max_rel_err_fe, 30)

# Next define the second-level filter.  
# We'll say: keep the column if (1) it wasn't already NaN from the first filter,
# and (2) it is <= 1.3 * rolling_median.
ERROR_CUTOFF_2_FACTOR = 1.3

final_mask_sk = ( ~np.isnan(filtered_max_rel_err_sk) &
                  (filtered_max_rel_err_sk <= rm_sk * ERROR_CUTOFF_2_FACTOR) )

final_mask_pk = ( ~np.isnan(filtered_max_rel_err_pk) &
                  (filtered_max_rel_err_pk <= rm_pk * ERROR_CUTOFF_2_FACTOR) )

final_mask_fe = ( ~np.isnan(filtered_max_rel_err_fe) &
                  (filtered_max_rel_err_fe <= rm_fe * ERROR_CUTOFF_2_FACTOR) )

# -------------------------------------------------------------------------
# 5) Create final T_i, f_i arrays after second-level filtering
#    We'll keep shape (5, 601) but set failing columns to NaN.
# -------------------------------------------------------------------------
f_i_sk_sorted2 = np.full_like(f_i_sk_sorted1, np.nan)
T_i_sk_sorted2 = np.full_like(T_i_sk_sorted1, np.nan)

f_i_pk_sorted2 = np.full_like(f_i_pk_sorted1, np.nan)
T_i_pk_sorted2 = np.full_like(T_i_pk_sorted1, np.nan)

f_i_fe_sorted2 = np.full_like(f_i_fe_sorted1, np.nan)
T_i_fe_sorted2 = np.full_like(T_i_fe_sorted1, np.nan)

for i in range(ncols):
    # Standard Kappa
    if final_mask_sk[i]:
        f_i_sk_sorted2[:, i] = f_i_sk_sorted1[:, i]
        T_i_sk_sorted2[:, i] = T_i_sk_sorted1[:, i]

    # Product Kappa
    if final_mask_pk[i]:
        f_i_pk_sorted2[:, i] = f_i_pk_sorted1[:, i]
        T_i_pk_sorted2[:, i] = T_i_pk_sorted1[:, i]

    # Fried-Egg
    if final_mask_fe[i]:
        f_i_fe_sorted2[:, i] = f_i_fe_sorted1[:, i]
        T_i_fe_sorted2[:, i] = T_i_fe_sorted1[:, i]

# And similarly for the final max_rel_err (second-level filtered)
final_max_rel_err_sk = np.where(final_mask_sk, filtered_max_rel_err_sk, np.nan)
final_max_rel_err_pk = np.where(final_mask_pk, filtered_max_rel_err_pk, np.nan)
final_max_rel_err_fe = np.where(final_mask_fe, filtered_max_rel_err_fe, np.nan)


# -------------------------------------------------------------------------
# 6) Plot the final filtered max_rel_err vs r
# -------------------------------------------------------------------------
plt.figure(figsize=(8,5))
plt.plot(r, final_max_rel_err_sk, label='Standard-Kappa', alpha=0.7)
plt.plot(r, final_max_rel_err_pk, label='Product-Kappa',  alpha=0.7)
plt.plot(r, final_max_rel_err_fe, label='Fried-Egg',      alpha=0.7)
plt.title("Final Filtered Max % Rel. Error vs. r (two-level filter: <10% + <=1.3*running_median)")
plt.xlabel("r [R_J]")
plt.ylabel("Max % Error")
plt.yscale('log')  # optional
plt.ylim(0.1, 10)
plt.legend()
plt.tight_layout()
plt.show()

# -------------------------------------------------------------------------
# 7) Plot T_i, f_i after second-level filtering
#    We'll do subplots: T_i in 3 rows, f_i in 3 rows
# -------------------------------------------------------------------------
fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(12, 10), sharex=True)
# axes[i, 0] => T_i plot for distribution i
# axes[i, 1] => f_i plot for distribution i

# ---------- Standard Kappa ----------
axes[0,0].set_title("T_i (sorted) vs r - Std Kappa - final filter")
for row in range(5):
    axes[0,0].plot(r, T_i_sk_sorted2[row,:], label=f"T{row+1}")
axes[0,0].set_ylabel("Temperature (eV)")
axes[0,0].set_yscale('log')
axes[0,0].legend()

axes[0,1].set_title("f_i (sorted) vs r - Std Kappa - final filter")
for row in range(5):
    axes[0,1].plot(r, f_i_sk_sorted2[row,:], label=f"f{row+1}")
axes[0,1].set_ylabel("fraction")
axes[0,1].legend()

# ---------- Product Kappa ----------
axes[1,0].set_title("T_i (sorted) vs r - Prod Kappa - final filter")
for row in range(5):
    axes[1,0].plot(r, T_i_pk_sorted2[row,:], label=f"T{row+1}")
axes[1,0].set_ylabel("Temperature (eV)")
axes[1,0].set_yscale('log')
axes[1,0].legend()

axes[1,1].set_title("f_i (sorted) vs r - Prod Kappa - final filter")
for row in range(5):
    axes[1,1].plot(r, f_i_pk_sorted2[row,:], label=f"f{row+1}")
axes[1,1].set_ylabel("fraction")
axes[1,1].legend()

# ---------- Fried-Egg ----------
axes[2,0].set_title("T_i (sorted) vs r - Fried-Egg - final filter")
for row in range(5):
    axes[2,0].plot(r, T_i_fe_sorted2[row,:], label=f"T{row+1}")
axes[2,0].set_ylabel("Temperature (eV)")
axes[2,0].set_yscale('log')
axes[2,0].set_xlabel("r [R_J]")
axes[2,0].legend()

axes[2,1].set_title("f_i (sorted) vs r - Fried-Egg - final filter")
for row in range(5):
    axes[2,1].plot(r, f_i_fe_sorted2[row,:], label=f"f{row+1}")
axes[2,1].set_ylabel("fraction")
axes[2,1].set_xlabel("r [R_J]")
axes[2,1].legend()

plt.tight_layout()
plt.savefig("plots.pdf")
plt.show()

# -------------------------------------------------------------------------
# 8) (Optional) Save final arrays to CSV
# -------------------------------------------------------------------------
np.savetxt('f_i_fit_5Max_Standard-Kappa_filtered_sorted_2level.csv', f_i_sk_sorted2, delimiter=',')
np.savetxt('T_i_fit_5Max_Standard-Kappa_filtered_sorted_2level.csv', T_i_sk_sorted2, delimiter=',')
np.savetxt('f_i_fit_5Max_Product-Kappa_filtered_sorted_2level.csv',  f_i_pk_sorted2, delimiter=',')
np.savetxt('T_i_fit_5Max_Product-Kappa_filtered_sorted_2level.csv',  T_i_pk_sorted2, delimiter=',')
np.savetxt('f_i_fit_5Max_Fried-Egg_filtered_sorted_2level.csv',      f_i_fe_sorted2, delimiter=',')
np.savetxt('T_i_fit_5Max_Fried-Egg_filtered_sorted_2level.csv',      T_i_fe_sorted2, delimiter=',')

print("Done. Two-level filtering complete. Results saved and plotted.")
