#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 24 15:41:56 2024

@author: edne8319
"""

# -*- coding: utf-8 -*-
"""
Created on Sun Sep 29 16:59:39 2024

Modified to calculate kappa and T using analytical expressions derived from equating the 0th, 2nd, and 4th moments.

@author: Owner
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# Load data files (assuming they are in the same directory)
nop0 = np.loadtxt('nop_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv', delimiter=',')
no2p0 = np.loadtxt('no2p_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv', delimiter=',')
noph0 = np.loadtxt('noph_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv', delimiter=',')
nsp0 = np.loadtxt('nsp_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv', delimiter=',')
ns2p0 = np.loadtxt('ns2p_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv', delimiter=',')
ns3p0 = np.loadtxt('ns3p_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv', delimiter=',')
nhp0 = np.loadtxt('nhp_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv', delimiter=',')
nnap0 = np.loadtxt('nnap_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv', delimiter=',')
nec0 = np.loadtxt('new_nominal_model_nec_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv', delimiter=',')
neh0 = np.loadtxt('new_nominal_model_neh_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv', delimiter=',')

ti0 = np.loadtxt('Tic_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv', delimiter=',')
thp0 = np.loadtxt('Thp_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv', delimiter=',')
toph0 = np.loadtxt('Toph_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv', delimiter=',')
tec0 = np.loadtxt('Tec_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv', delimiter=',')
teh0 = np.loadtxt('new_nominal_model_Teh_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv', delimiter=',')

# Ensure all data arrays have the same length
num_points = len(tec0)

# Initialize arrays to store computed parameters
Te_electrons = np.zeros(num_points)
kappa_electrons = np.zeros(num_points)
Te_O = np.zeros(num_points)
kappa_O = np.zeros(num_points)
Te_S = np.zeros(num_points)
kappa_S = np.zeros(num_points)

# Radial positions from 4 to 10 RJ
R = np.linspace(4, 10, num_points)  # Radial distance array

# Loop over all radial positions
for i in range(num_points):
    # Electrons
    n_ec_cm3 = nec0[i]             # Core electron density in cm^-3
    n_eh_cm3 = neh0[i]             # Hot electron density in cm^-3
    T_ec_eV = tec0[i]              # Core electron temperature in eV
    T_eh_eV = teh0[i]              # Hot electron temperature in eV

    n_e = n_ec_cm3 + n_eh_cm3      # Total electron density in cm^-3

    if n_e > 0:
        f_c = n_ec_cm3 / n_e
        f_h = n_eh_cm3 / n_e

        # Calculate T and kappa for electrons
        numerator_T = -2 * (f_c * T_ec_eV + f_h * T_eh_eV) * (f_c * T_ec_eV**2 + f_h * T_eh_eV**2)
        denominator_T = (6 * f_c * f_h * T_ec_eV * T_eh_eV
                         + f_c * (3 * f_c - 5) * T_ec_eV**2
                         + f_h * (3 * f_h - 5) * T_eh_eV**2)
        if denominator_T != 0:
            T_e_eV = numerator_T / denominator_T
        else:
            T_e_eV = np.nan

        numerator_kappa = (6 * f_c * f_h * T_ec_eV * T_eh_eV
                           + f_c * (3 * f_c - 5) * T_ec_eV**2
                           + f_h * (3 * f_h - 5) * T_eh_eV**2)
        denominator_kappa = (4 * f_c * f_h * T_ec_eV * T_eh_eV
                             + 2 * (f_c - 1) * f_c * T_ec_eV**2
                             + 2 * (f_h - 1) * f_h * T_eh_eV**2)
        if denominator_kappa != 0:
            kappa_e = numerator_kappa / denominator_kappa
        else:
            kappa_e = np.nan

        # Check if T_e_eV and kappa_e are valid
        if T_e_eV > 0 and kappa_e > 1.5 and np.isfinite(T_e_eV) and np.isfinite(kappa_e):
            Te_electrons[i] = T_e_eV
            kappa_electrons[i] = kappa_e
        else:
            Te_electrons[i] = np.nan
            kappa_electrons[i] = np.nan
    else:
        Te_electrons[i] = np.nan
        kappa_electrons[i] = np.nan

    # Calculate hot densities for O⁺ and S²⁺
    nop = nop0[i]
    ns2p = ns2p0[i]
    noph = noph0[i]

    if ns2p > 0.0:
        ratio = nop / ns2p
        nophot = (noph * ratio) / (ratio + 1.0)      # Hot O⁺ density
        ns2phot = noph / (ratio + 1.0)               # Hot S²⁺ density
    else:
        nophot = noph      # Hot O⁺ density
        ns2phot = 0.0      # Hot S²⁺ density

    # For O⁺
    n_oc_cm3 = nop            # Cold O⁺ density in cm^-3
    n_oh_cm3 = nophot         # Hot O⁺ density in cm^-3
    n_O_cm3 = n_oc_cm3 + n_oh_cm3  # Total O⁺ density in cm^-3

    n_O = n_O_cm3

    if n_O > 0:
        f_c = n_oc_cm3 / n_O
        f_h = n_oh_cm3 / n_O

        T_oc_eV = ti0[i]           # Cold ion temperature in eV
        T_oh_eV = toph0[i]         # Hot ion temperature in eV

        # Calculate T and kappa for O⁺ ions
        numerator_T = -2 * (f_c * T_oc_eV + f_h * T_oh_eV) * (f_c * T_oc_eV**2 + f_h * T_oh_eV**2)
        denominator_T = (6 * f_c * f_h * T_oc_eV * T_oh_eV
                         + f_c * (3 * f_c - 5) * T_oc_eV**2
                         + f_h * (3 * f_h - 5) * T_oh_eV**2)
        if denominator_T != 0:
            T_O_eV = numerator_T / denominator_T
        else:
            T_O_eV = np.nan

        numerator_kappa = (6 * f_c * f_h * T_oc_eV * T_oh_eV
                           + f_c * (3 * f_c - 5) * T_oc_eV**2
                           + f_h * (3 * f_h - 5) * T_oh_eV**2)
        denominator_kappa = (4 * f_c * f_h * T_oc_eV * T_oh_eV
                             + 2 * (f_c - 1) * f_c * T_oc_eV**2
                             + 2 * (f_h - 1) * f_h * T_oh_eV**2)
        if denominator_kappa != 0:
            kappa_O_ion = numerator_kappa / denominator_kappa
        else:
            kappa_O_ion = np.nan

        # Check if T_O_eV and kappa_O_ion are valid
        if T_O_eV > 0 and kappa_O_ion > 1.5 and np.isfinite(T_O_eV) and np.isfinite(kappa_O_ion):
            Te_O[i] = T_O_eV
            kappa_O[i] = kappa_O_ion
        else:
            Te_O[i] = np.nan
            kappa_O[i] = np.nan
    else:
        Te_O[i] = np.nan
        kappa_O[i] = np.nan

    # For S²⁺ ions
    n_s2c_cm3 = ns2p0[i]           # Cold S²⁺ density in cm^-3
    n_s2h_cm3 = ns2phot            # Hot S²⁺ density in cm^-3
    n_S2_cm3 = n_s2c_cm3 + n_s2h_cm3

    n_S2 = n_S2_cm3

    if n_S2 > 0:
        f_c = n_s2c_cm3 / n_S2
        f_h = n_s2h_cm3 / n_S2

        T_s2c_eV = ti0[i]              # Cold ion temperature in eV
        T_s2h_eV = toph0[i]            # Hot ion temperature in eV

        # Calculate T and kappa for S²⁺ ions
        numerator_T = -2 * (f_c * T_s2c_eV + f_h * T_s2h_eV) * (f_c * T_s2c_eV**2 + f_h * T_s2h_eV**2)
        denominator_T = (6 * f_c * f_h * T_s2c_eV * T_s2h_eV
                         + f_c * (3 * f_c - 5) * T_s2c_eV**2
                         + f_h * (3 * f_h - 5) * T_s2h_eV**2)
        if denominator_T != 0:
            T_S2_eV = numerator_T / denominator_T
        else:
            T_S2_eV = np.nan

        numerator_kappa = (6 * f_c * f_h * T_s2c_eV * T_s2h_eV
                           + f_c * (3 * f_c - 5) * T_s2c_eV**2
                           + f_h * (3 * f_h - 5) * T_s2h_eV**2)
        denominator_kappa = (4 * f_c * f_h * T_s2c_eV * T_s2h_eV
                             + 2 * (f_c - 1) * f_c * T_s2c_eV**2
                             + 2 * (f_h - 1) * f_h * T_s2h_eV**2)
        if denominator_kappa != 0:
            kappa_S2 = numerator_kappa / denominator_kappa
        else:
            kappa_S2 = np.nan

        # Check if T_S2_eV and kappa_S2 are valid
        if T_S2_eV > 0 and kappa_S2 > 1.5 and np.isfinite(T_S2_eV) and np.isfinite(kappa_S2):
            Te_S[i] = T_S2_eV
            kappa_S[i] = kappa_S2
        else:
            Te_S[i] = np.nan
            kappa_S[i] = np.nan
    else:
        Te_S[i] = np.nan
        kappa_S[i] = np.nan

# Handle NaN values and interpolate where necessary
# For simplicity, we will perform linear interpolation over finite values
def interpolate_finite(R, values):
    finite_indices = np.isfinite(values)
    if np.sum(finite_indices) > 1:
        interp_func = interp1d(R[finite_indices], values[finite_indices], bounds_error=False, fill_value="extrapolate")
        return interp_func(R)
    else:
        return values

Te_electrons = interpolate_finite(R, Te_electrons)
kappa_electrons = interpolate_finite(R, kappa_electrons)
Te_O = interpolate_finite(R, Te_O)
kappa_O = interpolate_finite(R, kappa_O)
Te_S = interpolate_finite(R, Te_S)
kappa_S = interpolate_finite(R, kappa_S)

# Save the computed parameters to CSV files
np.savetxt('new_nominal_modelTe_electrons_024_analytic_momentsequal_bimax_to_kappa.csv', Te_electrons, delimiter=',')
np.savetxt('new_nominal_modelappa_electrons_024_analytic_momentsequal_bimax_to_kappa.csv', kappa_electrons, delimiter=',')
np.savetxt('Te_O_024_analytic_momentsequal_bimax_to_kappa.csv', Te_O, delimiter=',')
np.savetxt('kappa_O_024_analytic_momentsequal_bimax_to_kappa.csv', kappa_O, delimiter=',')
np.savetxt('Te_S_024_analytic_momentsequal_bimax_to_kappa.csv', Te_S, delimiter=',')
np.savetxt('kappa_S_024_analytic_momentsequal_bimax_to_kappa.csv', kappa_S, delimiter=',')

# Plotting
R = np.linspace(4.0, 10.0, num_points)
plt.figure(figsize=(10, 6))
plt.semilogy(R, tec0, label='Electron Cold Temperature', linestyle='--')
plt.semilogy(R, Te_electrons, label='Electron Kappa Temperature')
plt.semilogy(R, ti0, label='O⁺ & S²⁺ Cold Temperature', linestyle='--')
plt.semilogy(R, Te_O, label='O⁺ Kappa Temperature')
plt.semilogy(R, Te_S, label='S²⁺ Kappa Temperature')
plt.xlabel(r'Radial Distance ($R_J$)', fontsize=14)
plt.ylabel('Temperature (eV)', fontsize=14)
plt.title('Kappa Temperatures vs Radial Distance', fontsize=16)
plt.legend(fontsize=12)
plt.grid(True)
plt.tight_layout()
plt.savefig('new_nominal_model_024_analytic_kappa_temperatures_vs_radial_distance_momentsequal_bimax_to_kappa.png', dpi=300)
plt.show()

# If you have Steffl (2004b) data, load and plot it
try:
    r_vs_kappa_steffl2004b = np.loadtxt('r_vs_kappa_steffl2004b.csv', delimiter=',', skiprows=1)
    r_vs_kappa_steffl2004b_r = r_vs_kappa_steffl2004b[:, 0]
    r_vs_kappa_steffl2004b_kappa = r_vs_kappa_steffl2004b[:, 1]

    plt.figure(figsize=(10, 6))
    plt.semilogy(R, kappa_electrons, label='Electron κ')
    plt.semilogy(R, kappa_O, label='O⁺ κ')
    plt.semilogy(R, kappa_S, label='S²⁺ κ')
    plt.semilogy(r_vs_kappa_steffl2004b_r, r_vs_kappa_steffl2004b_kappa, label='Steffl (2004b) Electron κ')
    plt.xlabel('Radial Distance (Rj)', fontsize=14)
    plt.ylabel('Kappa Value (κ)', fontsize=14)
    plt.title('Kappa Values vs Radial Distance', fontsize=16)
    plt.ylim((1e0, 1e1))
    plt.legend(fontsize=12)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig('new_nominal_model_024_analytic_kappa_values_vs_radial_distance_momentsequal_bimax_to_kappa.png', dpi=300)
    plt.show()
except Exception as e:
    print("Could not load Steffl (2004b) data:", e)
    # Plot without Steffl data
    plt.figure(figsize=(10, 6))
    plt.semilogy(R, kappa_electrons, label='Electron κ')
    plt.semilogy(R, kappa_O, label='O⁺ κ')
    plt.semilogy(R, kappa_S, label='S²⁺ κ')
    plt.xlabel('Radial Distance (Rj)', fontsize=14)
    plt.ylabel('Kappa Value (κ)', fontsize=14)
    plt.title('Kappa Values vs Radial Distance', fontsize=16)
    plt.legend(fontsize=12)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig('new_nominal_model_024_analytic_kappa_values_vs_radial_distance_012_momentsequal_bimax_to_kappa.png', dpi=300)
    plt.show()

print('Computation completed and results saved to CSV files.')
