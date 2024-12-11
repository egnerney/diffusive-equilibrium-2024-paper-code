# -*- coding: utf-8 -*-
"""
Created on Sun Sep 29 16:59:39 2024

Modified to numerically solve for kappa and T for each species using the product kappa distribution.

@author: Owner
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gamma
from scipy.optimize import fsolve
from scipy.interpolate import interp1d


guess_T_O = np.loadtxt('Te_O_012_momentsequal_bimax_to_kappa.csv', delimiter=',')
guess_kappa_O = np.loadtxt('kappa_O_012_momentsequal_bimax_to_kappa.csv', delimiter=',')


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

# Physical constants
m_e = 9.10938356e-31        # Electron mass in kg
e = 1.602176634e-19         # Elementary charge in C (J/eV)
pi = np.pi

# Ion masses in kg
m_O = 16. * 1.66053906660e-27        # Mass of O+ ion (16 u)
m_S = 32. * 1.66053906660e-27        # Mass of S++ ion (32 u)

# Arrays to store computed parameters
Te_electrons = np.zeros(num_points)
kappa_electrons = np.zeros(num_points)
Te_O = np.zeros(num_points)
kappa_O = np.zeros(num_points)
Te_S = np.zeros(num_points)
kappa_S = np.zeros(num_points)

# Radial positions from 4 to 10 RJ
R = np.linspace(4, 10, num_points)  # Radial distance array

alpha_ions = -0.92
beta_ions = -17.655

alpha_e = -7.221406  # -4.07075
beta_e = -21.4182    # -31.3772
gamma_e = -40.3772

# Loop over all radial positions
#for i in range(num_points):
    #Loop backwards now!
for i in range(num_points - 1, -1, -1):
    print(i)
    if R[i] > 6.:
        kappa_guess_electrons = 100.*((R[i]/6.) ** alpha_e)
        kappa_guess_s2p_and_Op = 4.*((R[i]/6.) ** alpha_ions)
    else:
        if R[i] > 5.8:
            kappa_guess_electrons = 100.*((R[i]/6.) ** beta_e)
            kappa_guess_s2p_and_Op = 4.*((R[i]/6.) ** beta_ions)
        else:
            if R[i] > 5.5:
                kappa_guess_electrons = 100.*((R[i]/6.) ** gamma_e)
                kappa_guess_s2p_and_Op = 4.*((R[i]/6.) ** beta_ions)
            else:
                kappa_guess_electrons = 500.
                kappa_guess_s2p_and_Op = 300.

    # Electrons
    n_ec_cm3 = nec0[i]             # Core electron density in cm^-3
    n_eh_cm3 = neh0[i]             # Hot electron density in cm^-3
    T_ec_eV = tec0[i]              # Core electron temperature in eV
    T_eh_eV = teh0[i]              # Hot electron temperature in eV

    # Convert densities to m^-3
    n_ec = n_ec_cm3 * 1e6       # Core electron density in m^-3
    n_eh = n_eh_cm3 * 1e6       # Hot electron density in m^-3
    n_e = n_ec + n_eh           # Total electron density in m^-3

    # Number density conservation
    n = n_e

    # Define T_mean and S for electrons
    T_mean = (n_ec * T_ec_eV + n_eh * T_eh_eV) / n
    S = (n_ec * np.sqrt(T_ec_eV) + n_eh * np.sqrt(T_eh_eV)) / n

    # Function to solve for electrons using the product kappa distribution
    def equations_electron(vars):
        T_e_eV, kappa_e = vars
        numerator = 2 * kappa_e * (3 * kappa_e - 2)
        denominator = 6 * kappa_e**2 - 9 * kappa_e + 3
        f1 = T_mean - T_e_eV * numerator / denominator
        f2 = S - np.sqrt(T_e_eV * kappa_e) * gamma(kappa_e - 0.5) / gamma(kappa_e)
        return [f1, f2]

    # Initial guesses for electrons
    if R[i] > 7:
        initial_guess = [T_mean, 4.0]
    else:
        initial_guess = [T_ec_eV, 4.0]

    # Solve equations numerically
    sol_electron, info, ier, mesg = fsolve(equations_electron, initial_guess, full_output=True)

    T_e_eV, kappa_e = sol_electron

    # Check if solution is valid
    if ier == 1 and kappa_e > 1.505 and T_e_eV > 0:
        Te_electrons[i] = T_e_eV
        kappa_electrons[i] = kappa_e
    else:
        print('electron error', kappa_e, T_e_eV, 'set to nans')
        print(ier, mesg)
        Te_electrons[i] = np.nan
        kappa_electrons[i] = np.nan

    # Now perform the same calculations for O⁺ ions
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

    n_oc = n_oc_cm3 * 1e6     # Convert to m^-3
    n_oh = n_oh_cm3 * 1e6
    n_O = n_oc + n_oh

    T_oc_eV = ti0[i]           # Cold ion temperature in eV
    T_oh_eV = toph0[i]         # Hot ion temperature in eV

    # Number density conservation
    n = n_O

    # Define T_mean and S for O⁺ ions
    T_mean_O = (n_oc * T_oc_eV + n_oh * T_oh_eV) / n
    S_O = (n_oc * np.sqrt(T_oc_eV) + n_oh * np.sqrt(T_oh_eV)) / n

    # Function to solve for O⁺ ions using the product kappa distribution
    def equations_O(vars):
        T_e_eV, kappa_e = vars
        numerator = 2 * kappa_e * (3 * kappa_e - 2)
        denominator = 6 * kappa_e**2 - 9 * kappa_e + 3
        f1 = T_mean_O - T_e_eV * numerator / denominator
        f2 = S_O - np.sqrt(T_e_eV * kappa_e) * gamma(kappa_e - 0.5) / gamma(kappa_e)
        return [f1, f2]

    # Initial guesses for O⁺ ions
    #initial_guess_O = [guess_T_O[i] , guess_kappa_O[i]/1.84] #[ T_oc_eV, 100.]#10.*((R[i]/6.) ** (-2.3))] #[T_mean_O, 5.0]

    # Initial guesses for O⁺ ions
    if R[i] > 5.9:
        initial_guess_O = [guess_T_O[i] , guess_kappa_O[i]/1.84]
    else:
        initial_guess_O = [Te_O[i + 1] , kappa_O[ i + 1]]
    # Solve equations numerically
    sol_O, info, ier, mesg = fsolve(equations_O, initial_guess_O, full_output=True)

    T_O_eV, kappa_O_ion = sol_O

    # Check if solution is valid
    if ier == 1 and kappa_O_ion > 1.505 and T_O_eV > 0:
        Te_O[i] = T_O_eV
        kappa_O[i] = kappa_O_ion
    else:
        print('O+ error', kappa_e, T_e_eV, 'set to nans')
        Te_O[i] = np.nan
        kappa_O[i] = np.nan

    # Now perform the same calculations for S²⁺ ions
    n_s2c_cm3 = ns2p0[i]           # Cold S²⁺ density in cm^-3
    n_s2h_cm3 = ns2phot            # Hot S²⁺ density in cm^-3
    n_S2_cm3 = n_s2c_cm3 + n_s2h_cm3

    n_s2c = n_s2c_cm3 * 1e6        # Convert to m^-3
    n_s2h = n_s2h_cm3 * 1e6
    n_S2 = n_s2c + n_s2h

    T_s2c_eV = ti0[i]              # Cold ion temperature in eV
    T_s2h_eV = toph0[i]            # Hot ion temperature in eV

    # Number density conservation
    n = n_S2

    # Define T_mean and S for S²⁺ ions
    T_mean_S2 = (n_s2c * T_s2c_eV + n_s2h * T_s2h_eV) / n if n > 0 else T_s2c_eV
    S_S2 = (n_s2c * np.sqrt(T_s2c_eV) + n_s2h * np.sqrt(T_s2h_eV)) / n if n > 0 else np.sqrt(T_s2c_eV)

    # Function to solve for S²⁺ ions using the product kappa distribution
    def equations_S2(vars):
        T_e_eV, kappa_e = vars
        numerator = 2 * kappa_e * (3 * kappa_e - 2)
        denominator = 6 * kappa_e**2 - 9 * kappa_e + 3
        f1 = T_mean_S2 - T_e_eV * numerator / denominator
        f2 = S_S2 - np.sqrt(T_e_eV * kappa_e) * gamma(kappa_e - 0.5) / gamma(kappa_e)
        return [f1, f2]

    # Initial guesses for S²⁺ ions
    if R[i] > 5.9:
        initial_guess_S2 = [guess_T_O[i] , guess_kappa_O[i]/1.84]
    else:
        initial_guess_S2 = [Te_S[i + 1] , kappa_S[ i + 1]]
        
     #[ T_s2c_eV, 100.]#10.*((R[i]/6.) ** (-2.3))]#[T_mean_S2, 5.0]

    # Solve equations numerically
    if n > 0:
        sol_S2, info, ier, mesg = fsolve(equations_S2, initial_guess_S2, full_output=True)
        T_S2_eV, kappa_S2 = sol_S2

        # Check if solution is valid
        if ier == 1 and kappa_S2 > 1.505 and T_S2_eV > 0:
            Te_S[i] = T_S2_eV
            kappa_S[i] = kappa_S2
        else:
            print('S2+ error', kappa_e, T_e_eV, 'set to nans')
            Te_S[i] = np.nan
            kappa_S[i] = np.nan
    else:
        # If n_S2 is zero, assign NaN
        print('S2+ error', kappa_e, T_e_eV, 'set to nans')
        Te_S[i] = np.nan
        kappa_S[i] = np.nan
    
    if R[i] < 5.67:
        Te_S[i] = guess_T_O[i]
        Te_O[i] = guess_T_O[i]
        kappa_S[i] = guess_kappa_O[i] / 1.6254312561173532  # ratio of standard kappa to product iso kappa at last usable solution poiint of 5.67 RJ
        kappa_O[i] = guess_kappa_O[i] / 1.6254312561173532 
        
    

# Handle NaN values and interpolate where necessary
for i in range(len(R)):
    if R[i] > 5.8:
        Te_electrons[i] = interp1d(R, Te_electrons)(R[i])
        kappa_electrons[i] = interp1d(R, kappa_electrons)(R[i])
        Te_O[i] = interp1d(R, Te_O)(R[i])
        kappa_O[i] = interp1d(R, kappa_O)(R[i])
        Te_S[i] = interp1d(R, Te_S)(R[i])
        kappa_S[i] = interp1d(R, kappa_S)(R[i])
    else:
        if R[i] > 5.7:
            Te_electrons[i] = interp1d(R, Te_electrons)(R[i])
            kappa_electrons[i] = interp1d(R, kappa_electrons)(R[i])
            Te_O[i] = interp1d(R, Te_O)(R[i])
            kappa_O[i] = interp1d(R, kappa_O)(R[i])
            Te_S[i] = interp1d(R, Te_S)(R[i])
            kappa_S[i] = interp1d(R, kappa_S)(R[i])
        else:
            if R[i] > 5.5:
                Te_electrons[i] = tec0[i]
                kappa_electrons[i] = 300.0
                Te_O[i] = interp1d(R, Te_O)(R[i])
                kappa_O[i] = interp1d(R, kappa_O)(R[i])
                Te_S[i] = interp1d(R, Te_S)(R[i])
                kappa_S[i] = interp1d(R, kappa_S)(R[i])
            else:
                if R[i] > 5.0:
                    Te_electrons[i] = tec0[i]
                    kappa_electrons[i] = 300.0
                    Te_O[i] = ti0[i]
                    vali = interp1d(R, kappa_O)(5.5)
                    pvali = np.log(300.0 / vali) / np.log(5.0 / 5.5)
                    kappa_O[i] = vali * ((R[i] / 5.5) ** pvali)
                    Te_S[i] = ti0[i]
                    kappa_S[i] = vali * ((R[i] / 5.5) ** pvali)
                else:
                    Te_electrons[i] = tec0[i]
                    kappa_electrons[i] = 300.0
                    Te_O[i] = ti0[i]
                    kappa_O[i] = 300.0
                    Te_S[i] = ti0[i]
                    kappa_S[i] = 300.0

# Interpolate finite values
Te_electrons_finite = Te_electrons[np.isfinite(Te_electrons)]
R_finite = R[np.isfinite(Te_electrons)]

Te_electrons = interp1d(R_finite, Te_electrons_finite, fill_value='extrapolate')(R)
kappa_electrons_finite = kappa_electrons[np.isfinite(kappa_electrons)]
R_finite = R[np.isfinite(kappa_electrons)]
kappa_electrons = interp1d(R_finite, kappa_electrons_finite, fill_value='extrapolate')(R)
Te_O = interp1d(R, Te_O, fill_value='extrapolate')(R)
kappa_O = interp1d(R, kappa_O, fill_value='extrapolate')(R)
Te_S = interp1d(R, Te_S, fill_value='extrapolate')(R)
kappa_S = interp1d(R, kappa_S, fill_value='extrapolate')(R)

# Save the computed parameters to CSV files
np.savetxt('new_nominal_model_Te_electrons_012_momentsequal_bimax_to_product_kappa.csv', Te_electrons, delimiter=',')
np.savetxt('new_nominal_model_kappa_electrons_012_momentsequal_bimax_to_product_kappa.csv', kappa_electrons, delimiter=',')
np.savetxt('Te_O_012_momentsequal_bimax_to_product_kappa.csv', Te_O, delimiter=',')
np.savetxt('kappa_O_012_momentsequal_bimax_to_product_kappa.csv', kappa_O, delimiter=',')
np.savetxt('Te_S_012_momentsequal_bimax_to_product_kappa.csv', Te_S, delimiter=',')
np.savetxt('kappa_S_012_momentsequal_bimax_to_product_kappa.csv', kappa_S, delimiter=',')

# Plotting
R = np.linspace(4.0, 10.0, 601)
# Set default parameters for bold text and thicker tick marks
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.titleweight'] = 'bold'
plt.rcParams['axes.linewidth'] = 2
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['ytick.major.width'] = 2

# First plot
plt.figure(figsize=(10, 6))
plt.semilogy(R, tec0, label='Electron Cold Temperature', linestyle='--')
plt.semilogy(R, Te_electrons, label='Electron Kappa Temperature')
plt.semilogy(R, ti0, label='O⁺ & S²⁺ Cold Temperature', linestyle='--')
plt.semilogy(R, Te_O, label='O⁺ Kappa Temperature')
plt.semilogy(R, Te_S, label='S²⁺ Kappa Temperature')
plt.xlabel(r'$\rho_c$ ($R_J$)', fontsize=14, fontweight='bold')
plt.ylabel('Temperature (eV)', fontsize=14, fontweight='bold')
plt.title('Nominal Product Isotropic Kappa Temperatures', fontsize=16, fontweight='bold')
plt.grid(True)
plt.tight_layout()

# Get current axes
ax = plt.gca()
# Set tick parameters
ax.tick_params(axis='both', which='major', labelsize=12, width=2)
# Set tick label font weight
for label in ax.get_xticklabels():
    label.set_fontweight('bold')
for label in ax.get_yticklabels():
    label.set_fontweight('bold')

# Set legend with bold font
plt.legend(fontsize=12, prop={'weight': 'bold'})

plt.savefig('new_nominal_model_product_kappa_temperatures_vs_radial_distance_012_momentsequal_bimax_to_product_kappa.png', dpi=600)
plt.show()

# Load data for the second plot
r_vs_kappa_steffl2004b = np.loadtxt('r_vs_kappa_steffl2004b.csv', delimiter=',', skiprows=1)
r_vs_kappa_steffl2004b_r = r_vs_kappa_steffl2004b[:, 0]
r_vs_kappa_steffl2004b_kappa = r_vs_kappa_steffl2004b[:, 1]

# Second plot
plt.figure(figsize=(10, 6))
plt.semilogy(R, kappa_electrons, label='Electron κ')
plt.semilogy(R, kappa_O, label='O⁺ κ')
plt.semilogy(R, kappa_S, label='S²⁺ κ')
plt.semilogy(r_vs_kappa_steffl2004b_r, r_vs_kappa_steffl2004b_kappa, label='Steffl (2004b) Electron κ')
plt.xlabel(r'$\rho_c$ ($R_J$)', fontsize=14, fontweight='bold')
plt.ylabel('Kappa Value (κ)', fontsize=14, fontweight='bold')
plt.title('Nominal Product Isotropic Kappa Values', fontsize=16, fontweight='bold')
plt.grid(True)
plt.tight_layout()

# Get current axes
ax = plt.gca()
# Set tick parameters
ax.tick_params(axis='both', which='major', labelsize=12, width=2)
# Set tick label font weight
for label in ax.get_xticklabels():
    label.set_fontweight('bold')
for label in ax.get_yticklabels():
    label.set_fontweight('bold')

# Set legend with bold font
plt.legend(fontsize=12, prop={'weight': 'bold'})

plt.savefig('new_nominal_model_product_kappa_values_vs_radial_distance_012_momentsequal_bimax_to_product_kappa.png', dpi=600)
plt.show()


print('Computation completed and results saved to CSV files.')
