#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 13:48:06 2024

@author: edne8319
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gamma
from scipy.optimize import differential_evolution

# Physical constants
m_e = 9.10938356e-31        # Electron mass in kg
e = 1.602176634e-19         # Elementary charge in C (J/eV)
pi = np.pi

# Given parameters for the double Maxwellian
n_ec_cm3 = 30.5             # Core electron density in cm^-3
n_eh_cm3 = 1.52             # Hot electron density in cm^-3
T_ec_eV = 17.79             # Core electron temperature in eV
T_eh_eV = 299.11            # Hot electron temperature in eV

# Convert densities to m^-3
n_ec = n_ec_cm3 * 1e6       # Core electron density in m^-3
n_eh = n_eh_cm3 * 1e6       # Hot electron density in m^-3
n_e = n_ec + n_eh           # Total electron density in m^-3

# Convert temperatures to Joules
T_ec = T_ec_eV * e          # Core temperature in J
T_eh = T_eh_eV * e          # Hot temperature in J

# Velocity grid
v_max = 2.5 * np.sqrt(2 * T_eh / m_e)   # Maximum velocity for the grid
v = np.linspace(0, v_max, 1000)       # Linearly spaced velocity array

# To avoid plotting log(0) for f_e(v=0), set the first velocity to a small non-zero value
v[0] = v[1]*1e-3  # Set v[0] to a small non-zero value

# Constants for the double Maxwellian distributions
A_ec = n_ec * (m_e / (2 * pi * T_ec))**(1.5)
A_eh = n_eh * (m_e / (2 * pi * T_eh))**(1.5)

# Double Maxwellian distribution
f_double = A_ec * np.exp(- m_e * v**2 / (2 * T_ec)) + A_eh * np.exp(- m_e * v**2 / (2 * T_eh))

# Assume 10% uncertainty on f_double
percent_err = np.linspace(0.08,0.1,len(f_double))
sigma = percent_err * f_double

# Define the kappa distribution function
def kappa_distribution(v, T_e_eV, kappa_e, n_e):
    T_e = T_e_eV * e  # Convert to Joules
    if kappa_e <= 1.5:
        return np.zeros_like(v)
    A_kappa = n_e * (m_e / (2 * pi * T_e))**(1.5) * gamma(kappa_e) / ( np.sqrt(kappa_e) * gamma(kappa_e - 0.5) )
    f_kappa = A_kappa * (1 + m_e * v**2 / (2 * kappa_e * T_e))**( - (kappa_e + 1))
    return f_kappa

# Define the chi-squared function to minimize
def chi_squared(params, v, f_double, sigma, n_e):
    T_e_eV, kappa_e = params
    # Enforce parameter bounds within function
    if not (1 <= T_e_eV <= 500 and 1.505 <= kappa_e <= 100):
        return np.inf
    f_kappa = kappa_distribution(v, T_e_eV, kappa_e, n_e)
    # Avoid division by zero or log of zero
    f_kappa = np.where(f_kappa > 0, f_kappa, 1e-300)
    f_double_nonzero = np.where(f_double > 0, f_double, 1e-300)
    sigma_nonzero = np.where(sigma > 0, sigma, 1e-300)
    chi2 = np.sum( ((f_kappa - f_double_nonzero) / sigma_nonzero)**2 )
    return chi2

# Set parameter bounds: T_e from 1 to 500 eV, kappa_e from 1.505 to 100
bounds = [(17.789, 17.791), (1.505, 15.)]

# Perform differential evolution optimization
result = differential_evolution(chi_squared, bounds, args=(v, f_double, sigma, n_e), strategy='best1bin', maxiter=1000, tol=1e-6)

# Extract best-fit parameters
best_T_e_eV, best_kappa_e = result.x

# Compute best-fit kappa distribution
f_kappa_fit = kappa_distribution(v, best_T_e_eV, best_kappa_e, n_e)

# Plotting
plt.figure(figsize=(10, 6))
plt.semilogy((0.5 * m_e * v ** 2.) / e, f_double, label='Double Maxwellian', linewidth=2)
plt.semilogy((0.5 * m_e * v ** 2.) / e, f_kappa_fit, label=f'Best-fit Kappa (κ={best_kappa_e:.2f}, Tₑ={best_T_e_eV:.2f} eV)', linestyle='--', linewidth=2)
plt.xlabel('Electron Velocity $v$ (m/s)', fontsize=14)
plt.ylabel('Distribution Function $f_e(v)$', fontsize=14)
plt.title('Fitting Kappa Distribution to Double Maxwellian', fontsize=16)
plt.legend(fontsize=12)
plt.grid(True, which="both", ls="--")
plt.tight_layout()
plt.savefig('kappa_fit_fixed_Te.png',dpi = 1000)
plt.show()

# Print best-fit parameters
print(f'Best-fit T_e = {best_T_e_eV:.2f} eV')
print(f'Best-fit kappa_e = {best_kappa_e:.2f}')
