#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  9 02:55:04 2024

@author: edne8319
"""


import numpy as np
from scipy.special import gamma
from scipy.optimize import minimize
import matplotlib.pyplot as plt


# Load the radial profiles from CSV files
nop = np.loadtxt('nop_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv', delimiter=',')
ns2p = np.loadtxt('ns2p_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv', delimiter=',')
noph = np.loadtxt('noph_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv', delimiter=',')


ns2ph = 0.5*noph
noph = 0.5*noph 



nc = nop

nh = noph 





tc = np.loadtxt('Tic_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv', delimiter=',')
th = np.loadtxt('Toph_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv', delimiter=',')
kappa_values = np.loadtxt('kappa_vals_current.csv', delimiter=',')



# Step 2: Smooth the Kappa Values
def smooth_kappa(kappa, window_size=5):
    smoothed_kappa = np.convolve(kappa, np.ones(window_size)/window_size, mode='same')
    smoothed_kappa = np.where(smoothed_kappa < 1.51, 1.51, smoothed_kappa)
    return smoothed_kappa

# Smooth kappa values
kappa_values_smoothed = smooth_kappa(kappa_values)

# Step 3: Define the Functions
def maxwellian_superposition(v, nc, nh, tc, th, m):
    part1 = (np.exp(-m * v**2 / (2 * tc)) * (nc * (m / tc)**1.5)) / (2 * np.sqrt(2) * np.pi**1.5)
    part2 = (np.exp(-m * v**2 / (2 * th)) * (nh * (m / th)**1.5)) / (2 * np.sqrt(2) * np.pi**1.5)
    return part1 + part2

def kappa_distribution(v, nc, nh, kappa, T, m):
    if T <= 0 or kappa <= 0:
        return np.inf
    n_total = nc + nh
    factor = (m**1.5 * n_total * gamma(kappa)) / (8 * np.pi**1.5 * np.sqrt(kappa) * T**1.5 * gamma(kappa - 0.5))
    dist = factor * (1 + (m * v**2) / (4 * kappa * T**2))**(-kappa - 1)
    if np.any(np.isnan(dist)) or np.any(np.isinf(dist)):
        print(f"Invalid values in kappa distribution: T={T}, kappa={kappa}, factor={factor}, dist={dist}")
    return dist

def residual(T, v, nc, nh, tc, th, kappa, m):
    maxwellian = maxwellian_superposition(v, nc, nh, tc, th, m)
    kappa_dist = kappa_distribution(v, nc, nh, kappa, T, m)
    if np.any(np.isnan(maxwellian)) or np.any(np.isinf(maxwellian)):
        print(f"Invalid values in maxwellian: nc={nc}, nh={nh}, tc={tc}, th={th}")
    return np.sum((maxwellian - kappa_dist)**2)

# Step 4: Perform the Fit with Constraints
m = 1.0  # mass of the particles (adjust as needed)
v = np.linspace(0, 10, 100)  # velocity range (adjust as needed)
T_fit = []

for i in range(len(nc)):
    nc_i = nc[i]
    nh_i = nh[i]
    tc_i = tc[i]
    th_i = th[i]
    kappa_i = kappa_values_smoothed[i]
    
    T0 = (tc_i + th_i) / 2
    bounds = [(tc_i, None)]  # T must be greater than or equal to tc_i
    
    result = minimize(residual, T0, args=(v, nc_i, nh_i, tc_i, th_i, kappa_i, m), bounds=bounds)
    T_fit.append(result.x[0])

T_fit = np.array(T_fit)
np.savetxt('fitted_temperatures_op.csv', T_fit, delimiter=',')

# Step 5: Generate and Save the Plots

# Generate the rho array
rho = np.linspace(4, 10, 601)

# Plot tc, th, and T_fit vs rho
plt.figure(figsize=(10, 6))
plt.plot(rho, tc, label='T_c', linestyle='-', marker='o')
plt.plot(rho, th, label='T_h', linestyle='-', marker='x')
plt.plot(rho, T_fit, label='T_fit', linestyle='-', marker='s')

# Adding labels and legend
plt.xlabel('rho')
plt.ylabel('Temperature')
plt.title('Temperature Profiles vs. rho')
plt.legend()
plt.grid(True)

# Save the plot to a file
plt.savefig('temperature_profiles_vs_rho.png')

# Show the plot
plt.show()

# Plot nc and nh vs rho
plt.figure(figsize=(10, 6))
plt.plot(rho, nc, label='n_c', linestyle='-', marker='o')
plt.plot(rho, nh, label='n_h', linestyle='-', marker='x')

# Adding labels and legend
plt.xlabel('rho')
plt.ylabel('Density')
plt.title('Density Profiles vs. rho')
plt.legend()
plt.grid(True)

# Save the plot to a file


# Save the plot to a file
plt.savefig('density_profiles_vs_rho.png')

# Show the plot
plt.show()


    
    # Generate the rho array
ho = np.linspace(4, 10, 601)

# Define the specific rho values to plot
rho_values_to_plot = [5.7, 6, 7, 9.4, 10]
indices_to_plot = [np.abs(rho - val).argmin() for val in rho_values_to_plot]

# Plot the Maxwellians and Kappa distributions for the specified rho values
for idx in indices_to_plot:
    nc_i = nc[idx]
    nh_i = nh[idx]
    tc_i = tc[idx]
    th_i = th[idx]
    kappa_i = kappa_values_smoothed[idx]
    T_i = T_fit[idx]
    
    maxwellian = maxwellian_superposition(v, nc_i, nh_i, tc_i, th_i, m)
    kappa_dist = kappa_distribution(v, nc_i, nh_i, kappa_i, T_i, m)
    
    plt.figure(figsize=(10, 6))
    plt.plot(v, maxwellian, label='Maxwellian Superposition', linestyle='-', marker='o')
    plt.plot(v, kappa_dist, label=f'Kappa Distribution (kappa={kappa_i:.2f}, T={T_i:.2f})', linestyle='-', marker='x')
    
    plt.xlabel('v')
    plt.ylabel('Distribution')
    plt.title(f'Distribution vs. v for rho={rho[idx]}')
    plt.legend()
    plt.grid(True)
    
    # Save the plot to a file
    plt.savefig(f'distribution_vs_v_rho_{rho[idx]}.png')
    
    # Show the plot
    plt.show()