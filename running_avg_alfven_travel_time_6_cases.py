# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 01:48:48 2024

@author: Owner
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.integrate import romb

# Constants
mu_0 = 4 * np.pi * 1e-7       # Vacuum permeability (H/m)
c = 299792458                 # Speed of light (m/s)
u = 1.66053906660e-27         # Atomic mass unit (kg)
R_J = 7.1492e7                # Jupiter radius in meters

# Species masses in kg
masses = {
    'op': 16 * u,
    'o2p': 16 * u,
    'sp': 32 * u,
    's2p': 32 * u,
    's3p': 32 * u,
    'hp': 1 * u,
    'nap': 23 * u,
    'oph': 16 * u,
    'eh': 9.10938356e-31,
    'elec': 9.10938356e-31
}

# Define the case names and filenames
case_names = ['A', 'B', 'C', 'D', 'E', 'F']
case_labels = ['Case A', 'Case B', 'Case C', 'Case D', 'Case E', 'Case F']



field_data_filenames = [
     'final_new_nominal_model_field_data_nominal_model_4-10_isoT_A=1_max_all_nominal_model.npz',
     'final_new_nominal_model_field_data_nominal_model_4-10_anisoT_A=2_max_all_nominal_model.npz',
    'final_new_nominal_model_field_data_nominal_model_4-10_isoT_isokappa_A=1_standard_kappa_all_nominal_model.npz',
     'final_new_nominal_model_field_data_nominal_model_4-10_anisoT_isokappa_A=2_standard_kappa_all_nominal_model.npz',
    'final_new_nominal_model_field_data_nominal_model_4-10_isoT_isokappa_A=1_product_kappa_all_nominal_model.npz',
     'final_new_nominal_model_field_data_nominal_model_4-10_isoT_isokappa_A=1_fried_egg_all_nominal_model.npz',
]

n_out_filenames = [
    'final_new_nominal_model_n_out_nominal_model_4-10_isoT_A=1_max_all_nominal_model.npz',
    'final_new_nominal_model_n_out_nominal_model_4-10_anisoT_A=2_max_all_nominal_model.npz',
     'final_new_nominal_model_n_out_nominal_model_4-10_isoT_isokappa_A=1_standard_kappa_all_nominal_model.npz',
     'final_new_nominal_model_n_out_nominal_model_4-10_anisoT_isokappa_A=2_standard_kappa_all_nominal_model.npz',
     'final_new_nominal_model_n_out_nominal_model_4-10_isoT_isokappa_A=1_product_kappa_all_nominal_model.npz',
    'final_new_nominal_model_n_out_nominal_model_4-10_isoT_isokappa_A=1_fried_egg_all_nominal_model.npz',
]



cases = {}

for case_name, n_out_filename, field_data_filename in zip(case_names, n_out_filenames, field_data_filenames):
    # Load n_out
    n_out_loaded = np.load(n_out_filename)
    n_out = {key: n_out_loaded[key] for key in n_out_loaded}

    # Load field_data
    field_data = np.load(field_data_filename)
    
    cases[case_name] = {'n_out': n_out, 'field_data': field_data}

# Initialize a figure for plotting
plt.figure(figsize=(10, 6))

# Loop over each case
for case_idx, case_name in enumerate(case_names):
    n_out = cases[case_name]['n_out']
    field_data = cases[case_name]['field_data']
    
    # Extract field data
    x_out = field_data['x_out']   # Shape: (601, N_points)
    y_out = field_data['y_out']
    z_out = field_data['z_out']
    B_out = field_data['B_out']   # Magnetic field magnitude in nanotesla
    
    num_field_lines = x_out.shape[0]
    T_total_array = []
    r_eq_array = []
    
    # Loop over each field line
    for i in range(num_field_lines):
        # Get positions along the field line
        x_line = x_out[i]
        y_line = y_out[i]
        z_line = z_out[i]
        B_line = B_out[i]  # B in nanotesla
        
        # Remove NaN values
        valid_indices = np.isfinite(x_line) & np.isfinite(y_line) & np.isfinite(z_line) & np.isfinite(B_line)
        x_line = x_line[valid_indices]
        y_line = y_line[valid_indices]
        z_line = z_line[valid_indices]
        B_line = B_line[valid_indices]
        
        if len(x_line) < 17:
            continue  # Skip if not enough points for romb (requires at least 17 points for k=4)
        
        # Convert positions from R_J to meters
        x_line_m = x_line * R_J
        y_line_m = y_line * R_J
        z_line_m = z_line * R_J
        
        # Convert B from nanotesla to tesla
        B_line_T = B_line * 1e-9  # Now in tesla
        
        # Compute ds along the field line
        dx = np.diff(x_line_m)
        dy = np.diff(y_line_m)
        dz = np.diff(z_line_m)
        ds = np.sqrt(dx**2 + dy**2 + dz**2)  # in meters
        
        # Compute s along the field line
        s = np.concatenate(([0], np.cumsum(ds)))  # in meters
        
        # Prepare s_uniform for Romberg integration
        len_s = len(s)
        max_k = int(np.floor(np.log2(len_s - 1)))
        if max_k < 4:
            continue  # Skip this field line if not enough points for romb with k>=4
        n_points = 2**max_k + 1
        s_uniform = np.linspace(s[0], s[-1], n_points)
        
        # Interpolate B_line_T onto s_uniform
        B_interp = interp1d(s, B_line_T, kind='linear', fill_value="extrapolate", bounds_error=False)
        B_line_T_uniform = B_interp(s_uniform)
        B_line_T_uniform[B_line_T_uniform < 0] = 0.0  # Ensure non-negative
        
        # Compute rho_total at positions
        rho_total = np.zeros_like(s_uniform)  # Initialize total mass density array
        
        # Sum over all species
        for species_key, mass in masses.items():
            # Get number density for the species along the field line
            n_line = n_out[species_key][i]
            n_line = n_line[valid_indices]
            n_line_m3 = n_line * 1e6  # Convert from cm^-3 to m^-3
            
            # Replace NaNs or negative values with zeros
            n_line_m3 = np.nan_to_num(n_line_m3, nan=0.0)
            n_line_m3[n_line_m3 < 0] = 0.0
            
            # Interpolate n_line_m3 onto s_uniform
            n_interp = interp1d(s, n_line_m3, kind='linear', fill_value="extrapolate", bounds_error=False)
            n_line_m3_uniform = n_interp(s_uniform)
            n_line_m3_uniform[n_line_m3_uniform < 0] = 0.0  # Ensure non-negative
            
            rho_species = n_line_m3_uniform * mass  # Mass density for the species
            rho_total += rho_species  # Add to total mass density
        
        # Replace zeros or negative values in rho_total with a small positive value
        rho_total[rho_total <= 0] = 1e-20  # kg/m^3
        
        # Compute Alfvén speed v_A at points
        v_A = B_line_T_uniform / np.sqrt(mu_0 * rho_total)  # v_A in m/s
        
        # Compute relativistic Alfvén speed v_rel
        v_rel = v_A / np.sqrt(1 + (v_A / c)**2)
        v_rel[v_rel <= 0] = 1e-20  # Avoid division by zero
        
        # Compute integrand dt/ds = 1 / v_rel
        integrand = 1 / v_rel  # s/m
        integrand = np.nan_to_num(integrand, nan=0.0, posinf=0.0, neginf=0.0)
        
        # Use romb to compute the integral
        dx = s_uniform[1] - s_uniform[0]  # Uniform spacing in meters
        T_line = romb(integrand, dx=dx)  # T_line in seconds
        
        # Total travel time for two bounces
        T_total = 2 * T_line  # in seconds
        
        # Convert T_total to minutes
        T_total_minutes = T_total / 60.0
        
        # Equatorial radial distance r_eq (assuming equator is at z = 0)
        equator_index = np.argmin(np.abs(z_line))
        x_eq = x_line[equator_index]
        y_eq = y_line[equator_index]
        r_eq = np.sqrt(x_eq**2 + y_eq**2)  # in R_J
        
        # Append results if T_total is finite and reasonable
        if np.isfinite(T_total_minutes) and np.isfinite(r_eq) and T_total_minutes > 0 and T_total_minutes < 1e4:
            T_total_array.append(T_total_minutes)
            r_eq_array.append(r_eq)
    
    # Convert lists to arrays
    T_total_array = np.array(T_total_array)
    r_eq_array = np.array(r_eq_array)
    
    # Sort by r_eq for plotting
    sorted_indices = np.argsort(r_eq_array)
    r_eq_array_sorted = r_eq_array[sorted_indices]
    T_total_array_sorted = T_total_array[sorted_indices]
    
    # Apply running average with different window sizes based on L
    # Split data based on L = 5.6 R_J
    L_split = 5.6
    indices_low_L = np.where(r_eq_array_sorted <= L_split)[0]
    indices_high_L = np.where(r_eq_array_sorted > L_split)[0]
    
    # Initialize lists to store smoothed data
    T_total_smoothed = []
    r_eq_smoothed = []
    
    # For L <= 5.6 R_J
    if len(indices_low_L) > 0:
        window_size_low_L = 25
        T_low_L = T_total_array_sorted[indices_low_L]
        r_eq_low_L = r_eq_array_sorted[indices_low_L]
        
        # Ensure window_size does not exceed data length
        window_size_low_L = min(window_size_low_L, len(T_low_L))
        window_size_low_L = max(1, window_size_low_L)
        
        # Compute running average
        T_low_L_smoothed = np.convolve(
            T_low_L, 
            np.ones(window_size_low_L)/window_size_low_L, 
            mode='valid'
        )
        # Adjust r_eq_array accordingly
        r_eq_low_L_smoothed = np.convolve(
            r_eq_low_L, 
            np.ones(window_size_low_L)/window_size_low_L, 
            mode='valid'
        )
        
        # Store smoothed data
        T_total_smoothed.append(T_low_L_smoothed)
        r_eq_smoothed.append(r_eq_low_L_smoothed)
    
    # For L > 5.6 R_J
    if len(indices_high_L) > 0:
        window_size_high_L = 4
        T_high_L = T_total_array_sorted[indices_high_L]
        r_eq_high_L = r_eq_array_sorted[indices_high_L]
        
        # Ensure window_size does not exceed data length
        window_size_high_L = min(window_size_high_L, len(T_high_L))
        window_size_high_L = max(1, window_size_high_L)
        
        # Compute running average
        T_high_L_smoothed = np.convolve(
            T_high_L, 
            np.ones(window_size_high_L)/window_size_high_L, 
            mode='valid'
        )
        # Adjust r_eq_array accordingly
        r_eq_high_L_smoothed = np.convolve(
            r_eq_high_L, 
            np.ones(window_size_high_L)/window_size_high_L, 
            mode='valid'
        )
        
        # Store smoothed data
        T_total_smoothed.append(T_high_L_smoothed)
        r_eq_smoothed.append(r_eq_high_L_smoothed)
    
    # Concatenate smoothed data
    T_total_smoothed = np.concatenate(T_total_smoothed)
    r_eq_smoothed = np.concatenate(r_eq_smoothed)
    
    # Sort the concatenated data
    sorted_indices_smoothed = np.argsort(r_eq_smoothed)
    r_eq_smoothed = r_eq_smoothed[sorted_indices_smoothed]
    T_total_smoothed = T_total_smoothed[sorted_indices_smoothed]
    
    # Plot the smoothed data
    plt.plot(r_eq_smoothed, T_total_smoothed, label=case_labels[case_idx])

legend_properties = {'weight':'bold'}
# Customize the plot
plt.xlabel(r' $\bf{\rho_{ceq}}$ ($\bf{R_J}$)', fontsize=14, fontweight='bold')
plt.ylabel('Alfvén Travel Time (Minutes)', fontsize=14, fontweight='bold')
plt.title(r'Alfvén Travel Time, $\bf{T = 2(T_N + T_S)}$', fontsize=16, fontweight='bold')
plt.legend(prop=legend_properties, fontsize=12)
plt.grid(True)

# Set x-limits from 4 to 10 $R_J$
plt.xlim(4, 10)

# Make tickmarks and tickmark labels bold
ax = plt.gca()
ax.tick_params(axis='both', which='both', direction='in', length=4, labelsize=12)
for label in ax.get_xticklabels() + ax.get_yticklabels():
    label.set_fontweight('bold')

plt.savefig(fname='new_nominal_model_twice_total_north_plus_south_alfven_travel_time_6cases_romberg.pdf', dpi=600, bbox_inches='tight')
plt.show()
