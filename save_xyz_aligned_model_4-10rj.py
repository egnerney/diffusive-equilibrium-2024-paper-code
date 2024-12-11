# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 00:34:22 2024

@author: Owner
"""

import numpy as np

# Define the case names and filenames
case_names = ['A', 'B', 'C', 'D', 'E', 'F']

field_data_filenames = [
    'new_nominal_model_field_data_new_nominal_model_4-10_iso_max_all.npz',                           # Case A
    'new_nominal_model_field_data_new_nominal_model_4-10_aniso2_max_all.npz',                        # Case B
    'new_nominal_model_field_data_nominal_model_4-10_iso_standard_kappa_all.npz',                    # Case C
    'new_nominal_model_field_data_nominal_model_4-10_aniso2_standard_kappa_all.npz',                 # Case D
    'new_nominal_model_field_data_nominal_model_4-10_isoT_isokappa_product_kappa_all_nominal_model.npz',  # Case E
    'new_nominal_model_field_data_nominal_model_4-10_aniso_fried_egg_all_nominal_model.npz'          # Case F
]

# Loop over each case
for case_idx, case_name in enumerate(case_names):
    print(f'Processing case {case_name}...')
    
    # Load field_data for the current case
    field_data_filename = field_data_filenames[case_idx]
    field_data = np.load(field_data_filename)
    
    # Extract x_out, y_out, z_out
    x_out = field_data['x_out']   # Shape: (601, N_points)
    y_out = field_data['y_out']
    z_out = field_data['z_out']
    
    # Ensure that the arrays are in the correct format
    # Replace any NaN values with np.nan (they should already be np.nan if they are missing)
    x_out = np.where(np.isfinite(x_out), x_out, np.nan)
    y_out = np.where(np.isfinite(y_out), y_out, np.nan)
    z_out = np.where(np.isfinite(z_out), z_out, np.nan)
    
    # Define output filenames
    x_csv_filename = f'case_{case_name}_x.csv'
    y_csv_filename = f'case_{case_name}_y.csv'
    z_csv_filename = f'case_{case_name}_z.csv'
    
    # Save to CSV files using numpy.savetxt
    # We need to handle the possibility of NaN values in the data
    # We will use '%.6e' format for floating point numbers
    np.savetxt(x_csv_filename, x_out, delimiter=',', fmt='%.6e', comments='', newline='\n')
    np.savetxt(y_csv_filename, y_out, delimiter=',', fmt='%.6e', comments='', newline='\n')
    np.savetxt(z_csv_filename, z_out, delimiter=',', fmt='%.6e', comments='', newline='\n')
    
    print(f'Saved x coordinates to {x_csv_filename}')
    print(f'Saved y coordinates to {y_csv_filename}')
    print(f'Saved z coordinates to {z_csv_filename}')

print('All cases processed successfully.')
