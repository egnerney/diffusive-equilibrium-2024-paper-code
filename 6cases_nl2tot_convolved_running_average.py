import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.integrate import romb
from multiprocessing import Pool
import multiprocessing

# Constants
R_J_m = 7.1492e7             # Jupiter radius in meters
R_J_cm = R_J_m * 100         # Convert R_J to centimeters

# Ion species keys (excluding electrons)
ion_species_keys = [
    'op',
    'o2p',
    'sp',
    's2p',
    's3p',
    'hp',
    'nap',
    'oph'
]

# Define the function to compute flux tube content for a field line
def compute_flux_tube_content(args):
    # Unpack the arguments
    i, x_out, y_out, z_out, n_out, ion_species_keys, R_J_cm = args
    
    # Get positions along the field line
    x_line = x_out[i]
    y_line = y_out[i]
    z_line = z_out[i]
    
    # Remove NaN values
    valid_indices = np.isfinite(x_line) & np.isfinite(y_line) & np.isfinite(z_line)
    x_line = x_line[valid_indices]
    y_line = y_line[valid_indices]
    z_line = z_line[valid_indices]
    
    if len(x_line) < 17:  # Need at least 17 points for romb with k=4
        return None  # Skip if not enough points
    
    # Find the equator index (z=0)
    equator_index = np.argmin(np.abs(z_line))
    
    # Take only from equator to northern extent
    x_line = x_line[equator_index:]
    y_line = y_line[equator_index:]
    z_line = z_line[equator_index:]
    
    if len(x_line) < 17:
        return None  # Skip if not enough points
    
    # Positions are in R_J units; keep ds in R_J units
    dx = np.diff(x_line)
    dy = np.diff(y_line)
    dz = np.diff(z_line)
    ds_RJ = np.sqrt(dx**2 + dy**2 + dz**2)  # in R_J units
    
    # Compute s along the field line in R_J units
    s_RJ = np.concatenate(([0], np.cumsum(ds_RJ)))  # in R_J units
    
    # Convert ds to centimeters
    ds_cm = ds_RJ * R_J_cm  # in cm
    
    # Compute s along the field line in cm
    s_cm = np.concatenate(([0], np.cumsum(ds_cm)))  # in cm
    
    # Prepare the s_cm and n_line_cm3 arrays for Romberg integration
    # We need s_cm to be evenly spaced, and number of samples to be 2^k + 1
    len_s_cm = len(s_cm)
    max_k = int(np.floor(np.log2(len_s_cm - 1)))
    if max_k < 4:  # Minimum k for romb is 4 (n=17)
        return None  # Skip if not enough points
    n_samples = 2**max_k + 1
    s_cm_uniform = np.linspace(s_cm[0], s_cm[-1], n_samples)
    
    # Initialize total number density integral
    total_integral = 0.0
    
    # Sum over all ion species (excluding electrons)
    for species_key in ion_species_keys:
        # Get number density for the species along the field line
        n_line = n_out[species_key][i]
        n_line = n_line[valid_indices]
        n_line = n_line[equator_index:]  # Take only northern extent
        n_line_cm3 = np.nan_to_num(n_line, nan=0.0)
        n_line_cm3[n_line_cm3 < 0] = 0.0  # Ensure non-negative
        
        # Interpolate n_line_cm3 onto s_cm_uniform
        interp_func = interp1d(s_cm, n_line_cm3, kind='linear', fill_value="extrapolate", bounds_error=False)
        n_line_cm3_uniform = interp_func(s_cm_uniform)
        
        # Use Romberg integration
        dx = (s_cm_uniform[1] - s_cm_uniform[0])/R_J_cm  # Uniform spacing
        integral_species = romb(n_line_cm3_uniform, dx=dx)
        
        # Sum over species
        total_integral += integral_species  # Units: cm⁻²
    
    # Equatorial radial distance r_eq (assuming equator is at z = 0)
    x_eq = x_line[0]  # At equator index
    y_eq = y_line[0]
    r_eq = np.sqrt(x_eq**2 + y_eq**2)  # in R_J units
    
    # L-shell value (in R_J units)
    L = r_eq  # Since r_eq is already in R_J units
    
    # Compute total flux tube content (NL²)_total
    flux_tube_content = 4. * np.pi * (R_J_cm)**3. * L**3. * total_integral  # Units: cm
    
    if L < 4.2:
        flux_tube_content = 0.0
        
    # Check if flux_tube_content is finite and reasonable
    if np.isfinite(flux_tube_content) and np.isfinite(r_eq) and flux_tube_content > 0:
        return (r_eq, flux_tube_content)
    else:
        return None

if __name__ == '__main__':
    # Include this line for Windows compatibility
    multiprocessing.freeze_support()
    
    # Define the case names and filenames
    case_names = ['A', 'B', 'C', 'D', 'E', 'F']
    case_labels = ['Case A', 'Case B', 'Case C', 'Case D', 'Case E', 'Case F']
    
    n_out_filenames = [
        'new_nominal_model_n_out_new_nominal_model_4-10_iso_max_all.npz',                           # Case A
        'new_nominal_model_n_out_new_nominal_model_4-10_aniso2_max_all.npz',                        # Case B
        'new_nominal_model_n_out_nominal_model_4-10_iso_standard_kappa_all.npz',                    # Case C
        'new_nominal_model_n_out_nominal_model_4-10_aniso2_standard_kappa_all.npz',                 # Case D
        'new_nominal_model_n_out_nominal_model_4-10_isoT_isokappa_product_kappa_all_nominal_model.npz',  # Case E
        'new_nominal_model_n_out_nominal_model_4-10_aniso_fried_egg_all_nominal_model.npz'          # Case F
    ]
    
    field_data_filenames = [
        'new_nominal_model_field_data_new_nominal_model_4-10_iso_max_all.npz',                      # Case A
        'new_nominal_model_field_data_new_nominal_model_4-10_aniso2_max_all.npz',                   # Case B
        'new_nominal_model_field_data_nominal_model_4-10_iso_standard_kappa_all.npz',               # Case C
        'new_nominal_model_field_data_nominal_model_4-10_aniso2_standard_kappa_all.npz',            # Case D
        'new_nominal_model_field_data_nominal_model_4-10_isoT_isokappa_product_kappa_all_nominal_model.npz',  # Case E
        'new_nominal_model_field_data_nominal_model_4-10_aniso_fried_egg_all_nominal_model.npz'     # Case F
    ]
    
    # Load data inside the main block
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
        
        num_field_lines = x_out.shape[0]
        
        # Prepare data for multiprocessing
        field_line_indices = range(num_field_lines)
        
        # Create a list of arguments for each field line
        args_list = [(i, x_out, y_out, z_out, n_out, ion_species_keys, R_J_cm) for i in field_line_indices]
        
        # Use multiprocessing Pool to compute flux tube content in parallel
        with Pool(processes=multiprocessing.cpu_count()) as pool:
            results = pool.map(compute_flux_tube_content, args_list)
        
        # Collect results
        r_eq_array = []
        flux_tube_content_array = []
        for result in results:
            if result is not None:
                r_eq, flux_tube_content = result
                r_eq_array.append(r_eq)
                flux_tube_content_array.append(flux_tube_content)
        
        # Convert lists to arrays
        r_eq_array = np.array(r_eq_array)
        flux_tube_content_array = np.array(flux_tube_content_array)
        
        # Sort by r_eq for plotting
        sorted_indices = np.argsort(r_eq_array)
        r_eq_array_sorted = r_eq_array[sorted_indices]
        flux_tube_content_array_sorted = flux_tube_content_array[sorted_indices]
        
        # Split data based on L = 5.6 R_J
        L_split = 5.6
        indices_low_L = np.where(r_eq_array_sorted <= L_split)[0]
        indices_high_L = np.where(r_eq_array_sorted > L_split)[0]
        
        # Ensure there are data points in both ranges
        flux_tube_content_smoothed = []
        r_eq_array_smoothed = []
        
        # For L <= 5.6 R_J
        if len(indices_low_L) > 0:
            window_size_low_L = 25
            flux_low_L = flux_tube_content_array_sorted[indices_low_L]
            r_eq_low_L = r_eq_array_sorted[indices_low_L]
            
            # Ensure window_size does not exceed data length
            window_size_low_L = min(window_size_low_L, len(flux_low_L))
            window_size_low_L = max(1, window_size_low_L)
            
            # Compute running average
            flux_low_L_smoothed = np.convolve(
                flux_low_L, 
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
            flux_tube_content_smoothed.append(flux_low_L_smoothed)
            r_eq_array_smoothed.append(r_eq_low_L_smoothed)
        
        # For L > 5.6 R_J
        if len(indices_high_L) > 0:
            window_size_high_L = 4
            flux_high_L = flux_tube_content_array_sorted[indices_high_L]
            r_eq_high_L = r_eq_array_sorted[indices_high_L]
            
            # Ensure window_size does not exceed data length
            window_size_high_L = min(window_size_high_L, len(flux_high_L))
            window_size_high_L = max(1, window_size_high_L)
            
            # Compute running average
            flux_high_L_smoothed = np.convolve(
                flux_high_L, 
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
            flux_tube_content_smoothed.append(flux_high_L_smoothed)
            r_eq_array_smoothed.append(r_eq_high_L_smoothed)
        
        # Concatenate smoothed data
        flux_tube_content_smoothed = np.concatenate(flux_tube_content_smoothed)
        r_eq_array_smoothed = np.concatenate(r_eq_array_smoothed)
        
        # Sort the concatenated data
        sorted_indices_smoothed = np.argsort(r_eq_array_smoothed)
        r_eq_array_smoothed = r_eq_array_smoothed[sorted_indices_smoothed]
        flux_tube_content_smoothed = flux_tube_content_smoothed[sorted_indices_smoothed]
        
        # Plot the smoothed data
        plt.plot(r_eq_array_smoothed, flux_tube_content_smoothed, label=case_labels[case_idx])
        
    legend_properties = {'weight':'bold'}
    # Customize the plot
    plt.xlabel(r'$\bf{\rho_c}$ ($\bf{R_J}$)', fontsize=14, fontweight='bold')
    plt.ylabel(r'Total Ion NL$\bf{^2}$', fontsize=14, fontweight='bold')
    plt.title(r'Total Ion NL$\bf{^2}$ vs $\bf{\rho_c}$', fontsize=16, fontweight='bold')
    plt.legend(prop=legend_properties, fontsize=12)
    plt.grid(True)

    # Set x-limits from 4 to 10 $R_J$
    plt.xlim(4, 10)

    # Set y-axis to logarithmic scale
    plt.yscale('log')  # Logarithmic y-axis

    # Make tickmarks and tickmark labels bold
    ax = plt.gca()
    ax.tick_params(axis='both', which='both', direction='in', length=4, labelsize=12)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight('bold')

    # Save the figure
    plt.savefig(fname='total_ion_nl2_vs_L_logscale.png', dpi=600, bbox_inches='tight')
    plt.show()
