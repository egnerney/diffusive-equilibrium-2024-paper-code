# -*- coding: utf-8 -*-
"""
Created on Fri May  2 01:53:44 2025

@author: Owner
"""

import numpy as np
import matplotlib.pyplot as plt

############################################################
# If you do not need interactive plots, you can optionally:
# plt.ioff()
############################################################

# ----------------------------------------------------------------------------------
# 1) Setup: filenames for each of the 6 cases (A..F) for both field_data and n_out/T_out
# ----------------------------------------------------------------------------------

case_names = ['A', 'B', 'C', 'D', 'E', 'F']

field_data_files = {
    'A': 'final_new_nominal_model_field_data_nominal_model_4-10_isoT_A=1_max_all_nominal_model.npz',
    'B': 'final_new_nominal_model_field_data_nominal_model_4-10_anisoT_A=2_max_all_nominal_model.npz',
    'C': 'final_new_nominal_model_field_data_nominal_model_4-10_isoT_isokappa_A=1_standard_kappa_all_nominal_model.npz',
    'D': 'final_new_nominal_model_field_data_nominal_model_4-10_anisoT_isokappa_A=2_standard_kappa_all_nominal_model.npz',
    'E': 'final_new_nominal_model_field_data_nominal_model_4-10_isoT_isokappa_A=1_product_kappa_all_nominal_model.npz',
    'F': 'final_new_nominal_model_field_data_nominal_model_4-10_isoT_isokappa_A=1_fried_egg_all_nominal_model.npz',
}

density_files = {
    'A': 'final_new_nominal_model_n_out_nominal_model_4-10_isoT_A=1_max_all_nominal_model.npz',
    'B': 'final_new_nominal_model_n_out_nominal_model_4-10_anisoT_A=2_max_all_nominal_model.npz',
    'C': 'final_new_nominal_model_n_out_nominal_model_4-10_isoT_isokappa_A=1_standard_kappa_all_nominal_model.npz',
    'D': 'final_new_nominal_model_n_out_nominal_model_4-10_anisoT_isokappa_A=2_standard_kappa_all_nominal_model.npz',
    'E': 'final_new_nominal_model_n_out_nominal_model_4-10_isoT_isokappa_A=1_product_kappa_all_nominal_model.npz',
    'F': 'final_new_nominal_model_n_out_nominal_model_4-10_isoT_isokappa_A=1_fried_egg_all_nominal_model.npz',
}

temperature_files = {
    'A': 'final_new_nominal_model_T_out_nominal_model_4-10_isoT_A=1_max_all_nominal_model.npz',
    'B': 'final_new_nominal_model_T_out_nominal_model_4-10_anisoT_A=2_max_all_nominal_model.npz',
    'C': 'final_new_nominal_model_T_out_nominal_model_4-10_isoT_isokappa_A=1_standard_kappa_all_nominal_model.npz',
    'D': 'final_new_nominal_model_T_out_nominal_model_4-10_anisoT_isokappa_A=2_standard_kappa_all_nominal_model.npz',
    'E': 'final_new_nominal_model_T_out_nominal_model_4-10_isoT_isokappa_A=1_product_kappa_all_nominal_model.npz',
    'F': 'final_new_nominal_model_T_out_nominal_model_4-10_isoT_isokappa_A=1_fried_egg_all_nominal_model.npz',
}

# ----------------------------------------------------------------------------------
# 2) Field line indices for Io and Europa
# ----------------------------------------------------------------------------------
field_line_indices = [189, 538]
field_line_labels = {189: 'Io Field Line', 538: 'Europa Field Line'}

# For the northern hemisphere portion of the array
lat_start = 700
lat_end = 1400

# Additional "s" cutoffs
s_max_io = 7.2
s_max_europa = 11.39

# ----------------------------------------------------------------------------------
# 3) Species information
# ----------------------------------------------------------------------------------
# According to your original post, the keys for the arrays are:
# ['op', 'o2p', 'sp', 's2p', 's3p', 'hp', 'nap', 'oph', 'eh', 'elec']
#
# We'll define an order that matches your plotting logic (10 species).
species_keys = [
    'elec',       # e- (cold)
    'eh',         # e- (hot)
    'op',         # O+
    'sp',         # S+
    'o2p',        # O++
    's2p',        # S2+
    'oph',        # O+ (hot)
    's3p',        # S3+
    'hp',         # H+
    'nap'         # Na+
]
species_labels = [
    'e$^{-}$',     # for 'elec'
    'e$^{-}$(hot)',# for 'eh'
    'O$^{+}$',     # 'op'
    'S$^{+}$',     # 'sp'
    'O$^{++}$',    # 'o2p'
    'S$^{2+}$',    # 's2p'
    'O$^{+}$(hot)',# 'oph'
    'S$^{3+}$',    # 's3p'
    'H$^{+}$',     # 'hp'
    'Na$^{+}$',    # 'nap'
]

# ----------------------------------------------------------------------------------
# 4) Load data for each case; store in a nested dictionary "results"
# ----------------------------------------------------------------------------------
#   results[case_name][field_line_index] = {
#       's': array(...),
#       'x': array(...),
#       'y': array(...),
#       'z': array(...),
#       'B': array(...),
#       'phis': array(...),
#       'n_out':   { 'elec': ..., 'eh': ..., ... },
#       'T_para':  { 'elec': ..., 'eh': ..., ... },
#       'T_perp':  { 'elec': ..., 'eh': ..., ... },
#   }
results = {c: {} for c in case_names}

for case in case_names:
    print(f"Loading data for case: {case}")
    # 4a) Load field data
    fd = np.load(field_data_files[case])
    x_all = fd['x_out']    # shape (601, 1401)
    y_all = fd['y_out']
    z_all = fd['z_out']
    B_all = fd['B_out']
    phi_all = fd['phi_out']       # shape (601, 1401)

    # 4b) Load densities and temperatures
    n_loaded = np.load(density_files[case])
    T_loaded = np.load(temperature_files[case])

    # 4c) For each field line of interest (Io & Europa)
    for fl_index in field_line_indices:

        # Slice out northern hemisphere
        x_line = x_all[fl_index, lat_start:lat_end]
        y_line = y_all[fl_index, lat_start:lat_end]
        z_line = z_all[fl_index, lat_start:lat_end]
        B_line = B_all[fl_index, lat_start:lat_end]
        phi_line = phi_all[fl_index, lat_start:lat_end]

        # Arc length s
        dx = np.diff(x_line)
        dy = np.diff(y_line)
        dz = np.diff(z_line)
        ds = np.sqrt(dx**2 + dy**2 + dz**2)
        s_line = np.concatenate(([0.0], np.cumsum(ds)))

        # Shift s so that the equator is s=0
        eq_index = np.argmin(np.abs(z_line))
        s_line = s_line - s_line[eq_index]

        # Prepare containers
        n_dict   = {}
        Tpara_dict = {}
        Tperp_dict = {}

        # 4d) Loop over species
        for skey, _ in zip(species_keys, species_labels):

            # If the species does not exist in the NPZ (for some cases), store zeros
            if (skey not in n_loaded.keys()) or ((skey + "_par") not in T_loaded.keys()):
                n_dict[skey] = np.zeros(len(x_line))
                Tpara_dict[skey] = np.zeros(len(x_line))
                Tperp_dict[skey] = np.zeros(len(x_line))
            else:
                n_full = n_loaded[skey][fl_index, :]
                T_par_full = T_loaded[skey + "_par"][fl_index, :]
                T_perp_full = T_loaded[skey + "_perp"][fl_index, :]

                # Slice to northern hemisphere
                n_dict[skey] = n_full[lat_start:lat_end]
                Tpara_dict[skey] = T_par_full[lat_start:lat_end]
                Tperp_dict[skey] = T_perp_full[lat_start:lat_end]

        # -----------------------------------------------------------------
        # NEW STEP: set 'eh' & 'oph' T_para and T_perp to np.nan if case != A or B
        # -----------------------------------------------------------------
        for hot_key in ['eh', 'oph']:
            if case not in ['A','B']:
                Tpara_dict[hot_key][:] = np.nan
                Tperp_dict[hot_key][:] = np.nan

        # Store results
        results[case][fl_index] = {
            's': s_line,
            'x': x_line,
            'y': y_line,
            'z': z_line,
            'B': B_line,
            'phis': phi_line,
            'n_out': n_dict,
            'T_para': Tpara_dict,
            'T_perp': Tperp_dict
        }

print("All data loaded.")

# ----------------------------------------------------------------------------------
# 5) Create plots
# ----------------------------------------------------------------------------------

# Helper function to choose s-limit for Io vs. Europa
def get_s_limit(fline_index):
    if fline_index == 189:
        return s_max_io
    elif fline_index == 538:
        return s_max_europa
    else:
        return 1e9  # or any large number

# A) Plot Number Density
for fl_index in field_line_indices:
    fl_label = field_line_labels[fl_index]
    fig, axs = plt.subplots(5, 2, figsize=(14, 16), sharex=True)
    fig.suptitle(f'Densities along {fl_label} (Northern Hemisphere)', fontsize=20, fontweight='bold')

    for i_sp, (skey, slabel) in enumerate(zip(species_keys, species_labels)):
        row = i_sp // 2
        col = i_sp % 2
        ax = axs[row, col]

        for case in case_names:
            s_arr = results[case][fl_index]['s']
            n_arr = results[case][fl_index]['n_out'][skey]

            # Only s >= 0
            mask = s_arr >= 0
            s_pos = s_arr[mask]
            n_pos = n_arr[mask]

            # Also limit by s <= (7.2 or 11.39)
            s_lim = get_s_limit(fl_index)
            mask2 = s_pos <= s_lim
            s_plot = s_pos[mask2]
            n_plot = n_pos[mask2]

            # Plot if there's any nonzero data
            if np.any(n_plot > 0):
                ax.plot(s_plot, n_plot, label=f'Case {case}')

        # Label on the subplot
        if skey == 'eh':
            ax.text(0.5, 0.7, slabel, transform=ax.transAxes, fontsize=20, fontweight='bold', ha='center')
        else:
            ax.text(0.5, 0.8, slabel, transform=ax.transAxes, fontsize=20, fontweight='bold', ha='center')

        ax.grid(True, which='both', linestyle='--')
        ax.tick_params(axis='both', which='both', labelsize=17, width=2)
        if row == 4:
            ax.set_xlabel('s (R$_J$)', fontsize=20, fontweight='bold')

    fig.text(0.04, 0.5, 'Density (cm$^{-3}$)', va='center', rotation='vertical',
             fontsize=20, fontweight='bold')
    axs[3,1].legend(loc='center right', frameon=False, prop={'weight': 'bold', 'size':16})
    plt.tight_layout(rect=[0.05, 0.03, 1, 0.95])

    # Save as PDF, 600 dpi
    pdf_name = f'Densities_vs_s_{fl_label.replace(" ", "_")}_NH.pdf'
    plt.savefig(pdf_name, dpi=600)
    plt.show()

# B) Plot Parallel Temperatures
for fl_index in field_line_indices:
    fl_label = field_line_labels[fl_index]
    fig, axs = plt.subplots(5, 2, figsize=(14, 20), sharex=True)
    fig.suptitle(f'Parallel Temperatures along {fl_label} (Northern Hemisphere)', fontsize=20, fontweight='bold')

    for i_sp, (skey, slabel) in enumerate(zip(species_keys, species_labels)):
        row = i_sp // 2
        col = i_sp % 2
        ax = axs[row, col]

        for case in case_names:
            s_arr = results[case][fl_index]['s']
            Tpar_arr = results[case][fl_index]['T_para'][skey]

            mask = s_arr >= 0
            s_pos = s_arr[mask]
            Tpar_pos = Tpar_arr[mask]

            s_lim = get_s_limit(fl_index)
            mask2 = s_pos <= s_lim
            s_plot = s_pos[mask2]
            Tpar_plot = Tpar_pos[mask2]

            if np.any(Tpar_plot > 0):
                ax.plot(s_plot, Tpar_plot, label=f'Case {case}')

        if (skey == 'eh') or (skey == 'oph'):
            ax.text(0.5, 0.75, slabel, transform=ax.transAxes, fontsize=20, fontweight='bold', ha='center')
        else:
            ax.text(0.5, 0.9, slabel, transform=ax.transAxes, fontsize=20, fontweight='bold', ha='center')

        ax.grid(True, which='both', linestyle='--')
        ax.tick_params(axis='both', which='both', labelsize=17, width=2)
        if row == 4:
            ax.set_xlabel('s (R$_J$)', fontsize=20, fontweight='bold')

    fig.text(0.04, 0.5, 'Parallel Temperature (eV)', va='center', rotation='vertical',
             fontsize=20, fontweight='bold')
    axs[3,1].legend(loc='upper left', frameon=False, prop={'weight': 'bold', 'size':16})
    plt.tight_layout(rect=[0.05, 0.03, 1, 0.99])
    pdf_name = f'mean_Parallel_Temperatures_vs_s_{fl_label.replace(" ", "_")}_NH.pdf'
    plt.savefig(pdf_name, dpi=600)
    plt.show()

# C) Plot Perpendicular Temperatures
for fl_index in field_line_indices:
    fl_label = field_line_labels[fl_index]
    fig, axs = plt.subplots(5, 2, figsize=(14, 20), sharex=True)
    fig.suptitle(f'Perpendicular Temperatures along {fl_label} (Northern Hemisphere)', fontsize=20, fontweight='bold')

    for i_sp, (skey, slabel) in enumerate(zip(species_keys, species_labels)):
        row = i_sp // 2
        col = i_sp % 2
        ax = axs[row, col]

        for case in case_names:
            s_arr = results[case][fl_index]['s']
            Tperp_arr = results[case][fl_index]['T_perp'][skey]

            mask = s_arr >= 0
            s_pos = s_arr[mask]
            Tperp_pos = Tperp_arr[mask]

            s_lim = get_s_limit(fl_index)
            mask2 = s_pos <= s_lim
            s_plot = s_pos[mask2]
            Tperp_plot = Tperp_pos[mask2]

            if np.any(Tperp_plot > 0):
                ax.plot(s_plot, Tperp_plot, label=f'Case {case}')

        if (skey == 'eh') or (skey == 'oph'):
            ax.text(0.5, 0.75, slabel, transform=ax.transAxes, fontsize=20, fontweight='bold', ha='center')
        elif (skey == 'hp'):
            ax.text(0.5, 0.8, slabel, transform=ax.transAxes, fontsize=20, fontweight='bold', ha='center')
        else:
            ax.text(0.5, 0.9, slabel, transform=ax.transAxes, fontsize=20, fontweight='bold', ha='center')

        ax.grid(True, which='both', linestyle='--')
        ax.tick_params(axis='both', which='both', labelsize=17, width=2)
        if row == 4:
            ax.set_xlabel('s (R$_J$)', fontsize=20, fontweight='bold')

    fig.text(0.04, 0.5, 'Perpendicular Temperature (eV)', va='center', rotation='vertical',
             fontsize=20, fontweight='bold')
    axs[3,1].legend(loc='upper left', frameon=False, prop={'weight': 'bold', 'size':16})
    plt.tight_layout(rect=[0.05, 0.03, 1, 0.99])
    pdf_name = f'mean_Perpendicular_Temperatures_vs_s_{fl_label.replace(" ", "_")}_NH.pdf'
    plt.savefig(pdf_name, dpi=600)
    plt.show()

# D) Plot Potential phi
for fl_index in field_line_indices:
    fl_label = field_line_labels[fl_index]
    fig, ax = plt.subplots(figsize=(10, 6))

    for case in case_names:
        s_arr = results[case][fl_index]['s']
        phi_arr = results[case][fl_index]['phis']

        mask = s_arr >= 0
        s_pos = s_arr[mask]
        phi_pos = phi_arr[mask]

        s_lim = get_s_limit(fl_index)
        mask2 = s_pos <= s_lim
        s_plot = s_pos[mask2]
        phi_plot = phi_pos[mask2]

        # Plot
        ax.plot(s_plot, phi_plot, label=f'Case {case}')

    ax.set_xlabel('s (R$_J$)', fontsize=20, fontweight='bold')
    ax.set_ylabel(r'$\phi$ (V)', fontsize=20, fontweight='bold')
    ax.set_title(f'Electrostatic Potential along {fl_label} (Northern Hemisphere)',
                 fontsize=20, fontweight='bold')

    ax.grid(True, which='both', linestyle='--')
    ax.tick_params(axis='both', which='both', labelsize=17, width=2)
    ax.legend(frameon=False, prop={'weight': 'bold', 'size':16})

    plt.tight_layout()
    pdf_name = f'Phi_vs_s_{fl_label.replace(" ", "_")}_NH.pdf'
    plt.savefig(pdf_name, dpi=600)
    plt.show()

print("All plots generated and saved as PDFs with 600 dpi.")
