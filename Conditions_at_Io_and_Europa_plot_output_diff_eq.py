import numpy as np
import matplotlib.pyplot as plt

# Define the case names and filenames
case_names = ['A', 'B', 'C', 'D', 'E', 'F']

n_out_filenames = [
    'new_nominal_model_n_out_new_nominal_model_4-10_iso_max_all.npz',                           # Case A
    'new_nominal_model_n_out_new_nominal_model_4-10_aniso2_max_all.npz',                        # Case B
    'new_nominal_model_n_out_nominal_model_4-10_iso_standard_kappa_all.npz',                # Case C
    'new_nominal_model_n_out_nominal_model_4-10_aniso2_standard_kappa_all.npz',             # Case D
    'new_nominal_model_n_out_nominal_model_4-10_isoT_isokappa_product_kappa_all_nominal_model.npz',  # Case E
    'new_nominal_model_n_out_nominal_model_4-10_aniso_fried_egg_all_nominal_model.npz'    # Case F
]

field_data_filenames = [
    'new_nominal_model_field_data_new_nominal_model_4-10_iso_max_all.npz',                           # Case A
    'new_nominal_model_field_data_new_nominal_model_4-10_aniso2_max_all.npz',                        # Case B
    'new_nominal_model_field_data_nominal_model_4-10_iso_standard_kappa_all.npz',                # Case C
    'new_nominal_model_field_data_nominal_model_4-10_aniso2_standard_kappa_all.npz',             # Case D
    'new_nominal_model_field_data_nominal_model_4-10_isoT_isokappa_product_kappa_all_nominal_model.npz',  # Case E
    'new_nominal_model_field_data_nominal_model_4-10_aniso_fried_egg_all_nominal_model.npz'    # Case F
]

cases = {}

for case_name, n_out_filename, field_data_filename in zip(case_names, n_out_filenames, field_data_filenames):
    # Load n_out
    n_out_loaded = np.load(n_out_filename)
    n_out = {key: n_out_loaded[key] for key in n_out_loaded}

    # Load field_data
    field_data = np.load(field_data_filename)
    
    cases[case_name] = {'n_out': n_out, 'field_data': field_data}

# For Case A, extract densities along field lines 190 and 539
case_A = cases['A']
n_out_A = case_A['n_out']
field_data_A = case_A['field_data']

x_out = field_data_A['x_out']
y_out = field_data_A['y_out']
z_out = field_data_A['z_out']

# Field line indices
field_line_index_190 = 189  # Indexing from 0
field_line_index_539 = 538

# Species keys and labels (including hot species)
species_keys = ['op', 'o2p', 'sp', 's2p', 's3p', 'hp', 'nap', 'oph', 'eh', 'elec']
species_labels = ['O$^{+}$', 'O$^{++}$', 'S$^{+}$', 'S$^{++}$', 'S$^{+++}$', 'H$^{+}$', 'Na$^{+}$', 'O$^{+}$ (hot)', 'e$^{-}$ (hot)', 'e$^{-}$']

# Extract data for field line 190
x_line_190 = x_out[field_line_index_190]
y_line_190 = y_out[field_line_index_190]
z_line_190 = z_out[field_line_index_190]

# Compute s along the field line
dx_190 = np.diff(x_line_190)
dy_190 = np.diff(y_line_190)
dz_190 = np.diff(z_line_190)
ds_190 = np.sqrt(dx_190**2 + dy_190**2 + dz_190**2)
s_190 = np.concatenate(([0], np.cumsum(ds_190)))

# Find equator index
equator_index_190 = np.argmin(np.abs(z_line_190))
s_190 = s_190 - s_190[equator_index_190]

# Only consider s >= 0
positive_s_indices_190 = s_190 >= 0
s_positive_190 = s_190[positive_s_indices_190]

# Limit s to 0 <= s <= 7.2 RJ for Io field line
s_limit_190 = s_positive_190 <= 7.2
s_positive_190 = s_positive_190[s_limit_190]

# Extract densities for all species
densities_190 = {}
for key in species_keys:
    density_line = n_out_A[key][field_line_index_190]
    density_line_positive = density_line[positive_s_indices_190]
    density_line_positive = density_line_positive[s_limit_190]
    densities_190[key] = density_line_positive

# Similarly for field line 539
x_line_539 = x_out[field_line_index_539]
y_line_539 = y_out[field_line_index_539]
z_line_539 = z_out[field_line_index_539]

dx_539 = np.diff(x_line_539)
dy_539 = np.diff(y_line_539)
dz_539 = np.diff(z_line_539)
ds_539 = np.sqrt(dx_539**2 + dy_539**2 + dz_539**2)
s_539 = np.concatenate(([0], np.cumsum(ds_539)))

equator_index_539 = np.argmin(np.abs(z_line_539))
s_539 = s_539 - s_539[equator_index_539]

positive_s_indices_539 = s_539 >= 0
s_positive_539 = s_539[positive_s_indices_539]

# Limit s to s <= 10.5 RJ for Europa field line
s_limit_539 = s_positive_539 <= 10.5
s_positive_539 = s_positive_539[s_limit_539]

densities_539 = {}
for key in species_keys:
    density_line = n_out_A[key][field_line_index_539]
    density_line_positive = density_line[positive_s_indices_539]
    density_line_positive = density_line_positive[s_limit_539]
    densities_539[key] = density_line_positive

# Now, get electron densities for all cases along field lines 190 and 539
electron_densities_190 = {}
electron_densities_539 = {}

for case_name in case_names:
    n_out = cases[case_name]['n_out']
    field_data = cases[case_name]['field_data']
    
    x_out = field_data['x_out']
    y_out = field_data['y_out']
    z_out = field_data['z_out']
    
    # Field line 190
    x_line = x_out[field_line_index_190]
    y_line = y_out[field_line_index_190]
    z_line = z_out[field_line_index_190]
    
    dx = np.diff(x_line)
    dy = np.diff(y_line)
    dz = np.diff(z_line)
    ds = np.sqrt(dx**2 + dy**2 + dz**2)
    s = np.concatenate(([0], np.cumsum(ds)))
    
    equator_index = np.argmin(np.abs(z_line))
    s = s - s[equator_index]
    
    positive_s_indices = s >= 0
    s_positive = s[positive_s_indices]
    
    s_limit = s_positive <= 7.2  # Limit s to 0 <= s <= 7.2 RJ
    s_positive = s_positive[s_limit]
    
    density_line = n_out['elec'][field_line_index_190]
    density_line_positive = density_line[positive_s_indices]
    density_line_positive = density_line_positive[s_limit]
    
    electron_densities_190[case_name] = {'s': s_positive, 'density': density_line_positive}
    
    # Field line 539
    x_line = x_out[field_line_index_539]
    y_line = y_out[field_line_index_539]
    z_line = z_out[field_line_index_539]
    
    dx = np.diff(x_line)
    dy = np.diff(y_line)
    dz = np.diff(z_line)
    ds = np.sqrt(dx**2 + dy**2 + dz**2)
    s = np.concatenate(([0], np.cumsum(ds)))
    
    equator_index = np.argmin(np.abs(z_line))
    s = s - s[equator_index]
    
    positive_s_indices = s >= 0
    s_positive = s[positive_s_indices]
    
    s_limit = s_positive <= 10.5  # Limit s to s <= 10.5 RJ
    s_positive = s_positive[s_limit]
    
    density_line = n_out['elec'][field_line_index_539]
    density_line_positive = density_line[positive_s_indices]
    density_line_positive = density_line_positive[s_limit]
    
    electron_densities_539[case_name] = {'s': s_positive, 'density': density_line_positive}

# Now, create the plots
fig, axs = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(12, 10))
plt.subplots_adjust(hspace=0, wspace=0)

ax_top_left = axs[0, 0]
ax_top_right = axs[0, 1]
ax_bottom_left = axs[1, 0]
ax_bottom_right = axs[1, 1]
legend_properties = {'weight':'bold'}

# Plot densities of all species for Case A along field line 190
ax = ax_top_left
for key, label in zip(species_keys, species_labels):
    density = densities_190[key]
    ax.plot(s_positive_190, density, label=label)
ax.set_yscale('log')
ax.set_ylim([0.1, 4000])
ax.set_xlim([0, 7.2])
ax.text(0.1, 0.975,'Io Field Line Iso-Maxwellian Densities', transform=ax.transAxes, fontsize=12, verticalalignment='top', fontweight='bold')
ax.legend(prop = legend_properties,fontsize=10, loc='center right')
ax.grid(True)
ax.tick_params(axis='both', which='both', direction='in', length=4, labelsize=12)
ax.tick_params(axis='both', which='major', labelsize=12)
for label in ax.get_xticklabels() + ax.get_yticklabels():
    label.set_fontweight('bold')

# Plot densities of all species for Case A along field line 539
ax = ax_top_right
for key, label in zip(species_keys, species_labels):
    density = densities_539[key]
    ax.plot(s_positive_539, density, label=label)
ax.set_yscale('log')
#ax.set_xscale('log')
ax.set_ylim([0.1, 4000])
ax.set_xlim([0.1, 10.5])
ax.text(0.05, 0.975, 'Europa Field Line Iso-Maxwellian Densities', transform=ax.transAxes, fontsize=12, verticalalignment='top', fontweight='bold')
# ax.legend(fontsize=10)
ax.grid(True)
ax.tick_params(axis='both', which='both', direction='in', length=4, labelsize=12)
ax.tick_params(axis='both', which='major', labelsize=12)
for label in ax.get_xticklabels() + ax.get_yticklabels():
    label.set_fontweight('bold')
legend_properties = {'weight':'bold'}
# Plot electron densities for all cases along field line 190
ax = ax_bottom_left
case_labels = ['Case A', 'Case B', 'Case C', 'Case D', 'Case E', 'Case F']
for case_name, label in zip(case_names, case_labels):
    s = electron_densities_190[case_name]['s']
    density = electron_densities_190[case_name]['density']
    ax.plot(s, density, label=label)
ax.set_yscale('log')
#ax.set_xscale('log')
ax.set_ylim([0.1, 4000])
ax.set_xlim([0.1, 7.2])
ax.text(0.2, 0.975, r'n$\bf{_{e}}$ Along Io Field Line', transform=ax.transAxes, fontsize=12, verticalalignment='top', fontweight='bold')
ax.legend(prop = legend_properties,fontsize=10, loc='lower left')
ax.grid(True)
ax.tick_params(axis='both', which='both', direction='in', length=4, labelsize=12)
ax.tick_params(axis='both', which='major', labelsize=12)
for label in ax.get_xticklabels() + ax.get_yticklabels():
    label.set_fontweight('bold')

# Plot electron densities for all cases along field line 539
ax = ax_bottom_right
for case_name, label in zip(case_names, case_labels):
    s = electron_densities_539[case_name]['s']
    density = electron_densities_539[case_name]['density']
    ax.plot(s, density, label=label)
ax.set_yscale('log')
#ax.set_xscale('log')
ax.set_ylim([0.1, 4000])
ax.set_xlim([0.1, 10.5])
ax.text(0.2, 0.975, r'n$\bf{_{e}}$ Along Europa Field Line', transform=ax.transAxes, fontsize=12, verticalalignment='top', fontweight='bold')
# ax.legend(fontsize=10)
ax.grid(True)
ax.tick_params(axis='both', which='both', direction='in', length=4, labelsize=12)
ax.tick_params(axis='both', which='major', labelsize=12)
for label in ax.get_xticklabels() + ax.get_yticklabels():
    label.set_fontweight('bold')

# Add common x and y axis labels
fig.text(0.5, 0.04, r's ($\bf{R_J}$)', ha='center', fontsize=16, fontweight='bold')
fig.text(0.04, 0.5, r'Density (cm$\bf{^{-3}}$)', va='center', rotation='vertical', fontsize=16, fontweight='bold')
plt.savefig(fname='new_nominal_model_log_Io_and_Europa_Comparison_plot.png', dpi=600, bbox_inches='tight')
plt.show()
