# Import necessary libraries
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec

r_00 = 0.01 * np.arange(601) + 4.0  
# Load more data files (ensure these files exist)
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


# Load your data (ensure that these files are in your working directory)
n_out_loaded = np.load('n_out_nominal_model_4-10_iso_standard_kappa_all_nominal_model.npz')
n_out = {key: n_out_loaded[key] for key in n_out_loaded}

T_out_loaded = np.load('T_out_nominal_model_4-10_iso_standard_kappa_all_nominal_model.npz')
T_out = {key: T_out_loaded[key] for key in T_out_loaded}

field_data = np.load('field_data_nominal_model_4-10_iso_standard_kappa_all_nominal_model.npz')
x_out = field_data['x_out']
y_out = field_data['y_out']
z_out = field_data['z_out']
B_out = field_data['B_out']

# Compute coordinate arrays
r_out = np.sqrt(x_out ** 2. + y_out ** 2. + z_out ** 2.)
rho_out = np.sqrt(x_out ** 2. + y_out ** 2.)
lat_out_deg = np.degrees(np.arcsin(z_out / r_out))

# Species keys and labels, including 'oph' and 'eh'
species_keys = ['op', 'o2p', 'sp', 's2p', 's3p', 'hp', 'nap', 'oph', 'eh', 'elec']
species_labels = ['O$^{+}$', 'O$^{++}$', 'S$^{+}$', 'S$^{++}$', 'S$^{+++}$', 'H$^{+}$', 'Na$^{+}$', 'O$^{+}$ (hot)', 'e$^{-}$ (hot)', 'e$^{-}$']
numdens_species_labels = [r'$\bf{n_{O^{+}}}$', r'$\bf{n_{O^{++}}}$', r'$\bf{n_{S^{+}}}$', r'$\bf{n_{S^{++}}}$', r'$\bf{n_{S^{+++}}}$', r'$\bf{n_{H^{+}}}$', r'$\bf{n_{Na^{+}}}$', r'$\bf{n_{O^{+}}}$ (hot)', r'$\bf{n_{e^{-}}}$ (hot)', r'$\bf{n_{e^{-}}}$']

# Flatten the coordinate arrays
rho_flat = rho_out.flatten()
z_flat = z_out.flatten()
r_flat = r_out.flatten()

# Index mapping
idxs = [9, 8, 0, 2, 1, 3, 7, 4, 5, 6]


fig = plt.figure(figsize=(9, 10))  # Adjust the height as needed
nrows = 5
ncols = 2
gs = fig.add_gridspec(nrows, ncols,left=0.11,bottom=0.09,hspace=0)#, left=0.075, right=1.0, top=0.95, bottom=0.05, wspace=0.0, hspace=0.0)
#gs = fig.add_gridspec(nrows, ncols,hspace=0,wspace=0)#, left=0.075, right=1.0, top=0.95, bottom=0.05, wspace=0.0, hspace=0.0)


(ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8), (ax9, ax10) = gs.subplots(sharex='col')
axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8,ax9,ax10]


minnec=1.
maxnec=n_out['elec'].max()

minnsp=0.1
maxnsp=n_out['sp'].max()

minns2p=0.1
maxns2p=n_out['s2p'].max()

minns3p=0.1
maxns3p=n_out['s3p'].max()

minnop=0.1
maxnop=n_out['op'].max()

minno2p=0.1
maxno2p=n_out['o2p'].max()

minnnap=0.1
maxnnap=n_out['nap'].max()

minnhp=0.5
maxnhp=n_out['hp'].max()

minnoph=0.1
maxnoph=n_out['oph'].max()

minneh=1.0
maxneh=n_out['eh'].max()


#minss = [minnop,minno2p,minnsp,minns2p,minns3p,min]
#maxss = 






# Prepare the tick marks and labels for the colorbars
cbar_ticks_nec = [
    1, 2, 3, 4, 5, 6, 7, 8, 9,
    10, 20, 30, 40, 50, 60, 70, 80, 90,
    100, 200, 300, 400, 500, 600, 700, 800, 900,
    1000, 2000, 3000
]
cbar_label_ticks_nec = [ 1, 10, 100, 1000, 3000]
cbar_tick_labels_nec = {tick: str(tick) if tick in cbar_label_ticks_nec else '' for tick in cbar_ticks_nec}

cbar_ticks_nop = [
    0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
    1, 2, 3, 4, 5, 6, 7, 8, 9,
    10, 20, 30, 40, 50, 60, 70, 80, 90,
    100, 200, 300, 400, 500, 600, 700, 800, 900,
    1000
]
cbar_label_ticks_nop = [0.1, 1, 10, 100, 1000]
cbar_tick_labels_nop = {tick: str(tick) if tick in cbar_label_ticks_nop else '' for tick in cbar_ticks_nop}



cbar_ticks_no2p = [
    0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
    1, 2, 3, 4, 5, 6, 7, 8, 9,
    10, 20, 30, 40, 50, 60, 70, 80, 90
]
cbar_label_ticks_no2p = [0.1, 1, 10,90]
cbar_tick_labels_no2p = {tick: str(tick) if tick in cbar_label_ticks_no2p else '' for tick in cbar_ticks_no2p}

cbar_ticks_nsp = [
    0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
    1, 2, 3, 4, 5, 6, 7, 8, 9,
    10, 20, 30, 40, 50, 60, 70, 80, 90,
    100, 200, 300, 400, 500, 600, 700, 800, 900,
    1000
]
cbar_label_ticks_nsp = [0.1, 1, 10, 100, 1000]
cbar_tick_labels_nsp = {tick: str(tick) if tick in cbar_label_ticks_nsp else '' for tick in cbar_ticks_nsp}


cbar_ticks_ns2p = [
    0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
    1, 2, 3, 4, 5, 6, 7, 8, 9,
    10, 20, 30, 40, 50, 60, 70, 80, 90,
    100, 200, 300, 400, 500
]
cbar_label_ticks_ns2p = [0.1, 1, 10, 100, 500]
cbar_tick_labels_ns2p = {tick: str(tick) if tick in cbar_label_ticks_ns2p else '' for tick in cbar_ticks_ns2p}


cbar_ticks_ns3p = [
    0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
    1, 2, 3, 4, 5, 6, 7, 8, 9,
    10, 20, 30, 40, 50, 60
]
cbar_label_ticks_ns3p = [0.1, 1, 10,60]
cbar_tick_labels_ns3p = {tick: str(tick) if tick in cbar_label_ticks_ns3p else '' for tick in cbar_ticks_ns3p}


cbar_ticks_nhp = [
     0.5, 0.6, 0.7, 0.8, 0.9,
    1, 2, 3, 4, 5, 6, 7, 8, 9,
    10, 20, 30, 40, 50, 60
]
cbar_label_ticks_nhp = [0.5, 1, 10,60]
cbar_tick_labels_nhp = {tick: str(tick) if tick in cbar_label_ticks_nhp else '' for tick in cbar_ticks_nhp}

cbar_ticks_nnap = [
    0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
    1, 2, 3, 4, 5, 6, 7, 8, 9,
    10, 20, 30, 40, 50, 60, 70, 80, 90,
    100
]
cbar_label_ticks_nnap = [0.1, 1, 10, 100]
cbar_tick_labels_nnap = {tick: str(tick) if tick in cbar_label_ticks_nnap else '' for tick in cbar_ticks_nnap}


cbar_ticks_noph = [
    0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
    1, 2, 3, 4, 5, 6, 7, 8, 9,
    10, 20, 30, 40, 50, 60, 70, 80, 90,
    100, 200
]
cbar_label_ticks_noph = [0.1, 1, 10, 100,200]
cbar_tick_labels_noph = {tick: str(tick) if tick in cbar_label_ticks_noph else '' for tick in cbar_ticks_noph}


cbar_ticks_neh = [
    1, 2, 3, 4, 6 ,8 
]
cbar_label_ticks_neh = [ 1, 4,8]
cbar_tick_labels_neh = {tick: str(tick) if tick in cbar_label_ticks_neh else '' for tick in cbar_ticks_neh}

xx_ticks =  [4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10]
xx_label_ticks =  [4,5,6,7,8,9,10]
xaxis_tick_labels = {tick: str(tick) if tick in xx_label_ticks else '' for tick in xx_ticks}


cbar_tickss = [cbar_ticks_nop,cbar_ticks_no2p,cbar_ticks_nsp,cbar_ticks_ns2p,cbar_ticks_ns3p,cbar_ticks_nhp,cbar_ticks_nnap,cbar_ticks_noph,cbar_ticks_neh,cbar_ticks_nec]
cbar_label_tickss = [cbar_label_ticks_nop,cbar_label_ticks_no2p,cbar_label_ticks_nsp,cbar_label_ticks_ns2p,cbar_label_ticks_ns3p,cbar_label_ticks_nhp,cbar_label_ticks_nnap,cbar_label_ticks_noph,cbar_label_ticks_neh,cbar_label_ticks_nec]
cbar_tick_labelss = [cbar_tick_labels_nop,cbar_tick_labels_no2p,cbar_tick_labels_nsp,cbar_tick_labels_ns2p,cbar_tick_labels_ns3p,cbar_tick_labels_nhp,cbar_tick_labels_nnap,cbar_tick_labels_noph,cbar_tick_labels_neh,cbar_tick_labels_nec]
# Plot each species

nss = [nop0,no2p0,nsp0,ns2p0,ns3p0,nhp0,nnap0,noph0,neh0,nec0]

# Loop over each species
for iii in range(10):
    idx = idxs[iii]
    ax = axes[iii]
    key = species_keys[idx]
    density = nss[idx]
    cbar_ticks = cbar_tickss[idx]
    cbar_label_ticks = cbar_label_tickss[idx]
    cbar_tick_labels = cbar_tick_labelss[idx]
    mins, maxs = cbar_ticks[0], cbar_ticks[-1]*2
    if key == 'eh':
        mins, maxs = cbar_ticks[0], 7
        


    # Plot the density vs. rho_c
    ax.plot(r_00,density, linewidth=2)

    # Set logarithmic y-axis
    ax.set_yscale('log')

    # Set plot limits
    ax.set_xlim([4, 10])
    ax.set_ylim([mins, maxs])

    # Set tick marks and labels
    ax.tick_params(axis='both', which='both', direction='in', length=4)

    # Manage x-axis tick labels
    if ((iii == 8) or (iii == 9)):
        # Bottom row plots: keep x-axis tick labels
        x_ticks = xx_ticks  # Adjust as needed
        x_tick_labels_list = [xaxis_tick_labels.get(tick, '') for tick in x_ticks]
        ax.set_xticks(x_ticks)
        ax.set_xticklabels(x_tick_labels_list)
        ax.tick_params(labelbottom=True)
    else:
        x_ticks = xx_ticks
        ax.set_xticks(x_ticks)
        ax.tick_params(labelbottom=False)

    # Manage y-axis tick labels
    y_ticks = cbar_ticks  # Use the same ticks as colorbars
    y_tick_labels_list = [cbar_tick_labels.get(tick, '') for tick in y_ticks]
    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_tick_labels_list)
    ax.tick_params(labelleft=True)

    for label in ax.get_xticklabels():
        label.set_fontweight('bold')
    for label in ax.get_yticklabels():
        label.set_fontweight('bold')

    # Plot the species label inside the plot at a suitable position
    ax.text(
        0.95, 0.95, numdens_species_labels[idx], fontsize=16, fontweight='bold',
        ha='right', va='top', color='black', transform=ax.transAxes
    )

# Add common x and y axis labels
fig.text(0.5, 0.04, '$\\rho_c$ ($R_J$)', ha='center', fontsize=18, fontweight='bold')
fig.text(0.04, 0.5, r'$n_\alpha$ (cm$^{-3}$)', va='center', rotation='vertical', fontsize=18, fontweight='bold')

# Adjust layout
#plt.tight_layout(rect=[0.05, 0.05, 1, 1])  # Adjust rect to make space for common labels

# Save the figure
plt.savefig('new_nominal_model_radial_profiles_species_densities_nominal_model_4-10RJ.pdf', dpi=600, bbox_inches='tight', pad_inches=0)

# Show the plot
plt.show()