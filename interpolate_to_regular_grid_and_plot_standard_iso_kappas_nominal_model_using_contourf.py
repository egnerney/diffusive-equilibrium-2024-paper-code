# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 14:54:00 2024

@author: Owner
"""
import numpy as np
from scipy.optimize import bisect
from scipy.interpolate import interp1d
import sys
import matplotlib.pyplot as plt
import matplotlib.colors as colors
#from scipy.ndimage import uniform_filter1d  # For smoothing
import scipy.special
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.image as mpimg
import matplotlib.tri as tri
from scipy.interpolate import griddata

def isnotfinite(arr):
    res = np.isfinite(arr)
    np.bitwise_not(res, out=res)  # in-place
    return res

# Load n_out

n_out_loaded_isos = np.load('new_nominal_model_n_out_nominal_model_4-10_iso_standard_kappa_all.npz')
n_out_isos = {key: n_out_loaded_isos[key] for key in n_out_loaded_isos}

n_out_loaded = np.load('new_nominal_model_n_out_nominal_model_4-10_aniso2_standard_kappa_all.npz')
n_out = {key: n_out_loaded[key] for key in n_out_loaded}

# Load T_out
T_out_loaded = np.load('new_nominal_model_T_out_nominal_model_4-10_aniso2_standard_kappa_all.npz')
T_out = {key: T_out_loaded[key] for key in T_out_loaded}

# Load field data
field_data = np.load('new_nominal_model_field_data_nominal_model_4-10_aniso2_standard_kappa_all.npz')
x_out = field_data['x_out']
y_out = field_data['y_out']
z_out = field_data['z_out']
B_out = field_data['B_out']

r_out = np.sqrt(x_out ** 2. + y_out ** 2. + z_out ** 2.)
rho_out = np.sqrt(x_out ** 2. + y_out ** 2. )
lat_out_deg = np.degrees(np.arcsin(z_out / r_out))
    
# Define the species keys and labels in the desired order
#species_keys = ['op', 'o2p', 'sp', 's2p', 's3p', 'hp', 'nap', 'oph', 'eh', 'elec']
#species_labels = ['O$^{+}$', 'O$^{++}$', 'S$^{+}$', 'S$^{++}$', 'S$^{+++}$', 'H$^{+}$', 'Na$^{+}$', 'O$^{+}$ (hot)', 'e$^{-}$ (hot)', 'e$^{-}$']

#numdens_species_labels = [r'$\bf{n_{O^{+}}}$', r'$\bf{n_{O^{++}}}$', r'$\bf{n_{S^{+}}}$', r'$\bf{n_{S^{++}}}$', r'$\bf{n_{S^{+++}}}$', r'$\bf{n_{H^{+}}}$', r'$\bf{n_{Na^{+}}}$', r'$\bf{n_{O^{+}}}$ (hot)', r'$\bf{n_{e^{-}}}$ (hot)', r'$\bf{n_{e^{-}}}$']

species_keys = ['op', 'o2p', 'sp', 's2p', 's3p', 'hp', 'nap', 'elec']
species_labels = ['O$^{+}$', 'O$^{++}$', 'S$^{+}$', 'S$^{++}$', 'S$^{+++}$', 'H$^{+}$', 'Na$^{+}$', 'e$^{-}$']

numdens_species_labels = [r'$\bf{n_{O^{+}}}$', r'$\bf{n_{O^{++}}}$', r'$\bf{n_{S^{+}}}$', r'$\bf{n_{S^{++}}}$', r'$\bf{n_{S^{+++}}}$', r'$\bf{n_{H^{+}}}$', r'$\bf{n_{Na^{+}}}$', r'$\bf{n_{e^{-}}}$']

# Flatten the coordinate and density arrays
rho_flat = rho_out.flatten()
z_flat = z_out.flatten()
r_flat = r_out.flatten()


# Load the image of Jupiter
jupiter_img = mpimg.imread('hubble-captures-vivid-auroras-in-jupiters-atmosphere_28000029525_o~large.png')

# Check if the image has an alpha channel; if not, add one
if jupiter_img.shape[2] == 3:
    # Add an alpha channel
    alpha_channel = np.ones((jupiter_img.shape[0], jupiter_img.shape[1], 1))
    jupiter_img = np.concatenate((jupiter_img, alpha_channel), axis=2)

# Adjust the image for Jupiter's oblateness
# Equatorial radius: 1 R_J = 71492 km
# Polar radius: approximately 0.93512 R_J = 66,854 km
oblateness_factor = 0.93512  # Adjust as needed
jupiter_height = 2 * oblateness_factor
# Set the extent of the image: [left, right, bottom, top] in data coordinates
img_extent = [-1.0, 1.0, -oblateness_factor, oblateness_factor]

# Create a figure with adjusted size and margins to accommodate labels and aspect ratio
fig = plt.figure(figsize=(11.25, 10))  # Adjust the height as needed
nrows = 4
ncols = 2
gs = fig.add_gridspec(nrows, ncols,left=0.075,bottom=0.05,hspace=0,wspace=0)#, left=0.075, right=1.0, top=0.95, bottom=0.05, wspace=0.0, hspace=0.0)
#gs = fig.add_gridspec(nrows, ncols,hspace=0,wspace=0)#, left=0.075, right=1.0, top=0.95, bottom=0.05, wspace=0.0, hspace=0.0)


(ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8) = gs.subplots(sharex='col', sharey='row')
axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8]

nlevels = 1000

minnec=1.
maxnec=n_out_isos['elec'].max()

minnsp=0.1
maxnsp=n_out_isos['sp'].max()

minns2p=0.1
maxns2p=n_out_isos['s2p'].max()

minns3p=0.1
maxns3p=n_out_isos['s3p'].max()

minnop=0.1
maxnop=n_out_isos['op'].max()

minno2p=0.1
maxno2p=n_out_isos['o2p'].max()

minnnap=0.1
maxnnap=n_out_isos['nap'].max()

minnhp=0.5
maxnhp=n_out_isos['hp'].max()









# Set the levels for contour plots
levelsnec = np.logspace(np.log10(minnec), np.log10(maxnec), nlevels)

# Create a Normalize object for logarithmic scale
normnec = colors.LogNorm(vmin=minnec, vmax=maxnec)


levelsnsp = np.logspace(np.log10(minnsp), np.log10(maxnsp), nlevels)
normnsp = colors.LogNorm(vmin=minnsp, vmax=maxnsp)


levelsns2p = np.logspace(np.log10(minns2p), np.log10(maxns2p), nlevels)
normns2p = colors.LogNorm(vmin=minnsp, vmax=maxnsp)


levelsns3p = np.logspace(np.log10(minns3p), np.log10(maxns3p), nlevels)
normns3p = colors.LogNorm(vmin=minns3p, vmax=maxns3p)


levelsnop = np.logspace(np.log10(minnop), np.log10(maxnop), nlevels)
normnop = colors.LogNorm(vmin=minnop, vmax=maxnop)

levelsno2p = np.logspace(np.log10(minno2p), np.log10(maxno2p), nlevels)
normno2p = colors.LogNorm(vmin=minno2p, vmax=maxno2p)



levelsnap = np.logspace(np.log10(minnnap), np.log10(maxnnap), nlevels)
normnap = colors.LogNorm(vmin=minnnap, vmax=maxnnap)

levelsnhp = np.logspace(np.log10(minnhp), np.log10(maxnhp), nlevels)
normnhp = colors.LogNorm(vmin=minnhp, vmax=maxnhp)

levelss = [levelsnop,levelsno2p,levelsnsp, levelsns2p, levelsns3p, levelsnhp, levelsnap, levelsnec]
normss = [normnop,normno2p,normnsp, normns2p, normns3p,normnhp, normnap, normnec]


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

"""
cbar_ticks_nhp = [
    0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
    1, 2, 3, 4, 5, 6, 7, 8, 9,
    10, 20, 30, 40, 50, 60
]
cbar_label_ticks_nhp = [0.1, 1, 10,60]
"""
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




cbar_tickss = [cbar_ticks_nop,cbar_ticks_no2p,cbar_ticks_nsp,cbar_ticks_ns2p,cbar_ticks_ns3p,cbar_ticks_nhp,cbar_ticks_nnap,cbar_ticks_nec]
cbar_label_tickss = [cbar_label_ticks_nop,cbar_label_ticks_no2p,cbar_label_ticks_nsp,cbar_label_ticks_ns2p,cbar_label_ticks_ns3p,cbar_label_ticks_nhp,cbar_label_ticks_nnap,cbar_label_ticks_nec]
cbar_tick_labelss = [cbar_tick_labels_nop,cbar_tick_labels_no2p,cbar_tick_labels_nsp,cbar_tick_labels_ns2p,cbar_tick_labels_ns3p,cbar_tick_labels_nhp,cbar_tick_labels_nnap,cbar_tick_labels_nec]
# Plot each species


#for idx, key in enumerate(species_keys):
#old_species_keys = ['op', 'o2p', 'sp', 's2p', 's3p', 'hp', 'nap', 'oph', 'eh', 'elec']
#species_keys = ['op', 'o2p', 'sp', 's2p', 's3p', 'hp', 'nap', 'elec']
idxs =[7,2,0,3,1,4,5,6]
# Define contour levels for black contour lines for each species
contour_levelss_old = [
    [1,  100],        # For 'op'
    [ 1, 30],        # For 'o2p'
    [1,  100],        # For 'sp'
    [1,  100],        # For 's2p'
    [1, 30],        # For 's3p'
    [5, 20],        # For 'hp'
    [1,  30],        # For 'nap'
    [10, 100, 1000]      # For 'elec'
]

contour_levelss_old2 = [
    [1, 50, 100,500],        # For 'op'
    [ 1,10, 30,50],        # For 'o2p'
    [1, 10,50, 100,500],        # For 'sp'
    [1, 10,50,300],        # For 's2p'
    [0.5,5, 20,50],        # For 's3p'
    [5,15,35],        # For 'hp'
    [1,  30],        # For 'nap'
    [10, 100, 500,1000,2000]      # For 'elec'
]


contour_levelss = [
    [1,100],        # For 'op'
    [1, 30],        # For 'o2p'
    [1,100],        # For 'sp'
    [1,100],        # For 's2p'
    [1,30],        # For 's3p'
    [5,20],        # For 'hp'
    [1,  30],        # For 'nap'
    [10, 100,1000]      # For 'elec'
]




    

for iii in range(8):
    idx = idxs[iii]
    ax = axes[iii]
    key = species_keys[idx]
    density = n_out[key].flatten()
    levels = levelss[idx]
    norm = normss[idx]
    cbar_ticks = cbar_tickss[idx]
    cbar_label_ticks = cbar_label_tickss[idx]
    cbar_tick_labels = cbar_tick_labelss[idx] 
    contour_levels = contour_levelss[idx]

    # Remove NaN or negative values and apply additional conditions
 
    valid = (
        np.isfinite(density) & (density > 0) &
        (r_flat > 1.6) &
        (z_flat >= -3) & (z_flat <= 3) & (rho_flat > 4.4*(rho_flat ** 3.) / ((z_flat ** 2. + rho_flat ** 2.) ** (3./2.))) & (rho_flat < 9.9*(rho_flat ** 3.) / ((z_flat ** 2. + rho_flat ** 2.) ** (3./2.))))
    
    rho_valid = rho_flat[valid]
    z_valid = z_flat[valid]
    #density[nvalid] = 0.00001
    density_valid = density[valid]
    
    # Define the grid
    x_min, x_max = 0., 12.0
    y_min, y_max = -3.0, 3.0
    num_rho_points = 12000
    num_z_points = 600
    rho_lin = np.linspace(x_min, x_max, num_rho_points)
    z_lin = np.linspace(y_min, y_max, num_z_points)
    rho_grid, z_grid = np.meshgrid(rho_lin, z_lin)

    # Compute r_grid
    r_grid = np.sqrt(rho_grid ** 2 + z_grid ** 2)

    # Interpolate the density onto the grid (use absolute value of rho for symmetry)
    points = np.column_stack((rho_valid, z_valid))
    values = density_valid
    density_grid = griddata(points, values, (np.abs(rho_grid), z_grid), method='linear')

    # Replace NaN values with zero
    density_grid = np.nan_to_num(density_grid, nan=0.0)

    # Compute valid_grid
    valid_grid = (
        (r_grid > 1.6) &
        (z_grid >= -3) & (z_grid <= 3) &
        (np.abs(rho_grid) > 4.4 * (np.abs(rho_grid) ** 3) / ((z_grid ** 2 + np.abs(rho_grid) ** 2) ** (1.5))) &
        (np.abs(rho_grid) < 9.9 * (np.abs(rho_grid) ** 3) / ((z_grid ** 2 + np.abs(rho_grid) ** 2) ** (1.5)))
    )

    # Set density to zero in non-valid areas
    density_grid[~valid_grid] = 0.0

    # Plot using contourf
    contour = ax.contourf(rho_grid, z_grid, density_grid, levels=levels, norm=norm, cmap='viridis')
    # Define contour levels for black contour lines for each species
    # Add black contour lines with labels
    contour_lines = ax.contour(rho_grid, z_grid, density_grid, levels=contour_levels, colors='black', linewidths=1)
    ax.clabel(contour_lines, fmt='%g', inline=True, fontsize=8, colors='black')


        
    # Set limits
    x_min, x_max = -1.5, 12.0
    y_min, y_max = -3.0, 3.0
    ax.set_xlim([x_min, x_max])
    ax.set_ylim([y_min, y_max])
    # Set aspect ratio to 'equal' to make 1 R_J the same length on both axes
    ax.set_aspect('equal')
    

    # Enable tick marks
    ax.tick_params(axis='both', which='both', direction='in', length=4)

    # Manage x-axis tick labels
    if iii // ncols == nrows - 1:
        # Bottom row plots: keep x-axis tick labels
        x_ticks = np.array([0,1,  2, 3, 4,5,  6,7, 8,9,10])#np.arange(-1.5, 13.5, 3)  # Adjust as needed
        ax.set_xticks(x_ticks)
        ax.tick_params(labelbottom=True)
    else:
        x_ticks = np.array([0,1,  2, 3, 4,5,  6,7, 8,9,10])#np.arange(-1.5, 13.5, 3)  # Adjust as needed
        ax.set_xticks(x_ticks)
        ax.tick_params(labelbottom=False)
        
    # Manage y-axis tick labels
    if iii % ncols == 0:
        y_ticks = np.array([-3,-2,-1,0,1,2])#np.arange(-3, 4, 1)  # Adjust as needed
        ax.set_yticks(y_ticks)
        ax.tick_params(labelleft=True)
    else:
        y_ticks = np.array([-3,-2,-1,0,1,2])#np.arange(-3, 4, 1)  # Adjust as needed
        ax.set_yticks(y_ticks)
        # Other plots: remove y-axis tick labels
        ax.tick_params(labelleft=False)
        
    for label in ax.get_xticklabels():
        label.set_fontweight('bold')
    for label in ax.get_yticklabels():
        label.set_fontweight('bold')
        
        

    # Plot the species label inside the plot at rho=3, z=0
    ax.text(
        3.0, 0.0,numdens_species_labels[idx], fontsize=12 , fontweight='bold' ,
        ha='center', va='center', color='black')#,
    ax.text(
        3.0, -0.7, '(cm$^{-3}$)', fontsize=12, fontweight='bold',
        ha='center', va='center', color='black')#,
       # bbox=dict(facecolor='black', alpha=0.6, edgecolor='none'))
    # Calculate the position and size of the colorbar in axes fraction coordinates
    # Desired position: centered at rho=11, z=0, width=1.5 RJ, height=5.5 RJ
    x_span = x_max - x_min
    y_span = y_max - y_min

    x_center = (10.3 - x_min) / x_span  # Center of colorbar in axes fraction
    x_width = 0.4 / x_span  # Width of colorbar in axes fraction
    x_left = x_center - x_width / 2

    y_center = (0 - y_min) / y_span  # Center of colorbar in axes fraction
    y_height = 5.5 / y_span  # Height of colorbar in axes fraction
    y_bottom = y_center - y_height / 2

    # Create an inset_axes for the colorbar
    cax = ax.inset_axes([x_left, y_bottom, x_width, y_height], transform=ax.transAxes)

    # Add colorbar
    cbar = fig.colorbar(contour, cax=cax, orientation='vertical')

    # Set the colorbar ticks and labels
    cbar.set_ticks(cbar_ticks)
    cbar_tick_labels_list = [cbar_tick_labels.get(tick, '') for tick in cbar_ticks]
    cbar.ax.set_yticklabels(cbar_tick_labels_list)
    # Set colorbar tick labels to bold
    for label in cbar.ax.get_yticklabels():
        label.set_fontweight('bold')

    # Adjust the colorbar appearance
    cbar.ax.tick_params(labelsize=9)
    #cbar.outline.set_visible(False)  # Optional: hide the colorbar outline
    cbar.set_label('')  # Remove colorbar label to avoid overlap
    
    # Overlay the image
    ax.imshow(
        jupiter_img,
        extent=img_extent,
        aspect='auto',
        origin='upper',  # Adjust if the image appears upside down
        zorder=10  # Draw the image on top
    )

    # Optional: Adjust the transparency of the image
    # You can set alpha parameter if needed, e.g., alpha=0.8

# Add common x and y axis labels
fig.text(0.49, 0.02, '$ρ_c$ ($R_J$)', ha='center', fontsize=16, fontweight='bold')
fig.text(0.02, 0.5, '$z_c$ ($R_J$)', va='center', rotation='vertical', fontsize=16, fontweight='bold')



# Save the figure with resolution 500 dpi
plt.savefig('caseD_final.png', dpi=600, bbox_inches='tight', pad_inches=0)

# Show the plot if running interactively
plt.show()