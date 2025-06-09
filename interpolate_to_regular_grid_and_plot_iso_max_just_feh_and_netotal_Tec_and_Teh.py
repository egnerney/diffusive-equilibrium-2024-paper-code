# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 14:54:00 2024

@author: Owner
"""
import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.image as mpimg
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


def isnotfinite(arr):
    res = np.isfinite(arr)
    np.bitwise_not(res, out=res)  # in-place
    return res

# -------------------------------------------------------------------------
# 1) LOAD DATA
# -------------------------------------------------------------------------
# Load n_out
n_out_loaded_iso_max = np.load('final_new_nominal_model_n_out_nominal_model_4-10_isoT_A=1_max_all_nominal_model.npz')
n_out_iso_max = {key: n_out_loaded_iso_max[key] for key in n_out_loaded_iso_max}

n_out_loaded = np.load('final_new_nominal_model_n_out_nominal_model_4-10_isoT_A=1_max_all_nominal_model.npz')
n_out = {key: n_out_loaded[key] for key in n_out_loaded}

T_out_loaded = np.load('final_new_nominal_model_T_out_nominal_model_4-10_isoT_A=1_max_all_nominal_model.npz')
T_out = {key: T_out_loaded[key] for key in T_out_loaded}

# Load field data
field_data = np.load('final_new_nominal_model_field_data_nominal_model_4-10_isoT_A=1_max_all_nominal_model.npz')
x_out = field_data['x_out']
y_out = field_data['y_out']
z_out = field_data['z_out']
B_out = field_data['B_out']

r_out = np.sqrt(x_out**2 + y_out**2 + z_out**2)
rho_out = np.sqrt(x_out**2 + y_out**2)
# lat_out_deg = np.degrees(np.arcsin(z_out / r_out))  # (not directly used here)

# Flatten coordinate arrays
r_flat   = r_out.flatten()
rho_flat = rho_out.flatten()
z_flat   = z_out.flatten()

# -------------------------------------------------------------------------
# 2) DEFINE YOUR NEW QUANTITIES
#    netotal = n_elec + n_eh
#    feh     = n_eh / (n_elec + n_eh)
# -------------------------------------------------------------------------
n_elec_flat = n_out['elec'].flatten()
n_eh_flat   = n_out['eh'].flatten()

Tpar_f  = T_out["elec_par"].flatten()   # shape (601, n_lat)
Tperp_f = T_out["elec_perp"].flatten()
A_f_ = Tperp_f[Tpar_f>0]/Tpar_f[Tpar_f>0]
A_f = A_f_[A_f_ > 0]

netotal_flat = n_elec_flat + n_eh_flat
# Avoid any divide-by-zero
feh_flat = np.zeros_like(netotal_flat)
valid_div = (netotal_flat > 0)
feh_flat[valid_div] = n_eh_flat[valid_div] / netotal_flat[valid_div]


print('netot max and min, 2.5, 97.5 percentile =',np.nanmax(netotal_flat[netotal_flat>0]),np.nanmin(netotal_flat[netotal_flat>0]), np.nanpercentile(netotal_flat[netotal_flat>0], 2.5), np.nanpercentile(netotal_flat[netotal_flat>0], 97.5) )

print('feh max and min, 2.5, 97.5 percentile =',np.nanmax(feh_flat[feh_flat>0]),np.nanmin(feh_flat[feh_flat>0]), np.nanpercentile(feh_flat[feh_flat>0], 2.5), np.nanpercentile(feh_flat[feh_flat>0], 97.5) )

print('Tpar_ec max and min, 2.5, 97.5 percentile =',np.nanmax(Tpar_f[Tpar_f>0]),np.nanmin(Tpar_f[Tpar_f>0]), np.nanpercentile(Tpar_f[Tpar_f>0], 2.5), np.nanpercentile(Tpar_f[Tpar_f>0], 97.5) )
print('Tperp_ec max and min,2.5, 97.5percentile =',np.nanmax(Tperp_f[Tperp_f>0]),np.nanmin(Tperp_f[Tperp_f>0]), np.nanpercentile(Tperp_f[Tperp_f>0], 2.5), np.nanpercentile(Tperp_f[Tperp_f>0], 97.5))

print('A_ec max and min, 2.5, 97.5 percentile =',np.nanmax(A_f),np.nanmin(A_f), np.nanpercentile(A_f, 2.5), np.nanpercentile(A_f, 97.5))


Tpar_f  = T_out["eh_par"].flatten()   # shape (601, n_lat)
Tperp_f = T_out["eh_perp"].flatten()
A_f_ = Tperp_f[Tpar_f>0]/Tpar_f[Tpar_f>0]
A_f = A_f_[A_f_ > 0]

print('Tpar_eh max and min, 2.5, 97.5 percentile =',np.nanmax(Tpar_f[Tpar_f>1]),np.nanmin(Tpar_f[Tpar_f>1]), np.nanpercentile(Tpar_f[Tpar_f>1], 2.5), np.nanpercentile(Tpar_f[Tpar_f>1], 97.5) )
print('Tperp_eh max and min,2.5, 97.5percentile =',np.nanmax(Tperp_f[Tperp_f>1]),np.nanmin(Tperp_f[Tperp_f>1]), np.nanpercentile(Tperp_f[Tperp_f>1], 2.5), np.nanpercentile(Tperp_f[Tperp_f>1], 97.5))

print('A_eh max and min, 2.5, 97.5 percentile =',np.nanmax(A_f),np.nanmin(A_f), np.nanpercentile(A_f, 2.5), np.nanpercentile(A_f, 97.5))

# -------------------------------------------------------------------------
# 3) JUPITER IMAGE (SAME AS ORIGINAL)
# -------------------------------------------------------------------------
jupiter_img = mpimg.imread('hubble-captures-vivid-auroras-in-jupiters-atmosphere_28000029525_o~large.png')
# Add alpha channel if missing
if jupiter_img.shape[2] == 3:
    alpha_channel = np.ones((jupiter_img.shape[0], jupiter_img.shape[1], 1))
    jupiter_img = np.concatenate((jupiter_img, alpha_channel), axis=2)

# Jupiter oblateness factor
oblateness_factor = 0.93512
img_extent = [-1.0, 1.0, -oblateness_factor, oblateness_factor]

# -------------------------------------------------------------------------
# 4) CREATE FIGURE AND 2 SUBPLOTS (INSTEAD OF 10)
# -------------------------------------------------------------------------
fig = plt.figure(figsize=(9, 10))
nrows = 2
ncols = 1
gs = fig.add_gridspec(nrows, ncols, left=0.075, bottom=0.05, hspace=0.0, wspace=0.0)
(ax1, ax2) = gs.subplots(sharex='col', sharey='row')  # two axes: top(ax1), bottom(ax2)

axes = [ax1, ax2]

# We'll define some shared parameters from your original code
nlevels = 1000
x_min, x_max = -1.5, 12.0
y_min, y_max = -3.0, 3.0

# -------------------------------------------------------------------------
# For netotal, use the SAME color scale as your original "n_elec" code
# -------------------------------------------------------------------------
minnec = 1.0  # same as you used for "elec"
maxnec = np.nanmax(netotal_flat) #n_out_iso_max['elec'].max() #+ n_out_iso_max['eh'].max()  # Keep same upper limit as the cold-electron max
levelsnec = np.logspace(np.log10(minnec), np.log10(maxnec), nlevels)
normnec = colors.LogNorm(vmin=minnec, vmax=maxnec)

# Use black contour lines you previously used for n_elec
contour_levels_netotal = [10, 100, 1000]  # from your original 'elec' lines

# Ticks and tick labels for the colorbar (same as you had for n_elec)
cbar_ticks_nec = [
    1, 2, 3, 4, 5, 6, 7, 8, 9,
    10, 20, 30, 40, 50, 60, 70, 80, 90,
    100, 200, 300, 400, 500, 600, 700, 800, 900,
    1000, 2000, 3000
]
cbar_label_ticks_nec = [1, 10, 100, 1000, 3000]
cbar_tick_labels_nec = {tick: str(tick) if tick in cbar_label_ticks_nec else '' for tick in cbar_ticks_nec}

# -------------------------------------------------------------------------
# For feh, new color scale: [5e-4 to 0.5], black contours at
# [0.0005, 0.001, 0.005, 0.01, 0.1, 0.5]
# -------------------------------------------------------------------------
minfeh = np.nanmin(feh_flat[feh_flat > 1e-30])#5e-4
maxfeh = np.nanmax(feh_flat)#0.5
print('maxfeh = ',maxfeh)
levelsfeh = np.logspace(np.log10(minfeh), np.log10(maxfeh), nlevels)
normfeh = colors.LogNorm(vmin=minfeh, vmax=maxfeh)

contour_levels_feh = [0.0025,0.0035, 0.005, 0.01, 0.1, 0.5]

# We can define colorbar ticks for feh similarly
cbar_ticks_feh = [ 0.0015, 0.005, 0.01, 0.1, 0.5]
cbar_tick_labels_feh = {tick: str(tick) for tick in cbar_ticks_feh}

# -------------------------------------------------------------------------
# HELPER FUNCTION FOR GRID & MASK (same as your approach)
# -------------------------------------------------------------------------
def get_density_grid(density_flat):
    """
    1) Only take valid = your original radial constraints
    2) Interpolate via griddata
    3) Zero out anything outside the region
    """
    valid = (
        np.isfinite(density_flat) & (density_flat > 0)
        & (r_flat > 1.6)
        & (z_flat >= -3) & (z_flat <= 3)
        & (rho_flat > 4.4*(rho_flat**3)/((z_flat**2 + rho_flat**2)**1.5))
        & (rho_flat < 9.9*(rho_flat**3)/((z_flat**2 + rho_flat**2)**1.5))
    )
    rho_val = rho_flat[valid]
    z_val   = z_flat[valid]
    d_val   = density_flat[valid]

    # Define same large grid in rho,z
    num_rho_points = 12000
    num_z_points   = 600
    rho_lin = np.linspace(0.0, 12.0, num_rho_points)
    z_lin   = np.linspace(-3.0, 3.0, num_z_points)
    RHO_grid, Z_grid = np.meshgrid(rho_lin, z_lin)
    R_grid = np.sqrt(RHO_grid**2 + Z_grid**2)

    # Interpolate
    points = np.column_stack((rho_val, z_val))
    density_grid = griddata(points, d_val, (np.abs(RHO_grid), Z_grid), method='linear')
    density_grid = np.nan_to_num(density_grid, nan=0.0)

    # Mask outside region
    region_mask = (
        (R_grid > 1.6)
        & (Z_grid >= -3) & (Z_grid <= 3)
        & (np.abs(RHO_grid) > 4.4*(np.abs(RHO_grid)**3)/((Z_grid**2 + RHO_grid**2)**1.5))
        & (np.abs(RHO_grid) < 9.9*(np.abs(RHO_grid)**3)/((Z_grid**2 + RHO_grid**2)**1.5))
    )
    density_grid[~region_mask] = 0.0

    return RHO_grid, Z_grid, density_grid

# -------------------------------------------------------------------------
# 5) PLOT 1: n_total = n_elec + n_eh (top subplot, ax1)
# -------------------------------------------------------------------------
rho_grid_netotal, z_grid_netotal, netotal_grid = get_density_grid(netotal_flat)

# Contourf
contour_netotal = ax1.contourf(
    rho_grid_netotal, z_grid_netotal, netotal_grid,
    levels=levelsnec, norm=normnec, cmap='viridis'
)
for c in contour_netotal.collections:
    c.set_rasterized(True)

# Black contour lines
clines_netotal = ax1.contour(
    rho_grid_netotal, z_grid_netotal, netotal_grid,
    levels=contour_levels_netotal, colors='black', linewidths=1
)
ax1.clabel(clines_netotal, fmt='%g', inline=True, fontsize=8, colors='black')

# Axes limits, aspect
ax1.set_xlim([x_min, x_max])
ax1.set_ylim([y_min, y_max])
ax1.set_aspect('equal')

# Turn on ticks
ax1.tick_params(axis='both', which='both', direction='in', length=4)
# For top subplot, you can hide the x-labels if desired
ax1.tick_params(labelbottom=False)
ax1.set_yticks([-3, -2, -1, 0, 1, 2])

# Label inside the plot
ax1.text(3.0, 0.0, r'$n_{\mathrm{total}}$', fontsize=12, fontweight='bold',
         ha='center', va='center', color='black')
ax1.text(3.0, -0.7, '(cm$^{-3}$)', fontsize=12, fontweight='bold',
         ha='center', va='center', color='black')

# Overlay Jupiter image
ax1.imshow(jupiter_img, extent=img_extent, aspect='auto', origin='upper', zorder=10)

# Inset colorbar for n_total
cax1 = ax1.inset_axes([0.85, 0.3, 0.03, 0.4])  # [x0, y0, width, height] in axes fraction
cb1 = fig.colorbar(contour_netotal, cax=cax1, orientation='vertical')
cb1.set_ticks(cbar_ticks_nec)
cb1.ax.set_yticklabels([cbar_tick_labels_nec.get(tick, '') for tick in cbar_ticks_nec])
for label in cb1.ax.get_yticklabels():
    label.set_fontweight('bold')
cb1.ax.tick_params(labelsize=9)
cb1.set_label('', labelpad=1)

# -------------------------------------------------------------------------
# 6) PLOT 2: f_ehot = n_eh / (n_elec + n_eh) (bottom subplot, ax2)
# -------------------------------------------------------------------------
rho_grid_feh, z_grid_feh, feh_grid = get_density_grid(feh_flat)

# Contourf
contour_feh = ax2.contourf(
    rho_grid_feh, z_grid_feh, feh_grid,
    levels=levelsfeh, norm=normfeh, cmap='viridis'
)
for c in contour_feh.collections:
    c.set_rasterized(True)

# Black contour lines
clines_feh = ax2.contour(
    rho_grid_feh, z_grid_feh, feh_grid,
    levels=contour_levels_feh, colors='black', linewidths=1
)
ax2.clabel(clines_feh, fmt='%g', inline=True, fontsize=8, colors='black')

# Axes limits, aspect
ax2.set_xlim([x_min, x_max])
ax2.set_ylim([y_min, y_max])
ax2.set_aspect('equal')

# Tick params
ax2.tick_params(axis='both', which='both', direction='in', length=4)
ax2.set_yticks([-3, -2, -1, 0, 1, 2])
ax2.set_xticks([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12])
# Label inside the plot
ax2.text(3.0, 0.0, r'$f_{\mathrm{e(hot)}}$', fontsize=12, fontweight='bold',
         ha='center', va='center', color='black')
ax2.text(3.0, -0.7, r'$\frac{n_{e^{-}(\text{hot})}}{n_{e^{-}}+n_{e^{-}(\text{hot})}}$', 
         fontsize=10, fontweight='bold', ha='center', va='center', color='black')

# Overlay Jupiter
ax2.imshow(jupiter_img, extent=img_extent, aspect='auto', origin='upper', zorder=10)

# Inset colorbar for f_ehot
cax2 = ax2.inset_axes([0.85, 0.3, 0.03, 0.4])
cb2 = fig.colorbar(contour_feh, cax=cax2, orientation='vertical')
cb2.set_ticks(cbar_ticks_feh)
cb2.ax.set_yticklabels([cbar_tick_labels_feh.get(tick, '') for tick in cbar_ticks_feh])
for label in cb2.ax.get_yticklabels():
    label.set_fontweight('bold')
cb2.ax.tick_params(labelsize=9)
cb2.set_label('', labelpad=1)

# -------------------------------------------------------------------------
# 7) ADD COMMON X AND Y LABELS
# -------------------------------------------------------------------------
fig.text(0.49, 0.02, r'$\rho_c$ ($R_J$)', ha='center', fontsize=16, fontweight='bold')
fig.text(0.02, 0.5, r'$z_c$ ($R_J$)', va='center', rotation='vertical', fontsize=16, fontweight='bold')

# -------------------------------------------------------------------------
# 8) SAVE AND SHOW
# -------------------------------------------------------------------------
plt.savefig('new_caseA_2subplots_netotal_feh.pdf', dpi=600, bbox_inches='tight', pad_inches=0)
plt.show()
