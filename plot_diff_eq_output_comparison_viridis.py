#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 13:15:18 2024

@author: edne8319
"""

#import matplotlib
#import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib.colors import LogNorm

# Load x, y, z values
x = np.transpose(np.loadtxt('x_out_RJ_lat-70to70_1401latvals_only9.38_integralproperly_aligned_phi0=75.85_perfectaligned_for_4.70-9.71.txt', delimiter=','))
y = np.transpose(np.loadtxt('y_out_RJ_lat-70to70_1401latvals_only9.38_integralproperly_aligned_phi0=75.85_perfectaligned_for_4.70-9.71.txt', delimiter=','))
z = np.transpose(np.loadtxt('z_out_RJ_lat-70to70_1401latvals_only9.38_integralproperly_aligned_phi0=75.85_perfectaligned_for_4.70-9.71.txt', delimiter=','))


# Load densities
iso_kappas_nec = np.transpose(np.loadtxt('aniso_kappas_mv_firast_trynelec_out_mymodel1_diffeq_jrm33+Con2020_0fillslat-70to70_601_3D_integral_including_feh_properly_aligned_phi0=75.85_perfectaligned_for_4.70-9.71.txt', delimiter=','))
aniso_kappas_nec = np.transpose(np.loadtxt('true_aniso_kappas_mv_firast_trynelec_out_mymodel1_diffeq_jrm33+Con2020_0fillslat-70to70_601_3D_integral_including_feh_properly_aligned_phi0=75.85_perfectaligned_for_4.70-9.71.txt', delimiter=','))
fried_egg_nec = np.transpose(np.loadtxt('fried_agg_nelec_out_mymodel1_diffeq_jrm33+Con2020_0fillslat-70to70_601_3D_integral_including_feh_properly_aligned_phi0=75.85_perfectaligned_for_4.70-9.71.txt', delimiter=','))
nec =np.transpose( np.loadtxt('nelec_out_mymodel1_diffeq_jrm33+Con2020_0fillslat-70to70_601_3D_integral_including_feh_properly_aligned_phi0=75.85_perfectaligned_for_4.70-9.71.txt', delimiter=','))


x = x.reshape(842001)
y = y.reshape(842001)
z = z.reshape(842001)
iso_kappas_nec=iso_kappas_nec.reshape(842001)
aniso_kappas_nec = aniso_kappas_nec.reshape(842001)
fried_egg_nec = fried_egg_nec.reshape(842001)
nec = nec.reshape(842001)

rho = np.sqrt(x ** 2. + y ** 2.)
r = np.sqrt(x ** 2. + y ** 2. + z ** 2.)


mask = r >= 1.01
x = x[mask]
y = y[mask]
z = z[mask]
r=r[mask]
rho = rho[mask]

iso_kappas_nec = iso_kappas_nec[mask]
aniso_kappas_nec=aniso_kappas_nec[mask]
fried_egg_nec = fried_egg_nec[mask]
nec = nec[mask]


print(fried_egg_nec.shape)
print(fried_egg_nec.size)


fig, axs = plt.subplots(4, 1, figsize=(10, 15))

triangog = tri.Triangulation(rho, z)

# Define the number of levels
num_levels = 20  # for example

# Generate the levels
#levels = np.linspace(0., 3500., num_levels)
mmin = 0.1
mmax = 3500.
levels = np.logspace(np.log10(mmin), np.log10(mmax), num=num_levels)
norm = LogNorm(vmin=mmin, vmax=mmax)  # Set these limits to fit your data range
# Plot aniso_kappas
cs1 = axs[0].tricontourf(triangog,iso_kappas_nec , levels = levels,cmap='viridis', norm=norm)
fig.colorbar(cs1, ax=axs[0])
axs[0].set_title(r'Isotropic Kappas n$_{ec}$ cm$^{-3}$')
axs[0].set_xlabel(r'$\rho$ ($R_J$)')
axs[0].set_ylabel(r'z ($R_J$)')
axs[0].set_xlim([-1.5,12])
axs[0].set_ylim([-5,5])
#axs[0].set_zlim([0,3500])


cs2 = axs[1].tricontourf(triangog, aniso_kappas_nec, levels = levels, cmap='viridis', norm=norm)
fig.colorbar(cs2, ax=axs[1])
axs[1].set_title(r'Anisotropic Kappas n$_{ec}$ cm$^{-3}$')
axs[1].set_xlabel(r'$\rho$ ($R_J$)')
axs[1].set_ylabel(r'z ($R_J$)')
axs[1].set_xlim([-1.5,12])
axs[1].set_ylim([-5,5])
#axs[1].set_zlim([0,3500])

mask = np.isfinite(fried_egg_nec)

rho_clean = rho[mask]
z_clean = z[mask]
fried_egg_nec_clean = fried_egg_nec[mask]
 
triang = tri.Triangulation(rho_clean ,z_clean  )

"""

# Plot fried_agg_nelec
cs3 = axs[2].tricontourf(triang, fried_egg_nec_clean, levels = levels, cmap='viridis')
fig.colorbar(cs3, ax=axs[2])
axs[2].set_title(r'Fried-Egg n$_{ec}$ cm$^{-3}$' )
axs[2].set_xlabel(r'$\rho$ ($R_J$)')
axs[2].set_ylabel(r'z ($R_J$)')
axs[2].set_xlim([-1.5,12])
axs[2].set_ylim([-5,5])
#axs[2].set_zlim([0,3500])
"""

# Plot fried_agg_nelec
cs4 = axs[2].tricontourf(triangog, nec, levels = levels, cmap='viridis', norm=norm)
fig.colorbar(cs4, ax=axs[2])
axs[2].set_title(r'Iso-Max n$_{ec}$ cm$^{-3}$' )
axs[2].set_xlabel(r'$\rho$ ($R_J$)')
axs[2].set_ylabel(r'z ($R_J$)')
axs[2].set_xlim([-1.5,12])
axs[2].set_ylim([-5,5])
#axs[2].set_zlim([0,3500])

plt.tight_layout()
plt.savefig('diffeq_plot_comparison.png', dpi =600)
plt.show()
