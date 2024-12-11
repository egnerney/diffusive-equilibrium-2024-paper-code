# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 11:33:11 2024

@author: Owner
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import pi
from scipy.special import gamma
import matplotlib.gridspec as gridspec
from matplotlib.font_manager import FontProperties
from matplotlib.colors import LogNorm

k = 1.60218e-19  # Boltzmann constant in Joules/eV
m_O = 16 * 1.66054e-27  # mass of Oxygen in kg (16 times the atomic mass unit)

# Parameters
n_alpha = 1750 * 1e6  # converting from cm^-3 to m^-3
T_parallel = 35 * k  # converting from eV to Joules
T_perp = 70 * k  # converting from eV to Joules
kappa_alpha = 2
val = 65.
vsamples = 1000
# Velocity grid in m/s
v_parallel = np.linspace(-val, val, vsamples) * 1e3
v_perp = np.linspace(-val, val, vsamples) * 1e3
V_perp, V_parallel = np.meshgrid(v_perp, v_parallel)

# Create the figure and specify GridSpec
fig = plt.figure(figsize=(12, 10))
gs = gridspec.GridSpec(2, 2, figure=fig)

# Define bold font properties
bold_font = FontProperties(weight='bold')

# Function calculations
ffriedegg = (m_O**1.5 * n_alpha / (2 * np.sqrt(2) * pi**(3/2) * np.sqrt(T_parallel) * T_perp)) * \
    np.exp(-m_O * V_parallel**2 / (2 * T_parallel)) * \
    (1 + m_O * V_perp**2 / (2 * kappa_alpha * T_perp))**(-(1 + kappa_alpha))

fanisokappa = (m_O**1.5 * n_alpha * gamma(kappa_alpha) / 
     (2 * np.sqrt(2) * pi**(3/2) * T_perp * gamma(kappa_alpha - 0.5) * np.sqrt(kappa_alpha * T_parallel))) * \
    (1 + m_O * V_perp**2 / (2 * kappa_alpha * T_perp) + m_O * V_parallel**2 / (2 * kappa_alpha * T_parallel))**(-(1 + kappa_alpha))

fanisokappa2 = (m_O**1.5 * n_alpha * gamma(kappa_alpha + 1) / (2 * np.sqrt(2) * pi**(3/2) * T_perp * gamma(kappa_alpha + 0.5) * np.sqrt(kappa_alpha * T_parallel))) * \
       (1 + m_O * V_parallel**2 / (2 * kappa_alpha * T_parallel))**(-(1 + kappa_alpha)) * \
       (1 + m_O * V_perp**2 / (2 * kappa_alpha * T_perp))**(-(1 + kappa_alpha))

fanisomax = (m_O**1.5 * n_alpha / (2 * np.sqrt(2) * pi**1.5 * np.sqrt(T_parallel) * T_perp)) * \
    np.exp(-m_O * V_parallel**2 / (2 * T_parallel)) * \
    np.exp(-m_O * V_perp**2 / (2 * T_perp))

fs = 16
fseqn = 10

unitconv = 1.  # 1e12
nlevels = 1000
vvmin = min([fanisomax.min()/unitconv, ffriedegg.min()/unitconv, fanisokappa.min()/unitconv, fanisokappa2.min()/unitconv])
vvmax = max([fanisomax.max()/unitconv, ffriedegg.max()/unitconv, fanisokappa.max()/unitconv, fanisokappa2.max()/unitconv])
print(vvmin, vvmax)
vvmin = 0.05e-17*(1e12/unitconv)
vvmax = 2e-17*(1e12/unitconv)

# Use logarithmically spaced levels
levels = np.logspace(np.log10(vvmin), np.log10(vvmax), num=nlevels)
# Define LogNorm normalization
norm = LogNorm(vmin=vvmin, vmax=vvmax)  # Set these limits to fit your data range

xlabel = 'Plasma Density [cm$^{-3}$]'

# Function calculations and plotting
titles = ['Maxwellian', 'Standard-Kappa', 'Product-Kappa', '"Fried Egg"']

# Updated equations with \perp outside bold formatting
fanisomaxstr = r'$\dfrac{\boldsymbol{m}^{3/2} \boldsymbol{n}}{2^{3/2} \pi^{3/2} \sqrt{\boldsymbol{\text{T}}_{\parallel}} \boldsymbol{\text{T}}_{\perp}} \boldsymbol{e}^{-\dfrac{\boldsymbol{m} \boldsymbol{\text{v}}_{\parallel}^2}{2 \boldsymbol{\text{T}}_{\parallel}}} \boldsymbol{e}^{-\dfrac{\boldsymbol{m} \boldsymbol{\text{v}}_{\perp}^2}{2 \boldsymbol{\text{T}}_{\perp}}}$'

fanisokappastr = r'$\dfrac{\boldsymbol{m}^{3/2} \boldsymbol{n} \Gamma (\boldsymbol{\kappa})}{2 \sqrt{2} \pi^{3/2} \boldsymbol{\text{T}}_{\perp} \Gamma \left(\boldsymbol{\kappa} - \dfrac{1}{2}\right) \sqrt{\boldsymbol{\kappa} \boldsymbol{\text{T}}_{\parallel}}} \left( \boldsymbol{1} + \dfrac{\boldsymbol{m} \boldsymbol{\text{v}}_{\perp}^2}{2 \boldsymbol{\kappa} \boldsymbol{\text{T}}_{\perp}} + \dfrac{\boldsymbol{m} \boldsymbol{\text{v}}_{\parallel}^2}{2 \boldsymbol{\kappa} \boldsymbol{\text{T}}_{\parallel}} \right)^{ - \left( 1 + \boldsymbol{\kappa} \right)}$'

fanisokappa2str = r'$\dfrac{\boldsymbol{m}^{3/2} \boldsymbol{n} \Gamma (\boldsymbol{\kappa} + 1)}{2 \sqrt{2} \pi^{3/2} \boldsymbol{\text{T}}_{\perp} \Gamma \left(\boldsymbol{\kappa} + \dfrac{1}{2}\right) \sqrt{\boldsymbol{\kappa} \boldsymbol{\text{T}}_{\parallel}}} \left( \boldsymbol{1} + \dfrac{\boldsymbol{m} \boldsymbol{\text{v}}_{\parallel}^2}{2 \boldsymbol{\kappa} \boldsymbol{\text{T}}_{\parallel}} \right)^{ - \left( 1 + \boldsymbol{\kappa} \right)} \left( \boldsymbol{1} + \dfrac{\boldsymbol{m} \boldsymbol{\text{v}}_{\perp}^2}{2 \boldsymbol{\kappa} \boldsymbol{\text{T}}_{\perp}} \right)^{ - \left( 1 + \boldsymbol{\kappa} \right)}$'

ffriedeggstr = r'$\dfrac{\boldsymbol{m}^{3/2} \boldsymbol{n}}{2 \sqrt{2} \pi^{3/2} \sqrt{ \boldsymbol{\text{T}}_{\parallel}} \boldsymbol{\text{T}}_{\perp}} \boldsymbol{e}^{ - \dfrac{ \boldsymbol{m} \boldsymbol{\text{v}}_{\parallel}^{2} }{2 \boldsymbol{\text{T}}_{\parallel} } } \left( \boldsymbol{1} + \dfrac{ \boldsymbol{m} \boldsymbol{\text{v}}_{\perp}^{2} }{2 \boldsymbol{\kappa} \boldsymbol{\text{T}}_{\perp} } \right)^{ - \left( 1 + \boldsymbol{\kappa} \right)}$'

eqns = [fanisomaxstr, fanisokappastr, fanisokappa2str, ffriedeggstr]
functions = [fanisomax, fanisokappa, fanisokappa2, ffriedegg]
for i, function in enumerate(functions):
    ax = fig.add_subplot(gs[i // 2, i % 2])
    cp = ax.contourf(V_perp * 1e-3, V_parallel * 1e-3, function / unitconv, levels=levels, cmap='viridis', norm=norm)
    ax.set_aspect('equal')
    ax.text(0.5, 0.92, titles[i], color='black', ha='center', va='bottom', transform=ax.transAxes, fontsize=fs, weight='bold')
    ax.text(0.5, 0.04, eqns[i], color='black', ha='center', va='bottom', transform=ax.transAxes, fontsize=fseqn, weight='bold')
    if i // 2 == 0:  # Top row
        ax.set_xticklabels([])
    if i % 2 == 1:  # Right column
        ax.set_yticklabels([])
    if i >= 2:  # Bottom row
        ax.set_xlabel('v$_{\perp}$ (km/s)', fontproperties=bold_font, fontsize=fs)
    if i % 2 == 0:  # Left column
        ax.set_ylabel(r'$\boldsymbol{v}_{\parallel}$ (km/s)', fontproperties=bold_font, fontsize=fs)

    # Set tick labels to bold
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontproperties(bold_font)
        label.set_fontsize(fs)
    ax.set_aspect('equal')

# Adjust subplot spacing and add colorbar
fig.subplots_adjust(left=0.05, right=0.8, top=0.95, bottom=0.05, wspace=0, hspace=0)
cbar_ax = fig.add_axes([0.85, 0.15, 0.03, 0.7])
cbar = fig.colorbar(cp, cax=cbar_ax, format='%.1e')
cbar.set_label(r'm$^{\boldsymbol{-6}}$ s$^{\boldsymbol{3}}$', weight='bold', fontsize=16)

# Define ticks and tick labels
ticks = [5e-7, 6e-7, 7e-7, 8e-7, 9e-7, 1e-6, 2e-6, 3e-6, 4e-6, 5e-6, 6e-6, 7e-6, 8e-6, 9e-6, 1e-5, 2e-5]
tick_labels = []
for t in ticks:
    if t in [5e-7, 1e-6, 1e-5]:
        tick_labels.append('{:.1e}'.format(t))
    else:
        tick_labels.append('')

# Set ticks and tick labels
cbar.set_ticks(ticks)
cbar.set_ticklabels(tick_labels)

# Set tick label properties
for l in cbar.ax.yaxis.get_ticklabels():
    l.set_weight("bold")
    l.set_fontsize(fs)

# Adding text box above the colorbar
textstr = r'Species: O$^{+}$' + '\n' + \
          r'n = 1750 cm$\boldsymbol{^{-3}}$' + '\n' + \
          r'T$_{\perp}$ = 70 eV' + '\n' + \
          r'T$\boldsymbol{_{\parallel}}$ = 35 eV' + '\n' + \
          r'$\boldsymbol{\kappa}$ = 2'

fig.text(0.9, 0.97, textstr, fontsize=fs, weight='bold', verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))

fig.set_size_inches(12, 10)

plt.savefig(fname='Fig1.png', dpi=300, bbox_inches='tight')
plt.show()
