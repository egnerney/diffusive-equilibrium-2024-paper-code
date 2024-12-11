import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import pi
from scipy.special import gamma
import matplotlib.gridspec as gridspec
from matplotlib.font_manager import FontProperties
#from matplotlib.ticker import ScalarFormatter

from matplotlib.colors import LogNorm


#plt.rcParams['text.usetex'] = True
#plt.rcParams['text.latex.preamble'] = [r'\usepackage{bm}']

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
#v_perp = np.linspace(0, 2*val, vsamples) * 1e3
v_perp = np.linspace(-val, val, vsamples) * 1e3
V_perp, V_parallel = np.meshgrid( v_perp, v_parallel)

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


fs =16
fseqn = 10

#plt.rcParams.update({
#    "text.usetex": True,
#    "text.latex.preamble": r"\usepackage{amsmath, amsfonts,float,gensymb,graphicx}",
#   "font.family": "DejaVu Sans"})

#plt.rcParams['mathtext.fontset'] = 'stix'

unitconv=1.#1e12
nlevels = 1000
vvmin = min([fanisomax.min()/unitconv,ffriedegg.min()/unitconv,fanisokappa.min()/unitconv,fanisokappa2.min()/unitconv])
vvmax = max([fanisomax.max()/unitconv,ffriedegg.max()/unitconv,fanisokappa.max()/unitconv,fanisokappa2.max()/unitconv])
print(vvmin,vvmax)
vvmin = 0.05e-17*(1e12/unitconv)
vvmax = 2e-17*(1e12/unitconv)


# Use logarithmically spaced levels
levels = np.logspace(np.log10(vvmin), np.log10(vvmax), num=nlevels)
# Define LogNorm normalization
norm = LogNorm(vmin=vvmin, vmax=vvmax)  # Set these limits to fit your data range


xlabel='Plasma Density [cm$^{-3}$]'

# Function calculations and plotting
titles = ['Aniso-Bi-Maxwellian', 'Aniso-Bi-Kappa', 'Aniso-Bi-Kappa Product', '"Fried Egg"']

#fanisomaxstr = r'$\bf{\dfrac{m^{3/2} n}{2 \sqrt{2} \pi ^{3/2} \sqrt{\text{T}_{\parallel}} \text{T}_{\perp}} \text{e}^{\left(-\dfrac{m \text{v}_{\parallel}^2}{2 \text{T}_{\parallel}}\right)}  \text{e}^{\left(-\dfrac{m \text{v}_{\perp}^2}{2 \text{T}_{\perp}}\right)}}$'
#fanisokappastr = r'$\bf{ \dfrac{m^{3/2} n \Gamma (\kappa )}{2 \sqrt{2} \pi ^{3/2} \text{T}_{\perp} \Gamma \left(\kappa -\dfrac{1}{2}\right) \sqrt{\kappa \text{T}_{ \parallel}}}  \left(1 + \dfrac{mv_{ \perp}^2}{2 \kappa T_{ \perp}} + \dfrac{mv_{\ \parallel}^2}{2 \kappa T_{ \parallel}}\right)^{-\left(1+\kappa\right)}}$'
#fanisokappa2str = r'$\bf{ \dfrac{m^{3/2} n\Gamma (\kappa +1)}{2 \sqrt{2} \pi ^{3/2} \text{T}_{\perp} \Gamma \left(\kappa +\dfrac{1}{2}\right) \sqrt{\kappa  \text{T}_{\parallel}}} \left(1 + \dfrac{m v_{\parallel}^2}{2 \kappa T_{ \parallel}}\right)^{-\left(1+\kappa\right)} \left(1 + \dfrac{m v_{\perp}^2}{2 \kappa  T_{\perp}}\right)^{-\left(1+\kappa\right)} }$'
#ffriedeggstr = r'$ \bf{\dfrac{m^{3/2} n}{2 \sqrt{2} \pi ^{3/2} \sqrt{\text{T}_{\parallel}} \text{T}_{\perp}} \text{e}^{-\dfrac{m v^2_{ \parallel }}{2  T_{ \parallel}}}\left(1 + \dfrac{m v_{\perp}^2}{2 \kappa  T_{\perp}}\right)^{-\left(1+\kappa\right)}} $'
fanisomaxstr = r'$\dfrac{\bf{m^{3/2} n}}{\bf{2^{3/2} \pi^{3/2} \sqrt{\text{T}_{\parallel}}} \bf{\text{T}}_{\bot}} \bf{\text{e}^{-\dfrac{m \text{v}_{\parallel}^2}{2 \text{T}_{\parallel}}}} \bf{\text{e}}^{\bf{-}\dfrac{\bf{m \text{v}}_{\bot}^{\bf{2}}}{\bf{2 \text{T}}_{\bot}}}$'
fanisokappastr = r'$\dfrac{\bf{m^{3/2} n \Gamma (\kappa )}}{\bf{2 \sqrt{2} \pi ^{3/2} \text{T}}_{\perp} \bf{\Gamma \left(\kappa -\dfrac{1}{2}\right) \sqrt{\kappa \text{T}_{ \parallel}}}} \left(\bf{1} + \dfrac{\bf{mv}_{\perp}\bf{^2}}{\bf{2 \kappa \text{T}}_{\perp}} + \bf{\dfrac{mv_{\parallel}^2}{2 \kappa \text{T}_{\parallel}}}\right)^{\bf{-\left(1+\kappa\right)}}$'
fanisokappa2str = r'$\dfrac{\bf{m^{3/2} n \Gamma (\kappa +1)}}{\bf{2 \sqrt{2} \pi ^{3/2} \text{T}}_{\perp} \bf{\Gamma \left(\kappa +\dfrac{1}{2}\right) \sqrt{\kappa \text{T}_{\parallel}}}} \bf{\left(1 + \dfrac{m v_{\parallel}^2}{2 \kappa \text{T}_{\parallel}}\right)^{-\left(1+\kappa\right)}} \left(\bf{1} + \dfrac{\bf{m v}_{\perp}\bf{^2}}{\bf{2 \kappa \text{T}}_{\perp}}\right)^{\bf{-\left(1+\kappa\right)}}$'
ffriedeggstr = r'$\dfrac{\bf{m^{3/2} n}}{\bf{2 \sqrt{2} \pi ^{3/2} \sqrt{\bf{\text{T}_{\parallel}}} \text{T}}_{\perp}} \bf{\text{e}^{-\dfrac{m v^2_{\parallel}}{2 \text{T}_{\parallel}}}} \left(\bf{1} + \dfrac{\bf{m v}_{\perp}^2}{\bf{2 \kappa \text{T}}_{\perp}}\right)^{\bf{-\left(1+\kappa\right)}}$'


eqns = [fanisomaxstr,fanisokappastr,fanisokappa2str,ffriedeggstr]
functions = [fanisomax, fanisokappa, fanisokappa2, ffriedegg]
for i, function in enumerate(functions):
    ax = fig.add_subplot(gs[i // 2, i % 2])
    cp = ax.contourf( V_perp * 1e-3,V_parallel * 1e-3, function / unitconv, levels=levels, cmap='viridis', norm=norm)
    ax.set_aspect('equal')
    ax.text(0.5, 0.92, titles[i], color='black', ha='center', va='bottom', transform=ax.transAxes, fontsize=fs, weight='bold')
    ax.text(0.5, 0.04, eqns[i], color='black', ha='center', va='bottom', transform=ax.transAxes, fontsize=fseqn, weight='bold')
    if i // 2 == 0:  # Top row
        ax.set_xticklabels([])
    if i % 2 == 1:  # Right column
        ax.set_yticklabels([])
    if i >= 2:  # Bottom row
        ax.set_xlabel('v$_{\perp}$ (km/s)', fontproperties=bold_font,fontsize=fs)
    if i % 2 == 0:  # Left column
        ax.set_ylabel(r'$\bf{v}_{\bf{\parallel}}$ (km/s)', fontproperties=bold_font,fontsize=fs)


    """
    # Special handling for i=2 to adjust tick labels
    if i == 2:
        ax.set_yticks(ax.get_xticks().tolist())  # Refresh tick setting
        ax.set_yticklabels([f"{y:.0f}" if y != 40 else "" for y in ax.get_yticks()])
        ax.set_xticks(ax.get_yticks().tolist())  # Refresh tick setting
        ax.set_xticklabels([f"{x:.0f}" if x != 80 else "" for x in ax.get_xticks()])
    """
   

    # Set tick labels to bold
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontproperties(bold_font)
        label.set_fontsize(fs)
    ax.set_aspect('equal')
# Adjust subplot spacing and add colorbar
fig.subplots_adjust(left=0.05, right=0.8, top=0.95, bottom=0.05, wspace=0, hspace=0)
cbar_ax = fig.add_axes([0.85, 0.15, 0.03, 0.7])
cbar = fig.colorbar(cp, cax=cbar_ax, format='%.2e')
cbar.set_label(r'm$^{\bf{-6}}$ s$^{\bf{3}}$',  weight='bold',fontsize=16)
#cbar.ax.tick_params(labelsize=10)  # Tick marks on colorbar
for l in cbar.ax.yaxis.get_ticklabels():
    l.set_weight("bold")
    l.set_fontsize(fs)
    
for l in cbar.ax.xaxis.get_ticklabels():
    l.set_weight("bold")
    l.set_fontsize(fs)
 
# Attempt to directly modify the exponent after rendering
def format_exponent(ax, fontweight='bold', fontsize=fs):
    """Format the exponent text to bold and set the fontsize."""
    ax.yaxis.get_offset_text().set_fontsize(fontsize)
    ax.yaxis.get_offset_text().set_weight(fontweight)
    #x.yaxis.get_offset_text().set_style('italic')  # You can also set to 'normal' if italics are not desired

# Apply custom formatting to the exponent
format_exponent(cbar.ax)

# Adding text box above the colorbar
textstr = r'Species: O$^{+}$' + '\n' + \
          r'n = 1750 cm$\bf{^{-3}}$' + '\n' + \
          r'T$_{\bot}$ = 70 eV' + '\n' + \
          r'T$\bf{_{\parallel}}$ = 35 eV' + '\n' + \
          r'$\bf{\kappa}$ = 2'
              
fig.text(0.9, 0.97, textstr, fontsize=fs,weight='bold', verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))



fig.set_size_inches(12, 10)


plt.savefig(fname='log_fmin=0.1e-17_fmax=2e-17_diff_Op_aniso_dists_2D_1000vlevels_1000vsamples_f_in_mks.pdf',dpi=600, bbox_inches='tight')
plt.show()