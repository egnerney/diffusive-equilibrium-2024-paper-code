# -*- coding: utf-8 -*-
"""
Created on Fri Apr 18 01:20:53 2025

@author: Owner
"""

#!/usr/bin/env python3
# -------------------------------------------------------------------------
# Contour maps of electron T‖ and T⊥ vs ρ–z   (Jovian magnetosphere model)
# -------------------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from pathlib import Path

# -------------------------------------------------------------------------
# file names for the three distribution models ----------------------------
# -------------------------------------------------------------------------


# final_new_nominal_model_n_out_nominal_model_4-10_isoT_isokappa_A=1_standard_kappa_all_nominal_model.npz
# final_new_nominal_model_n_out_nominal_model_4-10_isoT_isokappa_A=1_product_kappa_all_nominal_model.npz
# final_new_nominal_model_n_out_nominal_model_4-10_isoT_isokappa_A=1_fried_egg_all_nominal_model.npz

# final_new_nominal_model_T_out_nominal_model_4-10_isoT_isokappa_A=1_standard_kappa_all_nominal_model
# final_new_nominal_model_T_out_nominal_model_4-10_isoT_isokappa_A=1_product_kappa_all_nominal_model
# final_new_nominal_model_T_out_nominal_model_4-10_isoT_isokappa_A=1_fried_egg_all_nominal_model.npz

# final_new_nominal_model_field_data_nominal_model_4-10_isoT_isokappa_A=1_standard_kappa_all_nominal_model.npz
# final_new_nominal_model_field_data_nominal_model_4-10_isoT_isokappa_A=1_product_kappa_all_nominal_model.npz
# final_new_nominal_model_field_data_nominal_model_4-10_isoT_isokappa_A=1_fried_egg_all_nominal_model.npz

MODELS = {
    "Fried‑Egg": dict(
        Tfile="final_new_nominal_model_T_out_nominal_model_4-10_isoT_isokappa_A=1_fried_egg_all_nominal_model.npz",
        Ffile="final_new_nominal_model_field_data_nominal_model_4-10_isoT_isokappa_A=1_fried_egg_all_nominal_model.npz",
    ),
    "Standard‑Kappa": dict(
        Tfile="final_new_nominal_model_T_out_nominal_model_4-10_isoT_isokappa_A=1_standard_kappa_all_nominal_model.npz",
        Ffile="final_new_nominal_model_field_data_nominal_model_4-10_isoT_isokappa_A=1_standard_kappa_all_nominal_model.npz",
    ),
    "Product‑Kappa": dict(
        Tfile="final_new_nominal_model_T_out_nominal_model_4-10_isoT_isokappa_A=1_product_kappa_all_nominal_model.npz",
        Ffile="final_new_nominal_model_field_data_nominal_model_4-10_isoT_isokappa_A=1_product_kappa_all_nominal_model.npz",
    ),
    "Maxwellian": dict(
        Tfile="final_new_nominal_model_T_out_nominal_model_4-10_isoT_A=1_max_all_nominal_model.npz",
        Ffile="final_new_nominal_model_field_data_nominal_model_4-10_isoT_A=1_max_all_nominal_model.npz",
    )
}

# -------------------------------------------------------------------------
def log_levels(arr, n=100):
    """Return *n* log‑spaced contour levels covering [arr_min, arr_max]."""
    positives = arr[arr > 0]
    vmin, vmax = positives.min(), positives.max()
    v975 = np.nanpercentile( positives, 97.5)
    return np.logspace(np.log10(vmin), np.log10(v975), n)
    #return np.logspace(np.log10(vmin), np.log10(vmax), n)

def plot_field(rho, z, T, comp, model, levels):
    """Triangular contour plot of one temperature component."""
    # Keep only >0 temps for color mapping
    mask = T > 0
    tri = Triangulation(rho[mask], z[mask])

    fig, ax = plt.subplots(figsize=(7, 6))
    cf = ax.tricontourf(tri, T[mask], levels=levels, extend="both")
    cbar = fig.colorbar(cf)
    cbar.set_label("Temperature (eV)")
    ax.set_xlabel(r"$\rho\;(\mathrm{R_J})$")
    ax.set_ylabel(r"$z\;(\mathrm{R_J})$")
    ax.set_aspect("equal", adjustable="box")
    ax.set_title(f"{model}: electron {comp}")
    plt.tight_layout()

    # save figure
    fname_png = f"{model.replace(' ', '_')}_electron_{comp}.png"
    fig.savefig(fname_png, dpi=300)
    print(f"saved → {fname_png}")
    plt.show()

# -------------------------------------------------------------------------
# Main loop over the three models -----------------------------------------
# -------------------------------------------------------------------------
for model, files in MODELS.items():
    Tfile = Path(files["Tfile"])
    Ffile = Path(files["Ffile"])

    # 1. temperatures
    with np.load(Tfile) as tnpz:
        T_par  = tnpz["elec_par"]   # shape (601, n_lat)
        T_perp = tnpz["elec_perp"]

    # 2. coordinates
    with np.load(Ffile) as fnpz:
        x_out = fnpz["x_out"];  y_out = fnpz["y_out"];  z_out = fnpz["z_out"]

    rho_out = np.sqrt(x_out**2 + y_out**2)
    z_vals  = z_out

    # flatten for Triangulation
    rho_f  = rho_out.ravel()
    z_f    = z_vals.ravel()
    Tpar_f = T_par.ravel()
    Tperp_f = T_perp.ravel()
    
    A_f_ = Tperp_f[Tpar_f>0]/Tpar_f[Tpar_f>0]
    A_f = A_f_[A_f_ > 0]

    print('Tpar max and min, 2.5, 97.5 percentile =',np.nanmax(Tpar_f[Tpar_f>0]),np.nanmin(Tpar_f[Tpar_f>0]), np.nanpercentile(Tpar_f[Tpar_f>0], 2.5), np.nanpercentile(Tpar_f[Tpar_f>0], 95) )
    print('Tperp max and min,2.5, 97.5percentile =',np.nanmax(Tperp_f[Tperp_f>0]),np.nanmin(Tperp_f[Tperp_f>0]), np.nanpercentile(Tperp_f[Tperp_f>0], 2.5), np.nanpercentile(Tperp_f[Tperp_f>0], 95))
    
    print('A max and min, 2.5, 97.5 percentile =',np.nanmax(A_f),np.nanmin(A_f), np.nanpercentile(A_f, 2.5), np.nanpercentile(A_f, 97.5))
    
    # log‑spaced levels (same range for ‖ and ⊥ just for consistency)
    levels = log_levels(np.concatenate([Tpar_f[Tpar_f>0], Tperp_f[Tperp_f>0]]),n=100)

    # 3. plot
    plot_field(rho_f, z_f, Tpar_f, "T‖", model, levels)
    plot_field(rho_f, z_f, Tperp_f, "T⊥", model, levels)
