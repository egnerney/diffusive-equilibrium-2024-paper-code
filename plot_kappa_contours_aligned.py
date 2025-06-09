# -*- coding: utf-8 -*-
"""
Created on Fri Apr 18 01:20:53 2025

@author: Owner
"""

#!/usr/bin/env python3
# -------------------------------------------------------------------------
# Contour maps of electron kappa‖ and kappa⊥ vs ρ–z   (Jovian magnetosphere model)
# -------------------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from pathlib import Path

# -------------------------------------------------------------------------
# file names for the three distribution models ----------------------------
# -------------------------------------------------------------------------

# n output files
# final_new_nominal_model_n_out_nominal_model_4-10_isoT_isokappa_A=1_standard_kappa_all_nominal_model.npz
# final_new_nominal_model_n_out_nominal_model_4-10_isoT_isokappa_A=1_product_kappa_all_nominal_model.npz
# final_new_nominal_model_n_out_nominal_model_4-10_isoT_isokappa_A=1_fried_egg_all_nominal_model.npz

# T output files
# final_new_nominal_model_T_out_nominal_model_4-10_isoT_isokappa_A=1_standard_kappa_all_nominal_model
# final_new_nominal_model_T_out_nominal_model_4-10_isoT_isokappa_A=1_product_kappa_all_nominal_model
# final_new_nominal_model_T_out_nominal_model_4-10_isoT_isokappa_A=1_fried_egg_all_nominal_model.npz


# same field data files
# final_new_nominal_model_field_data_nominal_model_4-10_isoT_isokappa_A=1_standard_kappa_all_nominal_model.npz
# final_new_nominal_model_field_data_nominal_model_4-10_isoT_isokappa_A=1_product_kappa_all_nominal_model.npz
# final_new_nominal_model_field_data_nominal_model_4-10_isoT_isokappa_A=1_fried_egg_all_nominal_model.npz

# kappa output files
# final_new_nominal_model_kappa_out_nominal_model_4-10_isoT_isokappa_A=1_standard_kappa_all_nominal_model.npz
# final_new_nominal_model_kappa_out_nominal_model_4-10_isoT_isokappa_A=1_product_kappa_all_nominal_model.npz
# final_new_nominal_model_kappa_out_nominal_model_4-10_isoT_isokappa_A=1_fried_egg_all_nominal_model.npz
MODELS = {
    "Fried‑Egg": dict(
        kappafile="final_new_nominal_model_kappa_out_nominal_model_4-10_isoT_isokappa_A=1_fried_egg_all_nominal_model.npz",
        Ffile="final_new_nominal_model_field_data_nominal_model_4-10_isoT_isokappa_A=1_fried_egg_all_nominal_model.npz",
    ),
    "Standard‑Kappa": dict(
        kappafile="final_new_nominal_model_kappa_out_nominal_model_4-10_isoT_isokappa_A=1_standard_kappa_all_nominal_model.npz",
        Ffile="final_new_nominal_model_field_data_nominal_model_4-10_isoT_isokappa_A=1_standard_kappa_all_nominal_model.npz",
    ),
    "Product‑Kappa": dict(
        kappafile="final_new_nominal_model_kappa_out_nominal_model_4-10_isoT_isokappa_A=1_product_kappa_all_nominal_model.npz",
        Ffile="final_new_nominal_model_field_data_nominal_model_4-10_isoT_isokappa_A=1_product_kappa_all_nominal_model.npz",
    ),
}

# -------------------------------------------------------------------------
def log_levels(arr, n=100):
    """Return *n* log‑spaced contour levels covering [arr_min, arr_max]."""
    positives = arr[arr > 0]
    vmin, vmax = positives.min(), positives.max()
    v95 = np.nanpercentile(arr, 95)
    return np.logspace(np.log10(vmin), np.log10(v95), n)

def plot_field(rho, z, kappa, comp, model, levels):
    """Triangular contour plot of one temperature component."""
    # Keep only >0 temps for color mapping
    mask = kappa > 0
    tri = Triangulation(rho[mask], z[mask])

    fig, ax = plt.subplots(figsize=(7, 6))
    cf = ax.tricontourf(tri, kappa[mask], levels=levels, extend="both")
    cbar = fig.colorbar(cf)
    cbar.set_label(r"$\kappa$")
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
    kappafile = Path(files["kappafile"])
    Ffile = Path(files["Ffile"])

    # 1. temperatures
    with np.load(kappafile) as knpz:
        kappa_par  = knpz["elec_par"]   # shape (601, 1401) # n_lat=1401 (-70 to 70 degreees latitude at 0.1 degree dlat)
        kappa_perp = knpz["elec_perp"]

    print(kappa_par.shape,kappa_perp.shape) # returns (601,1401) which is correct...
    # 2. coordinates
    with np.load(Ffile) as fnpz:
        x_out = fnpz["x_out"];  y_out = fnpz["y_out"];  z_out = fnpz["z_out"]

    rho_out = np.sqrt(x_out**2 + y_out**2)
    z_vals  = z_out

    # flatten for Triangulation
    rho_f  = rho_out.ravel()
    z_f    = z_vals.ravel()

    kappa_perp_f = kappa_perp.ravel()
    
    
    print('kappa perp max and min, 2.5, 97.5 percentile =',np.nanmax(kappa_perp_f[kappa_perp_f>0]),np.nanmin(kappa_perp_f[kappa_perp_f>0]), np.nanpercentile(kappa_perp_f[kappa_perp_f>0], 2.5) , np.nanpercentile(kappa_perp_f[kappa_perp_f>0], 97.5))
   
    
    if model != 'Fried‑Egg':
        kappa_par_f = kappa_par.ravel()
        kratio = kappa_perp_f[kappa_par_f >0]/ kappa_par_f[kappa_par_f >0]
        print('kappa par max and min, 2.5, 97.5 percentile =',np.nanmax(kappa_par_f[kappa_par_f>0]),np.nanmin(kappa_par_f[kappa_par_f>0]), np.nanpercentile(kappa_par_f[kappa_par_f>0], 2.5), np.nanpercentile(kappa_par_f[kappa_par_f>0], 97.5) )
        print('kappa ratio max and min, 2.5, 97.5 percentile =',np.nanmax(kratio[kratio > 0]),np.nanmin(kratio[kratio >0]), np.nanpercentile(kratio, 2.5), np.nanpercentile(kratio, 97.5))
        levels = log_levels(np.concatenate([kappa_par_f[kappa_par_f>0], kappa_perp_f[kappa_perp_f>0]]),n=100)
        plot_field(rho_f, z_f, kappa_par_f, "kappa_par", model, levels)
    else:
        levels = log_levels( kappa_perp_f[kappa_perp_f>0],n=100)

    # log‑spaced levels (same range for ‖ and ⊥ just for consistency)
    
    # 3. plot
    plot_field(rho_f, z_f, kappa_perp_f, "kappa_perp", model, levels)

    
