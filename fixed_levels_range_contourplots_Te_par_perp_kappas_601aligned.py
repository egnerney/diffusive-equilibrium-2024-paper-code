# -*- coding: utf-8 -*-
"""
Created on Fri Apr 18 01:24:36 2025

@author: Owner
"""

#!/usr/bin/env python3
# -------------------------------------------------------------------------
# ρ‑z contour plots of electron T‖ and T⊥ (0.1–100 eV, 100 log levels)
# -------------------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from matplotlib.colors import LogNorm
from pathlib import Path

# -------------------------------------------------------------------------
# data‑file names ----------------------------------------------------------
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
}

# -------------------------------------------------------------------------
# fixed log‑spaced contour levels (0.1–100 eV, 100 steps) -----------------
LEVELS = np.logspace(np.log10(0.9), np.log10(100.0), 1000)
NORM   = LogNorm(vmin=LEVELS[0], vmax=LEVELS[-1])   # shared color scale

def plot_component(rho, z, T, model, comp):
    """One tricontourf plot (ρ, z) → T."""
    mask = T > 0                       # ignore zeros / negatives
    tri  = Triangulation(rho[mask], z[mask])

    fig, ax = plt.subplots(figsize=(7, 6))
    cf = ax.tricontourf(
        tri, T[mask],
        levels=LEVELS,
        norm=NORM,
        extend="both",
    )
    cbar = fig.colorbar(cf)
    cbar.set_label("Temperature (eV)   [log scale]")
    ax.set_xlabel(r"$\rho\;(\mathrm{R_J})$")
    ax.set_ylabel(r"$z\;(\mathrm{R_J})$")
    ax.set_title(f"{model}: electron {comp}")
    ax.set_aspect("equal", adjustable="box")
    plt.tight_layout()

    png = f"{model.replace(' ', '_')}_electron_{comp}.png"
    fig.savefig(png, dpi=300)
    print(f"saved → {png}")
    plt.show()

# -------------------------------------------------------------------------
# main loop ---------------------------------------------------------------
# -------------------------------------------------------------------------
for model, files in MODELS.items():
    Tfile, Ffile = Path(files["Tfile"]), Path(files["Ffile"])
    if not (Tfile.exists() and Ffile.exists()):
        print(f"** {model}: missing file(s) – skipping **")
        continue

    # 1) temperatures
    with np.load(Tfile) as Tnpz:
        T_par  = Tnpz["elec_par"]
        T_perp = Tnpz["elec_perp"]

    # 2) geometry
    with np.load(Ffile) as Fnpz:
        x_out, y_out, z_out = Fnpz["x_out"], Fnpz["y_out"], Fnpz["z_out"]

    rho = np.sqrt(x_out**2 + y_out**2).ravel()
    z   = z_out.ravel()

    # 3) plot T‖ and T⊥
    plot_component(rho, z, T_par.ravel(),  model, "T‖")
    plot_component(rho, z, T_perp.ravel(), model, "T⊥")
