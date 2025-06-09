import numpy as np
from pathlib import Path

# ---------------------------------------------------------------------------
# Helper to find extrema and their locations (non‑zero minima)
# ---------------------------------------------------------------------------
def extrema_with_location(T_par, T_perp, lat_deg, thresh=0.0):
    """
    Returns:
        (max_val, max_fl, max_lat),
        (min_val, min_fl, min_lat)  – min ignores values ≤ thresh
    """
    # ------- flatten helpers ---------
    flat_par  = T_par.ravel()
    flat_perp = T_perp.ravel()

    # ‣ Combined maximum (par or perp)
    if flat_par.max() >= flat_perp.max():
        max_val  = flat_par.max()
        max_idx  = np.argmax(flat_par)
        arr_shape = T_par.shape
    else:
        max_val  = flat_perp.max()
        max_idx  = np.argmax(flat_perp)
        arr_shape = T_perp.shape
    max_fl, max_lat_i = np.unravel_index(max_idx, arr_shape)

    # ‣ Combined minimum > thresh
    mask_par  = flat_par  > thresh
    mask_perp = flat_perp > thresh
    # concatenate masked values + original indices
    cand_vals  = np.concatenate([ flat_par[mask_par],  flat_perp[mask_perp] ])
    cand_idxs  = np.concatenate([ np.where(mask_par)[0], np.where(mask_perp)[0] ])

    if cand_vals.size == 0:                              # no positive temps
        min_val, min_fl, min_lat_i = np.nan, None, None
    else:
        j = np.argmin(cand_vals)
        min_val = cand_vals[j]
        min_idx = cand_idxs[j]
        # determine if idx belongs to par or perp for correct shape
        if min_idx < flat_par.size:
            arr_shape = T_par.shape
        else:
            min_idx -= flat_par.size
            arr_shape = T_perp.shape
        min_fl, min_lat_i = np.unravel_index(min_idx, arr_shape)

    # Latitude in degrees
    max_lat = lat_deg[max_fl, max_lat_i]
    min_lat = lat_deg[min_fl, min_lat_i] if min_fl is not None else np.nan

    return (max_val, max_fl, max_lat), (min_val, min_fl, min_lat)

# ---------------------------------------------------------------------------
# File sets for the three distribution models
# ---------------------------------------------------------------------------
MODELS = {
    "Fried‑Egg": dict(
        Tfile   = "final_new_nominal_model_T_out_nominal_model_4-10_isoT_isokappa_A=1_fried_egg_all_nominal_model.npz",
        Ffile   = "final_new_nominal_model_field_data_nominal_model_4-10_isoT_isokappa_A=1_fried_egg_all_nominal_model.npz",
    ),
    "Standard‑Kappa": dict(
        Tfile   = "final_new_nominal_model_T_out_nominal_model_4-10_isoT_isokappa_A=1_standard_kappa_all_nominal_model.npz",
        Ffile   = "final_new_nominal_model_field_data_nominal_model_4-10_isoT_isokappa_A=1_standard_kappa_all_nominal_model.npz",
    ),
    "Product‑Kappa": dict(
        Tfile   = "final_new_nominal_model_T_out_nominal_model_4-10_isoT_isokappa_A=1_product_kappa_all_nominal_model.npz",
        Ffile   = "final_new_nominal_model_field_data_nominal_model_4-10_isoT_isokappa_A=1_product_kappa_all_nominal_model.npz",
    ),
}

THRESH = 0.0            # ignore temps ≤ this value for the minima

print("Electron‑temperature extrema (eV) – 601 field lines, all latitudes")
print("================================================================")
print("{:<15s} {:>9s} {:>6s} {:>9s}   {:>9s} {:>6s} {:>9s}"
      .format("Model", "T_max", "FL#", "Lat°", "T_min>0", "FL#", "Lat°"))

for model, files in MODELS.items():
    Tfile, Ffile = Path(files["Tfile"]), Path(files["Ffile"])

    if not (Tfile.is_file() and Ffile.is_file()):
        print(f"{model:<15s}  **missing file(s)**")
        continue

    # --- 1. load temperatures ------------------------------------------------
    with np.load(Tfile) as tnpz:
        T_par  = tnpz["elec_par"]     # shape (601, n_lat)
        T_perp = tnpz["elec_perp"]

    # --- 2. load field‑line coordinates -------------------------------------
    with np.load(Ffile) as fnpz:
        x_out = fnpz["x_out"]
        y_out = fnpz["y_out"]
        z_out = fnpz["z_out"]

    # r, ρ, latitude
    r_out      = np.sqrt(x_out**2 + y_out**2 + z_out**2)
    rho_out    = np.sqrt(x_out**2 + y_out**2)
    lat_out_deg = np.degrees(np.arcsin(np.where(r_out > 0, z_out / r_out, 0.0)))

    # --- 3. extrema + location ----------------------------------------------
    (Tmax, fl_max, lat_max), (Tmin, fl_min, lat_min) = \
        extrema_with_location(T_par, T_perp, lat_out_deg, thresh=THRESH)

    # --- 4. print summary ----------------------------------------------------
    print("{:<15s} {:9.3g} {:6d} {:9.2f}   {:9.3g} {:6d} {:9.2f}"
          .format(model, Tmax, fl_max, lat_max, Tmin, fl_min, lat_min))
