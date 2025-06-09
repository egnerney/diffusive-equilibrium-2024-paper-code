# ───────────────────────────  export_all.py  ───────────────────────────
#!/usr/bin/env python3
"""
Convert three compressed NumPy archives to raw float‑64 Fortran‑order binaries
compatible with IDL / Fortran (column‑major).

Files produced
--------------
  standard‑κ set .......... *_sk.bin   (5 files)
  fried‑egg    set ........ *_fe.bin   (6 files)
  product‑κ    set ........ *_pk.bin   (7 files)

Each filename ends with the same tag used below (sk, fe, pk).
"""

import numpy as np
from pathlib import Path

# ------------------------------------------------------------------ specs
DATASETS = [
    # 1) standard‑kappa (T × κ)
    dict(
        tag="sk",
        npz="iso_std_kappa_to_5Mxw_phys_cutoff_Txkappa_14x24.npz",
        mapping=dict(
            T="T_grid",             # → T_grid_f64_sk.bin
            kappa="kappa_grid",
            T_i="T_i",
            f_i="f_i",
            max_err="mErr",
        ),
    ),
    # 2) fried‑egg (A × T⊥ × κ)
    dict(
        tag="fe",
        npz="fried_egg_5Mxw_fit_phys_cutoff_7x13x23_AxTperpxkappa.npz",
        mapping=dict(
            A="A_grid",
            Tperp="Tperp_grid",
            kappa="kappa_grid",
            T_i="T_i",
            f_i="f_i",
            max_err="mErr",
        ),
    ),
    # 3) product‑kappa (A × T⊥ × λ × κ⊥)
    dict(
        tag="pk",
        npz="product_kappa_to_5Mxw_phys_cutoff_A10_T14_L9_K17.npz",
        mapping=dict(
            A="A_grid",
            Tperp="Tperp_grid",
            lambda_="lambda_grid",
            kappa_perp="kperp_grid",
            T_i="T_i",
            f_i="f_i",
            max_err="mErr",
        ),
    ),
]

# ---------------------------------------------------------------- helper
def to_fortran_bin(arr: np.ndarray, fname: str) -> None:
    """Write *arr* as float‑64, column‑major, no header."""
    np.asfortranarray(arr).ravel(order="F").astype("float64").tofile(fname)


# ---------------------------------------------------------------- export
for ds in DATASETS:
    archive = Path(ds["npz"]).expanduser()
    if not archive.exists():
        raise FileNotFoundError(archive)
    data = np.load(archive)

    print(f"→  {archive}  ({ds['tag']})")
    for npz_key, base_name in ds["mapping"].items():
        arr = data[npz_key]
        out_file = f"{base_name}_f64_{ds['tag']}.bin"
        to_fortran_bin(arr, out_file)
        print(f"   wrote {out_file}")
print("All done.")
