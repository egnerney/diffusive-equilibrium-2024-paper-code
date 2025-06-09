# -*- coding: utf-8 -*-
"""
Created on Wed Apr 30 02:19:22 2025

@author: Owner
"""

def filter_data_for_distribution(
    dist_type,
    n_in, Tpar_in, Tperp_in,
    kpar_in, kperp_in
):
    """
    Apply distribution-specific constraints and unify the parallel/perp
    values as needed. Return masked arrays (with NaNs) for:
        n, Tpar, Tperp, kpar, kperp

    dist_type can be:
      - 'standard_kappa'   (isotropic => T = T_perp, kappa = k_perp)
      - 'product_kappa'    (fully anisotropic => ratio A, lambda)
      - 'fried_egg'        (Maxwellian parallel => no kpar, only kperp)
      - 'maxwellian'       (pure isotropic Maxwellian => no kappa)

    We always require:
       n > 0 & finite,
       Tpar/Tperp/kpar/kperp finite.
    Then distribution-specific constraints on T, kappa, anisotropy, etc.
    Anything failing these => set to NaN.
    """

    # Flatten everything
    n   = n_in.flatten()
    Tp  = Tpar_in.flatten()
    Tpp = Tperp_in.flatten()
    kp  = kpar_in.flatten()
    kpp = kperp_in.flatten()

    # Basic valid data: finite, n>0
    mask = np.isfinite(n) #& (n > 0)
    # Also require Tpar, Tperp, kpar, kperp to be finite
    mask &= np.isfinite(Tp) & (Tp > 0) &  np.isfinite(Tpp) & (Tpp > 0)
    # mask &= np.isfinite(kp) & np.isfinite(kpp)

    # Prepare final arrays (flattened) with NaN by default
    n_out   = np.full_like(n, np.nan)
    Tpar    = np.full_like(Tp, np.nan)
    Tperp   = np.full_like(Tpp, np.nan)
    kparf   = np.full_like(kp, np.nan)
    kperpf  = np.full_like(kpp, np.nan)

    # Copy n where valid
    #n_out[mask] = n[mask]

    if dist_type == 'standard_kappa':
        # Isotropic => T = T_perp, kappa = k_perp
        # 0.971 < T < 68.19, 1.941 < kappa < 299.9
        # Copy n where valid
        #n_temp = n[mask]
        #Ttemp = Tpp[mask]
        #ktemp = kpp[mask]
        mask &= (Tpp > 0.971) & (Tpp < 68.19)
        mask &= (kpp > 1.941) & np.isfinite(kpp) & (kpp < 299.9)
        
        Tpar[mask]   = Tp[mask]
        Tperp[mask]  = Tpp[mask]
        kparf[mask]  = kp[mask]
        kperpf[mask] = kpp[mask]
        n_out[mask] = n[mask]

    elif dist_type == 'product_kappa':

        with np.errstate(divide='ignore', invalid='ignore'):
            A   = Tpp / Tp
            lam = kpp / kp

        mask &= (Tpp > 0.971) & (Tpp < 143.3)
        mask &= (A > 0.281) & (A < 1.089) & np.isfinite(A)
        mask &= (lam >= 1.0) & (lam < 2.59) & np.isfinite(lam) & np.isfinite(lam)
        mask &= (kpp > 1.701) & (kpp <= 300.0) & np.isfinite(kpp)



        Tpar[mask]   = Tpp[mask]/A[mask]
        Tperp[mask]  = Tpp[mask]
        kparf[mask]  = kpp[mask]/lam[mask]
        kperpf[mask] = kpp[mask]
        n_out[mask] = n[mask]


    elif dist_type == 'fried_egg':


        with np.errstate(divide='ignore', invalid='ignore'):
            A = Tpp / Tp

        mask &= (A > 0.842) & (A < 1.119) & np.isfinite(A)
        mask &= (Tpp >= 1.0) & (Tpp < 12.469)
        mask &= (kpp > 1.231) & (kpp < 299.9) & np.isfinite(kpp)


        Tpar[mask]   = Tpp[mask]/A[mask]
        Tperp[mask]  = Tpp[mask]
        kparf[mask]  = kpp[mask]
        kperpf[mask] = kpp[mask]
        n_out[mask] = n[mask]

    elif dist_type == 'maxwellian':
        # Pure isotropic Maxwellian => Tpar = Tperp, no kappa

        #mask &= (Tpp > 1.0)  # just > 0

        #Tpar[mask]   = Tpp[mask]
        #Tperp[mask]  = Tpp[mask]
        # kparf[mask]  = kpp[mask]
        # kperpf[mask] = kpp[mask]
        # n_out[mask] = n[mask]
        
        Tpar  = Tpp
        Tperp  = Tpp
        kparf  = kpp
        kperpf = kpp
        n_out = n


    # Reshape to original shape
    sh = n_in.shape
    n_out   = n_out.reshape(sh)
    Tpar    = Tpar.reshape(sh)
    Tperp   = Tperp.reshape(sh)
    kparf   = kparf.reshape(sh)
    kperpf  = kperpf.reshape(sh)

    return n_out, Tpar, Tperp, kparf, kperpf
