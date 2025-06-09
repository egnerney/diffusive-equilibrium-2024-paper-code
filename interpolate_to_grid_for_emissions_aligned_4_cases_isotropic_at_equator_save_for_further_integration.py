import numpy as np
import os
from scipy.interpolate import griddata, RegularGridInterpolator



def main():
    """
    Updated final script:
      1) Uses rho1D_unique = [4..6 by 0.01, 6.1..10 by 0.1] => 241 total
      2) Uses phi1D_unique = linspace(270, 360, 11) => includes phi=360
      3) Skips 'oph' & 'eh' for standard_kappa, product_kappa, fried_egg
         but includes them for maxwellian if present
      4) z1D_unique remains piecewise (53 pts from 0..3)
      5) s in [0..10 by 0.01] => 1001 steps
      6) Avoids shape mismatch by building X_3D, Y_3D, Z_3D carefully
      7) Filters data by distribution constraints, interpolates, replicates
         in phi, line-of-sight samples => saves .npz
    """
    
    

    # Four distribution types we process:
    dist_list = ['standard_kappa','product_kappa','fried_egg','maxwellian']

    # Map each distribution to the relevant file names:
    file_map = {
        'standard_kappa': {
            'n_file':    'final_new_nominal_model_n_out_nominal_model_4-10_isoT_isokappa_A=1_standard_kappa_all_nominal_model.npz',
            'T_file':    'final_new_nominal_model_T_out_nominal_model_4-10_isoT_isokappa_A=1_standard_kappa_all_nominal_model.npz',
            'k_file':    'final_new_nominal_model_kappa_out_nominal_model_4-10_isoT_isokappa_A=1_standard_kappa_all_nominal_model.npz',
            'field_file':'final_new_nominal_model_field_data_nominal_model_4-10_isoT_isokappa_A=1_standard_kappa_all_nominal_model.npz'
        },
        'product_kappa': {
            'n_file':    'final_new_nominal_model_n_out_nominal_model_4-10_isoT_isokappa_A=1_product_kappa_all_nominal_model.npz',
            'T_file':    'final_new_nominal_model_T_out_nominal_model_4-10_isoT_isokappa_A=1_product_kappa_all_nominal_model.npz',
            'k_file':    'final_new_nominal_model_kappa_out_nominal_model_4-10_isoT_isokappa_A=1_product_kappa_all_nominal_model.npz',
            'field_file':'final_new_nominal_model_field_data_nominal_model_4-10_isoT_isokappa_A=1_product_kappa_all_nominal_model.npz'
        },
        'fried_egg': {
            'n_file':    'final_new_nominal_model_n_out_nominal_model_4-10_isoT_isokappa_A=1_fried_egg_all_nominal_model.npz',
            'T_file':    'final_new_nominal_model_T_out_nominal_model_4-10_isoT_isokappa_A=1_fried_egg_all_nominal_model.npz',
            'k_file':    'final_new_nominal_model_kappa_out_nominal_model_4-10_isoT_isokappa_A=1_fried_egg_all_nominal_model.npz',
            'field_file':'final_new_nominal_model_field_data_nominal_model_4-10_isoT_isokappa_A=1_fried_egg_all_nominal_model.npz'
        },
        'maxwellian': {
            'n_file':    'final_new_nominal_model_n_out_nominal_model_4-10_isoT_A=1_max_all_nominal_model.npz',
            'T_file':    'final_new_nominal_model_T_out_nominal_model_4-10_isoT_A=1_max_all_nominal_model.npz',
            'k_file':    'final_new_nominal_model_kappa_out_nominal_model_4-10_isoT_A=1_max_all_nominal_model.npz',
            'field_file':'final_new_nominal_model_field_data_nominal_model_4-10_isoT_A=1_max_all_nominal_model.npz'
        }
    }

    
    #species_names = ['op','o2p','sp','s2p','s3p','hp','nap','oph','eh','elec']
    #species_names = ['op','sp','s2p','oph','eh','elec']
    
    

    # 1) Build final (rho,z) grids
    #    241 points in rho:
    """
    rho1D_unique = np.concatenate([
        np.linspace(4.0, 6.0, 201),
        np.linspace(6.1, 10.0, 40)
    ])  # total 241
    # 53 points in z:
    z1D_unique = np.concatenate([
        np.linspace(0.0, 0.29, 30),
        np.linspace(0.3, 2.0, 18),
        np.linspace(2.2, 3.0, 5)
    ])
    """
    rho1D_unique =np.linspace(4.0,10.0,601)
    z1D_unique =np.linspace(0.0,2.5,251)
    
    # 22 points in phi => [270..360], inclusive
    phi1D_unique = np.linspace(270., 360., 22)

    Nx = len(rho1D_unique)  # 241
    Nz = len(z1D_unique)    # 53
    Nphi = len(phi1D_unique)  # 11

    # 2) Lines of sight s in [0..10.01] (0.01 step => ~1002 steps)
    ds = 0.01
    s_1D = np.arange(0., 10.0 + ds*0.5, ds)
    Ns = len(s_1D)

    # We create a 2D mesh for (x0,z0) => shape(241,53)
    X2D, Z2D = np.meshgrid(rho1D_unique, z1D_unique, indexing='ij')  # (241,53)

    # Expand to 3D by replicating along the s dimension => shape(241,53,Ns)
    X_3D = np.repeat(X2D[..., None], Ns, axis=-1)
    Z_3D = np.repeat(Z2D[..., None], Ns, axis=-1)

    # Y depends only on s => shape(241,53,Ns)
    Y_3D = -10.0 + np.tile(s_1D, (Nx*Nz, 1)).reshape(Nx, Nz, Ns)

    def interpolate_rhoz(rhof, zf, vals, rho_target, z_target):
        """
        2D interpolation from scattered (rho,z)->vals onto (rho_target,z_target).
        Fill out-of-range with NaN if no points remain after filtering.
        """
        rf = rhof.flatten()
        zff= zf.flatten()
        vf = vals.flatten()

        mask = np.isfinite(vf)
        rf   = rf[mask]
        zff  = zff[mask]
        vf   = vf[mask]

        if rf.size == 0:
            return np.full((len(rho_target), len(z_target)), np.nan)

        pts = np.column_stack((rf, zff))
        Rg, Zg = np.meshgrid(rho_target, z_target, indexing='ij')
        out = griddata(pts, vf, (Rg, Zg), method='linear', fill_value=np.nan)
        return out

    def replicate_phi(arr_2d, nphi):
        """
        Replicate (Nrho, Nz)->(Nrho, nphi, Nz) along phi dimension.
        """
        return np.repeat(arr_2d[:, None, :], repeats=nphi, axis=1)

    # 3) Main loop
    #dist_list = ['standard_kappa','product_kappa','fried_egg','maxwellian']
    #for dist_type in dist_list:
    #dist_type = 'standard_kappa'
    #dist_type = 'product_kappa'
    #dist_type = 'fried_egg'
    dist_type = 'maxwellian'
    species_names = ['op','sp','s2p','oph','eh','elec']
    #species_names = ['op','sp','s2p','elec']
    """
    print(f"\n=== Distribution: {dist_type} ===")
    if dist_type not in file_map:
        print(f"Unknown distribution {dist_type}, skipping.")
        continue
    """

    files = file_map[dist_type]
    n_file     = files['n_file']
    T_file     = files['T_file']
    k_file     = files['k_file']
    field_file = files['field_file']

    # Check existence
    """
    if not (os.path.exists(n_file) and os.path.exists(T_file)
            and os.path.exists(k_file) and os.path.exists(field_file)):
        print("One or more input files is missing for", dist_type)
        continue
    """

    print("Loading:", n_file, T_file, k_file, field_file)
    n_data = np.load(n_file)
    T_data = np.load(T_file)
    #k_data = np.load(k_file)
    f_data = np.load(field_file)

    # Field data => shape(601,1401)
    x_out = f_data['x_out']
    y_out = f_data['y_out']
    z_out = f_data['z_out']

    # Flatten => scattered (rho,z)
    Xf = x_out.ravel()
    Yf = y_out.ravel()
    Zf = z_out.ravel()

    rhof = np.sqrt(Xf**2. + Yf**2.)
    zff  = Zf

    # Prepare final dictionaries
    n_3D_dict      = {}
    Tpar_3D_dict   = {}
    Tperp_3D_dict  = {}
    kpar_3D_dict   = {}
    kperp_3D_dict  = {}

    #skip_hot_species = (dist_type in ['standard_kappa','product_kappa','fried_egg'])

    for sp in species_names:
        print(sp)
        # skip oph & eh for the first 3 distributions
        """
        if skip_hot_species and sp in ['oph','eh']:
            print(f"Skipping {sp} for {dist_type} (no hot species).")
            continue

        if sp not in n_data:
            print(f"{sp} not in n_data => skipping.")
            continue
        """

        n_in   = n_data[sp]
        #Tpar_in= T_data.get(f"{sp}_par",  np.zeros_like(n_in))
        Tperp_in= T_data.get(f"{sp}_perp", np.zeros_like(n_in))
        # kpar_in= k_data.get(f"{sp}_par",  np.zeros_like(n_in))
        #kperp_in= k_data.get(f"{sp}_perp", np.zeros_like(n_in))

        # Filter
        """
        n_filt, Tpar_filt, Tperp_filt, kpar_filt, kperp_filt = \
            filter_data_for_distribution(
                dist_type,
                n_in, Tpar_in, Tperp_in,
                kpar_in, kperp_in
            )

        """
        # n_filt, Tpar_filt, Tperp_filt, kpar_filt, kperp_filt = n_in, Tpar_in, Tperp_in, kpar_in, kperp_in
        #n_filt, Tperp_filt, kperp_filt = n_in,  Tperp_in, kperp_in
        #n_filt, Tpar_filt, Tperp_filt, kpar_filt, kperp_filt = n_in, Tpar_in, Tperp_in, kpar_in, kperp_in
        #n_filt, Tpar_filt, Tperp_filt, kperp_filt = n_in, Tpar_in, Tperp_in, kperp_in
        n_filt, Tperp_filt = n_in,  Tperp_in
        
        

        # Interpolate => shape(241,53)
        n_rz   = interpolate_rhoz(rhof, zff, n_filt,   rho1D_unique, z1D_unique)
        #Tpar_rz= interpolate_rhoz(rhof, zff, Tpar_filt,rho1D_unique, z1D_unique)
        Tperp_rz=interpolate_rhoz(rhof, zff, Tperp_filt,rho1D_unique, z1D_unique)
        #kpar_rz= interpolate_rhoz(rhof, zff, kpar_filt,rho1D_unique, z1D_unique)
        #kperp_rz=interpolate_rhoz(rhof, zff, kperp_filt,rho1D_unique, z1D_unique)

        # Replicate in phi => shape(241, nphi, 53) => (241,11,53)
        n_3D_sp      = replicate_phi(n_rz,     Nphi)
        #Tpar_3D_sp   = replicate_phi(Tpar_rz,  Nphi)
        Tperp_3D_sp  = replicate_phi(Tperp_rz, Nphi)
        #kpar_3D_sp   = replicate_phi(kpar_rz,  Nphi)
        #kperp_3D_sp  = replicate_phi(kperp_rz, Nphi)

        n_3D_dict[sp]      = n_3D_sp
        #Tpar_3D_dict[sp]   = Tpar_3D_sp
        Tperp_3D_dict[sp]  = Tperp_3D_sp
        #kpar_3D_dict[sp]   = kpar_3D_sp
        #kperp_3D_dict[sp]  = kperp_3D_sp

    # Build 3D interpolators in (rho,phi,z)
    interpolators_n     = {}
    interpolators_Tpar  = {}
    interpolators_Tperp = {}
    interpolators_kpar  = {}
    interpolators_kperp = {}

    for sp in n_3D_dict.keys():
        print(sp)
        # each array => shape(241,11,53)
        interpolators_n[sp] = RegularGridInterpolator(
            (rho1D_unique, phi1D_unique, z1D_unique),
            n_3D_dict[sp], method='linear', bounds_error=False, fill_value=np.nan
        )
        
        """
        interpolators_Tpar[sp] = RegularGridInterpolator(
            (rho1D_unique, phi1D_unique, z1D_unique),
            Tpar_3D_dict[sp], method='linear', bounds_error=False, fill_value=np.nan
        )
        """
        
        
        interpolators_Tperp[sp] = RegularGridInterpolator(
            (rho1D_unique, phi1D_unique, z1D_unique),
            Tperp_3D_dict[sp], method='linear', bounds_error=False, fill_value=np.nan
        )
        """
        
        interpolators_kpar[sp] = RegularGridInterpolator(
            (rho1D_unique, phi1D_unique, z1D_unique),
            kpar_3D_dict[sp], method='linear', bounds_error=False, fill_value=np.nan
        )
        
        interpolators_kperp[sp] = RegularGridInterpolator(
            (rho1D_unique, phi1D_unique, z1D_unique),
            kperp_3D_dict[sp], method='linear', bounds_error=False, fill_value=np.nan
        )
        """

    # Now sample lines-of-sight: shape(241,53,Ns)
    Xf = X_3D.flatten()  # => 241*53*Ns
    Yf = Y_3D.flatten()
    Zf = Z_3D.flatten()

    # Convert to cylindrical => (rho, phi, z)
    rho_f = np.sqrt(Xf**2 + Yf**2)
    phi_f = np.degrees(np.arctan2(Yf, Xf))
    phi_f[phi_f < 0] += 360.0

    pts_3D = np.column_stack((rho_f, phi_f, Zf))  # shape => (241*53*Ns, 3)

    # Interpolate each species
    n_final   = {}
    Tpar_final= {}
    Tperp_final= {}
    kpar_final= {}
    kperp_final= {}

    for sp in interpolators_n.keys():
        print(sp)
        i_n    = interpolators_n[sp](pts_3D)
        # i_Tpar = interpolators_Tpar[sp](pts_3D)
        i_Tperp= interpolators_Tperp[sp](pts_3D)
        #i_kpar = interpolators_kpar[sp](pts_3D)
        #i_kperp= interpolators_kperp[sp](pts_3D)

        # Reshape => (241,53,Ns)
        n_final[sp]      = i_n.reshape((Nx,Nz,Ns))
        #Tpar_final[sp]   = i_Tpar.reshape((Nx,Nz,Ns))
        Tperp_final[sp]  = i_Tperp.reshape((Nx,Nz,Ns))
        #kpar_final[sp]   = i_kpar.reshape((Nx,Nz,Ns))
        #kperp_final[sp]  = i_kperp.reshape((Nx,Nz,Ns))

    # Save final results
    outfilename = f"los_results_{dist_type}_with_constraints.npz"
    print("Saving =>", outfilename)
    save_dict = {
        'X_3D': X_3D,
        'Y_3D': Y_3D,
        'Z_3D': Z_3D,
        's_1D': s_1D,
        'rho1D': rho1D_unique,
        'z1D': z1D_unique,
        'phi1D': phi1D_unique
    }
    for sp in n_final.keys():
        print(sp)
        save_dict[f"n_{sp}"]      = n_final[sp]
        #save_dict[f"Tpar_{sp}"]   = Tpar_final[sp]
        save_dict[f"Tperp_{sp}"]  = Tperp_final[sp]
        #save_dict[f"kpar_{sp}"]   = kpar_final[sp]
        #save_dict[f"kperp_{sp}"]  = kperp_final[sp]

    np.savez_compressed(outfilename, **save_dict)

    print("\nAll done!")

if __name__ == "__main__":
    main()
