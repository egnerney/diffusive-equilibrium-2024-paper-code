import numpy as np
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator, griddata
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.integrate import simpson
import matplotlib.colors as colors

#import matplotlib.image as mpimg
#from matplotlib import cm, ticker
import matplotlib.ticker as mticker
import scipy.ndimage as ndi

##############################################################################
#               MATPLOTLIB CONFIG: BOLD AXES, LABELS, TIGHT LAYOUT
##############################################################################
#mpl.rcParams["figure.dpi"] = 600
#mpl.rcParams["savefig.dpi"] = 600
mpl.rcParams["axes.labelweight"] = "bold"
mpl.rcParams["axes.titleweight"] = "bold"
mpl.rcParams["axes.linewidth"]   = 1.5
mpl.rcParams["xtick.labelsize"]  = 10
mpl.rcParams["ytick.labelsize"]  = 10
mpl.rcParams["xtick.major.width"]= 1.3
mpl.rcParams["ytick.major.width"]= 1.3
mpl.rcParams["font.weight"]      = "bold"
mpl.rcParams["axes.labelsize"]   = 11
mpl.rcParams["axes.titlesize"]   = 11



##############################################################################
#               1) LOAD EMISSION TABLES & INTERPOLATORS
##############################################################################

def load_standard_kappa_emiss():
    """Standard Kappa => (n_e, T, kappa) => shape (9,14,24,3)."""
    netotal_sk_grid = np.array([1,10,100,1000,2000,2500,3000,3500,4000], dtype=float)
    T_sk_grid       = np.array([0.97,1.,2.5,5.,7.5,10.,15.,20.,30.,40.,50.,60.,65.,69.], dtype=float)
    kappa_sk_grid   = np.array([1.51,1.6,1.7,1.8,1.9,2.,2.5,3.,4.,5.,6.,7.,10.,15.,20.,
                                25.,30.,40.,50.,75.,100.,150.,200.,300.], dtype=float)

    n_netotal, n_T_sk, n_kappa_sk = 9,14,24
    fname = 'emiss_table_sk_f64.bin'
    if not os.path.exists(fname):
        return None
    with open(fname,'rb') as f:
        arr = np.fromfile(f, dtype=np.float64)
    arr = arr.reshape((n_netotal, n_T_sk, n_kappa_sk, 3), order='F')

    # Build the scaled log interpolator
    #coord_grids = [netotal_sk_grid, T_sk_grid, kappa_sk_grid]
    return RegularGridInterpolator(
        (netotal_sk_grid, T_sk_grid, kappa_sk_grid),
        arr, method='linear', bounds_error=False, fill_value=np.nan
    )

def load_product_kappa_emiss():
    """Product Kappa => (n_e, A, T_perp, lambda, kappa_perp) => shape(9,10,14,9,17,3)."""
    netotal_grid     = np.array([1,10,100,1000,2000,2500,3000,3500,4000], dtype=float)
    A_pk_grid        = np.array([0.28,0.5,0.7,0.9,1.0,1.1,1.3,1.5,1.7,1.9], dtype=float)
    T_perp_pk_grid   = np.array([0.97,1.0,2.5,5.,7.5,10.,15.,20.,30.,50.,75.,100.,125.,144.], dtype=float)
    lambda_pk_grid   = np.array([1.,1.15,1.35,1.55,1.75,1.95,2.15,2.35,2.6], dtype=float)
    kappa_perp_pk_grid = np.array([1.01,1.1,1.3,1.7,2.,2.5,3.,4.,5.,7.,10.,15.,20.,30.,50.,150.,300.], dtype=float)

    n_netotal, nA, nT, nlam, nkp = 9,10,14,9,17
    fname = 'emiss_table_pk_f64.bin'
    if not os.path.exists(fname):
        return None
    with open(fname,'rb') as f:
        arr = np.fromfile(f, dtype=np.float64)
    arr = arr.reshape((n_netotal,nA,nT,nlam,nkp,3), order='F')

    # Build the scaled log interpolator
    #coord_grids = [netotal_grid, A_pk_grid, T_perp_pk_grid, lambda_pk_grid, kappa_perp_pk_grid]
    return RegularGridInterpolator(
        (netotal_grid, A_pk_grid, T_perp_pk_grid, lambda_pk_grid, kappa_perp_pk_grid),
        arr, method='linear', bounds_error=False, fill_value=np.nan
    )
def load_fried_egg_emiss():
    """Fried-Egg => (n_e, A, T_perp, kappa_perp) => shape(9,7,13,23,3)."""
    netotal_grid = np.array([1,10,100,1000,2000,2500,3000,3500,4000], dtype=float)
    A_fe_grid    = np.array([0.84,0.9,0.95,1.0,1.05,1.1,1.12], dtype=float)
    T_fe_grid    = np.array([0.99,1.0,2.,2.5,3.,4.,5.,6.,7.5,10.,11.,12.,12.5], dtype=float)
    kappa_fe_grid= np.array([1.23,1.3,1.5,1.7,1.8,1.9,2.,2.5,3.,4.,5.,7.,10.,15.,20.,
                             25.,30.,50.,70.,100.,150.,200.,300.], dtype=float)

    n_netotal, nA, nT, nkap = 9,7,13,23
    fname = 'emiss_table_fe_f64.bin'
    if not os.path.exists(fname):
        return None
    with open(fname,'rb') as f:
        arr = np.fromfile(f, dtype=np.float64)
    arr = arr.reshape((n_netotal,nA,nT,nkap,3), order='F')

    # Build the scaled log interpolator
    #coord_grids = [netotal_grid, A_fe_grid, T_fe_grid, kappa_fe_grid]
    return RegularGridInterpolator(
        (netotal_grid, A_fe_grid, T_fe_grid, kappa_fe_grid),
        arr, method='linear', bounds_error=False, fill_value=np.nan
    )

def load_double_maxwellian_emiss():
    """Double Maxwellian => (n_total, feh, T_ec, T_eh) => shape(9,13,23,12,3)."""
    # 9x13x23x12
    netotal_dm_grid = np.array([1,10,100,1000,2000,2500,3000, 3287,3500], dtype=float)
    feh_dm = np.array([0.000001,0.001469,0.002,0.0025,0.005,0.0075,0.01,0.025,0.05,0.1,0.2,0.35,0.51])
    

    # Tec_dm => 23
    n_Tec_dm_1 =4 
    n_Tec_dm_1p5 =7 
    n_Tec_dm_1_cont =6
    n_Tec_dm_2 = 5
   
    Tec_dm = np.concatenate([
        0.5*np.arange(n_Tec_dm_1) + 1.,   # 4 points
        0.25*np.arange(n_Tec_dm_1p5) + 2.75, # 7 points
        0.5*np.arange(n_Tec_dm_1_cont) + 4.5, # 6 points
        2.0*np.arange(n_Tec_dm_2) + 7.5,  # 5 points
        np.array([20.]) # 1 point
    ])
   
  
    Teh_dm = np.array([1.,249.,260.,270.,300.,350.,400.,450.,500.,550.,587.,601.])
    

    n_netotal, n_feh, n_Tec, n_Teh = 9,13,23,12 #9x13x23x12
    fname = 'emiss_table_dm_f64.bin'
    if not os.path.exists(fname):
        return None
    with open(fname,'rb') as f:
        arr = np.fromfile(f, dtype=np.float64)
    arr = arr.reshape((n_netotal,n_feh,n_Tec,n_Teh,3), order='F')

    # Build the scaled log interpolator
    #coord_grids = [netotal_dm_grid, feh_dm, Tec_dm, Teh_dm]
    return RegularGridInterpolator(
        (netotal_dm_grid, feh_dm, Tec_dm, Teh_dm),
        arr, method='linear', bounds_error=False, fill_value=np.nan
    )

def load_single_maxwellian_emiss():
    """Single Cor Maxwellian => (n_total, T_ec) => shape(9,38,3)."""
    netotal_sm_grid = np.array([1,10,100,1000,2000,2500,3000,3500,4000], dtype=float)
    n_Te_single_max_1 = 26
    n_Te_single_max_2 = 11
    Tec_sm = np.concatenate((np.array([0.99]),0.25*np.arange(n_Te_single_max_1) + 1.,0.5*np.arange(n_Te_single_max_2) + 7.5  ))
    n_Te_sm = 38
    n_netotal = 9
    fname = 'emiss_table_max_f64.bin'
    if not os.path.exists(fname):
        return None
    with open(fname,'rb') as f:
        arr = np.fromfile(f, dtype=np.float64)
    arr = arr.reshape((n_netotal,n_Te_sm,3), order='F')

    # Build the scaled log interpolator
    return RegularGridInterpolator(
        (netotal_sm_grid, Tec_sm),
        arr, method='linear', bounds_error=False, fill_value=np.nan
    )



##############################################################################
#    2) SIMPSON LOS INTEGRATION + HELPER
##############################################################################

def integrate_along_los_scipy(integrand, s_axis, axis=2):
    """
    integrand: shape (..., Ns)
    s_axis: shape(Ns,)
    => simpson => shape(...)
    """
    return simpson(integrand, x=s_axis, axis=axis)

##############################################################################
#    3) MIRROR Z=0..3 => Z=-3..3
##############################################################################

def mirror_z_data(col2d, z_array):
    """
    col2d: shape(nrho, nz) with z_array in [0..3].
    Mirror about z=0 => shape(nrho, 2*nz-1).
    """
    z_neg = -z_array[::-1]   # e.g. if z_array=[0,0.1,...,3], then z_neg= [-3,...,-0.1, -0.0]
    z_neg = z_neg[:-1]       # remove the repeated 0
    z_mir = np.concatenate((z_neg, z_array))  # => shape(2*nz-1,)

    col_neg = col2d[:, ::-1]    # flips the z dimension
    col_neg = col_neg[:, :-1]   # remove repeated row for z=0
    col2d_mir = np.concatenate((col_neg, col2d), axis=1)
    return col2d_mir, z_mir

##############################################################################
#    4) PLOT 3 EMISSIONS (S+, S++, O+) AS 3 STACKED CONTOUR PLOTS
##############################################################################

def plot_emission_maps(rho_vals, z_vals, col0, col1, col2, figtitle, outpdf):
    """
    Single figure, 3 stacked subplots (S^+ top, S^{++} middle, O^+ bottom).
    Mirror z => [-2.5.2.5]. X= rho_vals in [4..10].
    Save as outpdf, tight layout, no extra space.
    """
    col0_mir, z_mir = mirror_z_data(col0, z_vals)
    col1_mir, _     = mirror_z_data(col1, z_vals)
    col2_mir, _     = mirror_z_data(col2, z_vals)

    # We'll do meshgrid(Rg, Zg) => shape(nrho,nzmir)
    # We'll call contourf(rho_vals, z_mir, colX_mir.T).
    RATIO = 1.2
    h_axis  = 3.5    # inches you want per panel in the vertical direction
    w_axis  = h_axis * RATIO        # = 4.2 inches here
    fig_h   = h_axis * 3            # three stacked panels
    fig_w   = w_axis
    
 
    fig, axes = plt.subplots(nrows=3, ncols=1, sharex=True, figsize=(fig_w,fig_h))
    plt.subplots_adjust(hspace=0.02)
    
    
    
    

    species_labels = [r"S$^{+}$ (673.1 nm) Rayleighs",
                      r"S$^{++}$ (68.0 nm) Rayleighs",
                      r"O$^{+}$ (83.3 nm) Rayleighs"]
    col_mir_list   = [col0_mir, col1_mir, col2_mir]
    
    nlevels=1000
    

    contour_levelss = [
        [10,100],        
        [5, 50],        
        [5,50],        
    ]
    rho_mesh0, z_mesh0 =  np.meshgrid(rho_vals, z_mir, indexing='ij')
    rho_mesh, z_mesh = rho_mesh0.T, z_mesh0.T
    
    rho_flat, z_flat =  rho_mesh.ravel(), z_mesh.ravel()
    print(col_mir_list[0].T.shape,rho_mesh.shape)
    
    for i, ax in enumerate(axes):
        ax.set_aspect("equal", adjustable="box")   # 1 data-unit = 1 data-unit
        
        val0 = col_mir_list[i].T
        val = val0.ravel()
        valmax = np.nanmax(val)
        levelss = np.logspace(np.log10(1), np.log10(valmax), nlevels)
        normss = colors.LogNorm(vmin=1, vmax=valmax)
        
        clevelss = contour_levelss[i]
        
        valid = (np.isfinite(val) & (val > 0) )
        rho_valid = rho_flat[valid]
        z_valid = z_flat[valid]
        #density[nvalid] = 0.00001
        val_valid = val[valid]
        
        # Define the grid
        x_min, x_max = 4.0, 10.0
        y_min, y_max = -2.5, 2.5
        num_rho_points = 3601
        num_z_points = 501
        rho_lin = np.linspace(x_min, x_max, num_rho_points)
        z_lin = np.linspace(y_min, y_max, num_z_points)
        rho_grid, z_grid = np.meshgrid(rho_lin, z_lin, indexing='ij')

        # Interpolate the density onto the grid (use absolute value of rho for symmetry)
        points = np.column_stack((rho_valid, z_valid))
        values = val_valid
        val_grid0 = griddata(points, values, (np.abs(rho_grid), z_grid), method='linear')

        # Replace NaN values with zero
        val_grid1 = np.nan_to_num(val_grid0, nan=0.0)
        
        # ⇩⇩ NEW: atmospheric-seeing convolution ⇩⇩
        sigma_RJ = 0.031       # physical σ in R_J
        drho     = rho_lin[1] - rho_lin[0]   # pixel size in ρ
        dz       = z_lin[1] - z_lin[0]       # pixel size in z
        sigma_pix = (sigma_RJ / drho, sigma_RJ / dz)   # (σ_ρ, σ_z) in pixels
        
        # Convolve; 'nearest' avoids edge artefacts, truncate=3 → kernel = 3 σ
        val_grid = ndi.gaussian_filter(val_grid1,
                                       sigma=sigma_pix,
                                       mode='nearest',
                                       truncate=3.0)

        # Compute valid_grid
        valid_grid = ( val_grid >0 & np.isfinite(val_grid))

        # Set density to nan in non-valid areas
        val_grid[~valid_grid] = np.nan#1e-30
        
        contour = ax.contourf(rho_grid, z_grid, val_grid, levels=levelss, norm=normss, cmap='viridis')
        for c in contour.collections:
            c.set_rasterized(True)
        # Define contour levels for black contour lines for each species
        # Add black contour lines with labels
        contour_lines = ax.contour(rho_grid, z_grid, val_grid, levels=clevelss, colors='black', linewidths=1)
        ax.clabel(contour_lines, fmt='%g', inline=True, fontsize=11, colors='black')
        # Manage y-axis tick labels
        
        y_ticks = np.array([-2.5,-2,-1,0,1,2])
        ax.set_yticks(y_ticks)
        ax.tick_params(labelleft=True)
        
            
        for label in ax.get_xticklabels():
            label.set_fontweight('bold')
        for label in ax.get_yticklabels():
            label.set_fontweight('bold')
            
        
        
        cbar = fig.colorbar(contour, ax=ax,
                    ticks=mticker.LogLocator(),         # major ticks 1,10,100…
                    format=mticker.LogFormatterSciNotation(), 
                    pad=0.02,orientation='vertical', 
                    fraction=0.05, )  # 1×10^n labels
        cbar.ax.minorticks_on()      
        cbar.ax.tick_params(labelsize=9)
        

        
        ax.set_xlim(x_min,x_max)
        ax.set_ylim(y_min,y_max)
        
        ax.text(
            0.99, 0.99,
            f"{species_labels[i]}\n{figtitle}",   # first line + newline + second line
            transform=ax.transAxes,
            ha='right', va='top',
            fontsize=10, fontweight='bold',
            linespacing=1.1,        # (optional) vertical spacing between the two lines
            color='black'
        )
        
        ax.set_ylabel(r'z (R$_J$)', fontsize=10, fontweight='bold')

        if i < 2:
            ax.tick_params(labelbottom=False)
        else:
            ax.set_xlabel(r'$\rho_{\mathrm{ceq}}$ (R$_J$)', fontsize=10, fontweight='bold')

        ax.tick_params(axis='both', which='major', labelsize=9)
        
    plt.savefig(outpdf, dpi=600, bbox_inches='tight', pad_inches=0)
    plt.close(fig)      

##############################################################################
#    5) MAIN: LOAD LOS, INTERPOLATE, PLOT
##############################################################################

def main():
    ds = 0.01
    s_1D = np.arange(0.0, 10.0+ds*0.5, ds)  # => shape(1001,)
    factor_1e6 = 1e-6
    symmetry_factor = 2.0

    # Build interpolators
    sk_interp  = load_standard_kappa_emiss()
    pk_interp  = load_product_kappa_emiss()
    fe_interp  = load_fried_egg_emiss()
    dm_interp  = load_double_maxwellian_emiss()
    sm_interp  = load_single_maxwellian_emiss()

    # ---- STANDARD KAPPA ----
    sk_file = "los_results_standard_kappa_with_constraints.npz"
    if os.path.exists(sk_file) and sk_interp is not None:
        print("\n--- Standard Kappa ---")
        data_sk = np.load(sk_file)

        n_elec = data_sk["n_elec"]       # shape(241,53, Ns)
        Tperp_e = data_sk["Tperp_elec"]    # T = Tperp
        kperp_e = data_sk["kperp_elec"]    # kappa = kperp
        n_sp   = data_sk["n_sp"]         # S^+
        n_s2p  = data_sk["n_s2p"]        # S^{++}
        n_op   = data_sk["n_op"]         # O^{+}

        rho1D  = data_sk["rho1D"]        # shape(241,)
        z1D    = data_sk["z1D"]          # shape(53,)

        Nx, Nz, Ns_ = n_elec.shape

        # Flatten
        n_e_flat = n_elec.ravel()
        T_flat   = Tperp_e.ravel()
        k_flat   = kperp_e.ravel()

        # Interpolate => shape(N,3)
        if sk_interp is not None:
            points_sk = np.column_stack((n_e_flat, T_flat, k_flat))
            emiss_sk_flat = sk_interp(points_sk)  # => shape(N,3)

            # Reshape => (241,53,Ns,3)
            emiss_4D = emiss_sk_flat.reshape((Nx,Nz,Ns_,3))

            # Multiply by densities
            eps0_Sp  = emiss_4D[...,0]* n_sp
            eps1_S2p = emiss_4D[...,1]* n_s2p
            eps2_Op  = emiss_4D[...,2]* n_op
            bad0 = (~np.isfinite(eps0_Sp)) | (eps0_Sp < 1.1e-30) 
            bad1 = (~np.isfinite(eps1_S2p)) | (eps1_S2p < 1.1e-30)
            bad2 = (~np.isfinite(eps2_Op)) | (eps2_Op < 1.1e-30)
            eps0_Sp[bad0] = 0.0
            eps1_S2p[bad1] = 0.0
            eps2_Op[bad2] = 0.0
            RJcm = 7.1492e9 # 1 RJ in cm
            col0 = integrate_along_los_scipy(eps0_Sp, s_1D, axis=2)*factor_1e6*symmetry_factor*RJcm
            col1 = integrate_along_los_scipy(eps1_S2p,s_1D, axis=2)*factor_1e6*symmetry_factor*RJcm
            col2 = integrate_along_los_scipy(eps2_Op, s_1D, axis=2)*factor_1e6*symmetry_factor*RJcm

            plot_emission_maps(rho1D, z1D, col0, col1, col2,
                               figtitle="Standard Kappa",
                               outpdf="standard_kappa_emission.pdf")

    # ---- PRODUCT KAPPA ----
    pk_file = "los_results_product_kappa_with_constraints.npz"
    if os.path.exists(pk_file) and pk_interp is not None:
        print("\n--- Product Kappa ---")
        data_pk = np.load(pk_file)

        n_elec   = data_pk["n_elec"]       # shape(241,53,Ns)
        Tpar_e   = data_pk["Tpar_elec"]    # for A => Tperp / Tpar
        Tperp_e  = data_pk["Tperp_elec"]
        kpar_e   = data_pk["kpar_elec"]    # for lambda => kperp / kpar
        kperp_e  = data_pk["kperp_elec"]
        n_sp     = data_pk["n_sp"]
        n_s2p    = data_pk["n_s2p"]
        n_op     = data_pk["n_op"]
        rho1D    = data_pk["rho1D"]
        z1D      = data_pk["z1D"]

        Nx, Nz, Ns_ = n_elec.shape

        # Flatten
        n_e_flat    = n_elec.ravel()
        Tpar_flat   = Tpar_e.ravel()
        Tperp_flat  = Tperp_e.ravel()
        kpar_flat   = kpar_e.ravel()
        kperp_flat  = kperp_e.ravel()

        # Compute A = T_perp / T_par, lambda = kperp / kpar
        # watch for zero Tpar or kpar
        A_arr   = np.zeros_like(Tpar_flat)
        lam_arr = np.zeros_like(kpar_flat)
        valid_tpar = (Tpar_flat>0)
        valid_kpar = (kpar_flat>0)
        A_arr[valid_tpar]   = Tperp_flat[valid_tpar] / Tpar_flat[valid_tpar]
        lam_arr[valid_kpar] = kperp_flat[valid_kpar] / kpar_flat[valid_kpar]

        # Interpolate => (n_e, A, Tperp, lam, kperp)
        points_pk = np.column_stack((n_e_flat, A_arr, Tperp_flat, lam_arr, kperp_flat))
        emiss_pk_flat = pk_interp(points_pk)  # => shape(N,3)

        emiss_4D = emiss_pk_flat.reshape((Nx,Nz,Ns_,3))

        eps0_Sp  = emiss_4D[...,0]* n_sp
        eps1_S2p = emiss_4D[...,1]* n_s2p
        eps2_Op  = emiss_4D[...,2]* n_op
        bad0 = (~np.isfinite(eps0_Sp)) | (eps0_Sp < 1.1e-30) 
        bad1 = (~np.isfinite(eps1_S2p)) | (eps1_S2p < 1.1e-30)
        bad2 = (~np.isfinite(eps2_Op)) | (eps2_Op < 1.1e-30)
        eps0_Sp[bad0] = 0.0
        eps1_S2p[bad1] = 0.0
        eps2_Op[bad2] = 0.0
        RJcm = 7.1492e9 # 1 RJ in cm
        col0 = simpson(eps0_Sp, x=s_1D, axis=2)*factor_1e6*symmetry_factor*RJcm
        col1 = simpson(eps1_S2p,x=s_1D, axis=2)*factor_1e6*symmetry_factor*RJcm
        col2 = simpson(eps2_Op, x=s_1D, axis=2)*factor_1e6*symmetry_factor*RJcm

        plot_emission_maps(rho1D, z1D, col0, col1, col2,
                           figtitle="Product Kappa",
                           outpdf="product_kappa_emission.pdf")

    # ---- FRIED-EGG ----
    fe_file = "los_results_fried_egg_with_constraints.npz"
    if os.path.exists(fe_file) and fe_interp is not None:
        print("\n--- Fried-Egg ---")
        data_fe = np.load(fe_file)

        n_elec  = data_fe["n_elec"]
        Tpar_e  = data_fe["Tpar_elec"]   # to compute A => Tperp / Tpar
        Tperp_e = data_fe["Tperp_elec"]
        kperp_e = data_fe["kperp_elec"]  # only kperp => "kappa"
        n_sp    = data_fe["n_sp"]
        n_s2p   = data_fe["n_s2p"]
        n_op    = data_fe["n_op"]
        rho1D   = data_fe["rho1D"]
        z1D     = data_fe["z1D"]

        Nx, Nz, Ns_ = n_elec.shape

        # Flatten
        n_e_flat   = n_elec.ravel()
        Tpar_flat  = Tpar_e.ravel()
        Tperp_flat = Tperp_e.ravel()
        kperp_flat = kperp_e.ravel()

        # A = Tperp / Tpar
        A_arr = np.zeros_like(Tpar_flat)
        valid_tpar = (Tpar_flat>0)
        A_arr[valid_tpar] = Tperp_flat[valid_tpar]/ Tpar_flat[valid_tpar]

        # Interpolate => (n_e, A, Tperp, kperp)
        points_fe = np.column_stack((n_e_flat, A_arr, Tperp_flat, kperp_flat))
        emiss_fe_flat = fe_interp(points_fe)  # => shape(N,3)
        emiss_4D = emiss_fe_flat.reshape((Nx,Nz,Ns_,3))

        eps0_Sp  = emiss_4D[...,0]* n_sp
        eps1_S2p = emiss_4D[...,1]* n_s2p
        eps2_Op  = emiss_4D[...,2]* n_op
        
        bad0 = (~np.isfinite(eps0_Sp)) | (eps0_Sp < 1.1e-30) 
        bad1 = (~np.isfinite(eps1_S2p)) | (eps1_S2p < 1.1e-30)
        bad2 = (~np.isfinite(eps2_Op)) | (eps2_Op < 1.1e-30)
        eps0_Sp[bad0] = 0.0
        eps1_S2p[bad1] = 0.0
        eps2_Op[bad2] = 0.0
        
        RJcm = 7.1492e9 # 1 RJ in cm
        col0 = simpson(eps0_Sp, x=s_1D, axis=2)*factor_1e6*symmetry_factor*RJcm
        col1 = simpson(eps1_S2p,x=s_1D, axis=2)*factor_1e6*symmetry_factor*RJcm
        col2 = simpson(eps2_Op, x=s_1D, axis=2)*factor_1e6*symmetry_factor*RJcm

        plot_emission_maps(rho1D, z1D, col0, col1, col2,
                           figtitle="Fried-Egg",
                           outpdf="fried_egg_emission.pdf")

    # ---- DOUBLE MAXWELLIAN ----
    dm_file = "los_results_maxwellian_with_constraints.npz"
    dm_file2 ="los_results_netot_feh_maxwellian.npz"
    if os.path.exists(dm_file) and dm_interp is not None:
        print("\n--- Double Maxwellian (DM) ---")
        data_dm = np.load(dm_file)
        data_dm2 = np.load(dm_file2)

        netotal = data_dm2["n_elec_tot"] # core+hot=total double maxwellian electron density
        feh = data_dm2["feh"] # feh=neh/netotal
        T_ec = data_dm["Tperp_elec"] # core single maxwellian electron Temperature along LOS
        T_eh = data_dm["Tperp_eh"]
        n_sp   = data_dm["n_sp"]
        n_s2p  = data_dm["n_s2p"]
        n_op   = data_dm["n_op"]
        rho1D  = data_dm["rho1D"]
        z1D    = data_dm["z1D"]

        Nx, Nz, Ns_ = netotal.shape
        #print(netotal.shape)
        
        
        

        #feh_z0 = feh[:,0,:]
        #feh_z0_flat = feh_z0.ravel()
        netotal_flat = netotal.flatten()
        feh_flat     = feh.flatten()
        Tec_flat     = T_ec.flatten()
        Teh_flat     = T_eh.flatten()
        
        
        
        #print("DM netotal max, min = ",np.nanmax(netotal_flat),np.nanmin(netotal_flat))
        #print("DM feh max, min = ",np.nanmax(feh_flat),np.nanmin(feh_flat[feh_flat > 1.1e-30 ]))
        #print("DM Tec max, min = ",np.nanmax(Tec_flat),np.nanmin(Tec_flat))
        #print("DM Teh max, min = ",np.nanmax(Teh_flat),np.nanmin(Teh_flat))
        
        #print("DM feh at z =0  max, min = ",np.nanmax(feh_z0_flat),np.nanmin(feh_z0_flat[ (feh_z0_flat > 1.1e-30) ]))
        feh_modified = feh_flat.copy()
        cond_mask = (feh_modified < 0.001469) #(feh_modified >= 0.0015) & (feh_modified < 0.002)
        feh_modified[cond_mask] = 0.000001 #0.002
        #mask_dm = (feh_modified >= 0.002)
        #mask_sm = ~mask_dm  # i.e., feh_modified < 0.002
        
       
        
        #points_dm = np.column_stack((netotal_flat, feh_flat, Tec_flat, Teh_flat))
        points_dm = np.column_stack((netotal_flat, feh_modified, Tec_flat, Teh_flat))
        emiss_dm_flat = dm_interp(points_dm)  # => shape(N,3)
        emiss_4D = emiss_dm_flat.reshape((Nx,Nz,Ns_,3))
        
        """
        # --- Double Maxwellian Emission ---
        if mask_dm.any():
            points_dm = np.column_stack((
                netotal_flat[mask_dm],
                feh_modified[mask_dm],
                Tec_flat[mask_dm],
                Teh_flat[mask_dm]
            ))
            emiss_dm_flat = dm_interp(points_dm)  # shape => (N_dm, 3)
        #else:
           # emiss_dm_flat = np.zeros((0,3), dtype=float)
    
        # --- Single Maxwellian Emission ---
        if mask_sm.any():
            # Single Maxwellian only needs (netotal_sm, Tec_sm)
            points_sm = np.column_stack((
                netotal_flat[mask_sm],
                Tec_flat[mask_sm]
            ))
            emiss_sm_flat = sm_interp(points_sm)  # shape => (N_sm, 3)
        #else:
            #emiss_sm_flat = np.zeros((0,3), dtype=float)
    
        # -------------------------------------------------------------------------
        # 3) Combine (dm + sm) emission values back into one array of shape (N, 3).
        #    We must place the results in their original positions.
        # -------------------------------------------------------------------------
        # Preallocate total array
        emiss_tot_flat = np.empty((netotal_flat.size, 3), dtype=float)
        
        # DM points in their original order
        emiss_tot_flat[mask_dm] = emiss_dm_flat
        # SM points in their original order
        emiss_tot_flat[mask_sm] = emiss_sm_flat
    
        # Finally reshape to (Nx, Nz, Ns_, 3)
        emiss_4D = emiss_tot_flat.reshape((Nx, Nz, Ns_, 3))
        """
        eps0_Sp  = emiss_4D[...,0]* n_sp
        eps1_S2p = emiss_4D[...,1]* n_s2p
        eps2_Op  = emiss_4D[...,2]* n_op

        

        bad0 = (~np.isfinite(eps0_Sp)) | (eps0_Sp < 1.1e-30) 
        bad1 = (~np.isfinite(eps1_S2p)) | (eps1_S2p < 1.1e-30)
        bad2 = (~np.isfinite(eps2_Op)) | (eps2_Op < 1.1e-30)
        eps0_Sp[bad0] = 0.0
        eps1_S2p[bad1] = 0.0
        eps2_Op[bad2] = 0.0




        RJcm = 7.1492e9 # 1 RJ in cm
        col0 = simpson(eps0_Sp, x=s_1D, axis=2)*factor_1e6*symmetry_factor*RJcm
        col1 = simpson(eps1_S2p,x=s_1D, axis=2)*factor_1e6*symmetry_factor*RJcm
        col2 = simpson(eps2_Op, x=s_1D, axis=2)*factor_1e6*symmetry_factor*RJcm
        
        #print("col_0",np.nanmin(col0),np.nanmax(col0))
        #print("col_1",np.nanmin(col1),np.nanmax(col1))
        #print("col_2",np.nanmin(col2),np.nanmax(col2))
        plot_emission_maps(rho1D, z1D, col0, col1, col2,
                           figtitle="Double Maxwellian",
                           outpdf="double_maxwellian_emission.pdf")

    print("\nDone all cases.")

if __name__=="__main__":
    main()
