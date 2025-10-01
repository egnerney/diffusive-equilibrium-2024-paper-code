#!/usr/bin/env python3
"""
NAME:
    sim_citep_2_diff_eq.py

PURPOSE:
    Calculate UV emission spectra from the Io Plasma Torus using a 3D
    diffusive equilibrium model aligned with Jupiter's magnetic field lines.
    Integrates emission along arbitrary lines of sight through the torus.

DESCRIPTION:
    This code simulates UV emission (550-2100 Å) using a physically realistic
    3D torus model based on diffusive equilibrium along JRM33+CON2020 magnetic
    field lines. The emission calculation properly integrates through the
    non-axisymmetric torus structure:

        B_λ = 10^(-6) × ∫ ε_λ(T_e(s), n_e(s)) × n_ion(s) ds

    where the plasma parameters vary along field lines according to
    diffusive equilibrium with centrifugal and pressure gradient forces.

MODEL STRUCTURE:
    The torus model is defined on an irregular grid following field lines:
    - 501 field lines crossing equator at ρ = 5.00 to 10.00 R_J (0.01 R_J steps)
    - 360 azimuthal positions φ = 0° to 359° (1° steps)
    - 1001 points along each field line from λ_III = -50° to +50° (0.1° steps)
    - System III coordinates (not centrifugal equator)
    - Field-aligned diffusive equilibrium determines density stratification

REQUIRED INPUT FILES:
    Either the individual model files:
        - nelec_out_mymodel1_*.txt  - Electron density [cm^-3]
        - nsp_out_mymodel1_*.txt    - S+ density [cm^-3]
        - ns2p_out_mymodel1_*.txt   - S++ density [cm^-3]
        - ns3p_out_mymodel1_*.txt   - S+++ density [cm^-3]
        - nop_out_mymodel1_*.txt    - O+ density [cm^-3]
        - no2p_out_mymodel1_*.txt   - O++ density [cm^-3]
        - Telec_out_mymodel1_*.txt  - Electron temperature [eV]
        - Tic_out_mymodel1_*.txt    - Ion core temperature [eV]
        - x_out*.txt, y_out*.txt, z_out*.txt - Positions [R_J]

    Or the pre-processed pickle file:
        - variables_rebin3D.pkl containing all arrays in 501x360x1001 format

    Plus CHIANTI emission tables:
        - CHIANTI_11.0.2_emiss_arrays_all_species_all_wavelengths_50x50_logspaced.npz

COORDINATE SYSTEM:
    System III Jovicentric coordinates:
    - X, Y in Jovian equatorial plane
    - Z along rotation axis
    - ρ = sqrt(X² + Y²) cylindrical radius
    - φ = atan(Y/X) east longitude
    - λ_III = asin(Z/r) System III latitude

KEY ASSUMPTIONS:
    - Diffusive equilibrium along magnetic field lines
    - JRM33+CON2020 magnetic field model
    - Optically thin emission
    - Single Maxwellian electron distribution
    - Steady-state torus (no temporal variations)
    - CHIANTI coronal approximation for atomic physics

USAGE:
    python sim_citep_2_diff_eq.py

OUTPUT:
    - UV spectra for different viewing geometries [Rayleighs/Å]
    - 2D emission maps
    - Diagnostic plots showing torus structure

AUTHOR:
    CITEP Team - Community Io Torus Emission Package
    For the Magnetospheres of Outer Planets (MOP) community
    Python version created for open community access

HISTORY:
    2024 - Initial release with field-aligned diffusive equilibrium model
    2025 - Python conversion for broader accessibility
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf
from scipy.interpolate import RegularGridInterpolator
import os
import pickle
from typing import Dict, Tuple, Optional, Any
import warnings

# Constants
RJ_CM = 7.1492e9  # Jupiter radius in cm

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def load_3d_torus_model() -> Dict[str, np.ndarray]:
    """
    Load the 3D field-aligned torus model.
    
    Returns:
        Dictionary with all plasma parameters on 501x360x1001 grid
    """
    
    # Check if pre-processed pickle file exists
    if os.path.exists('variables_rebin3D.pkl'):
        print('Loading pre-processed 3D model from variables_rebin3D.pkl...')
        with open('variables_rebin3D.pkl', 'rb') as f:
            model = pickle.load(f)
        print('Model loaded: 501 x 360 x 1001 grid')
        return model
    
    else:
        # Load from individual text files
        print('Loading 3D model from text files...')
        print('This may take several minutes...')
        
        # Load position arrays
        print('Loading position arrays...')
        x = np.loadtxt('x_outlat-50to50_0.1deginterps_501x360_3D_analytic_180360x1001.txt')
        y = np.loadtxt('y_outlat-50to50_0.1deginterps_501x360_3D_analytic_180360x1001.txt')
        z = np.loadtxt('z_outlat-50to50_0.1deginterps_501x360_3D_analytic_180360x1001.txt')
        
        # Load density arrays
        print('Loading density arrays...')
        nel = np.loadtxt('nelec_out_mymodel1_diffeq_jrm33+Con2020_0fillslat-50to50_0.1deginterps_501x360_3D_analytic_180360x1001.txt')
        nsp = np.loadtxt('nsp_out_mymodel1_diffeq_jrm33+Con2020_0fillslat-50to50_0.1deginterps_501x360_3D_analytic_180360x1001.txt')
        ns2p = np.loadtxt('ns2p_out_mymodel1_diffeq_jrm33+Con2020_0fillslat-50to50_0.1deginterps_501x360_3D_analytic_180360x1001.txt')
        ns3p = np.loadtxt('ns3p_out_mymodel1_diffeq_jrm33+Con2020_0fillslat-50to50_0.1deginterps_501x360_3D_analytic_180360x1001.txt')
        nop = np.loadtxt('nop_out_mymodel1_diffeq_jrm33+Con2020_0fillslat-50to50_0.1deginterps_501x360_3D_analytic_180360x1001.txt')
        no2p = np.loadtxt('no2p_out_mymodel1_diffeq_jrm33+Con2020_0fillslat-50to50_0.1deginterps_501x360_3D_analytic_180360x1001.txt')
        
        # Load temperature arrays
        print('Loading temperature arrays...')
        tec = np.loadtxt('Telec_out_mymodel1_diffeq_jrm33+Con2020_0fillslat-50to50_0.1deginterps_501x360_3D_analytic_180360x1001.txt')
        tico = np.loadtxt('Tic_out_mymodel1_diffeq_jrm33+Con2020_0fillslat-50to50_0.1deginterps_501x360_3D_analytic_180360x1001.txt')
        
        # Calculate derived quantities
        rho = np.sqrt(x**2 + y**2)
        r = np.sqrt(x**2 + y**2 + z**2)
        lat = np.degrees(np.arcsin(z/r))
        phi = np.degrees(np.arctan2(y, x))
        
        # Reshape from 2D (180360x1001) to 3D (501x360x1001)
        print('Reshaping arrays to 3D format...')
        
        nrfls = 501
        nphifls = 360
        nlat = 1001
        
        x_rebin = np.zeros((nrfls, nphifls, nlat))
        y_rebin = np.zeros((nrfls, nphifls, nlat))
        z_rebin = np.zeros((nrfls, nphifls, nlat))
        r_rebin = np.zeros((nrfls, nphifls, nlat))
        rho_rebin = np.zeros((nrfls, nphifls, nlat))
        phi_rebin = np.zeros((nrfls, nphifls, nlat))
        lat_rebin = np.zeros((nrfls, nphifls, nlat))
        ne_rebin = np.zeros((nrfls, nphifls, nlat))
        nsp_rebin = np.zeros((nrfls, nphifls, nlat))
        ns2p_rebin = np.zeros((nrfls, nphifls, nlat))
        ns3p_rebin = np.zeros((nrfls, nphifls, nlat))
        nop_rebin = np.zeros((nrfls, nphifls, nlat))
        no2p_rebin = np.zeros((nrfls, nphifls, nlat))
        Tec_rebin = np.zeros((nrfls, nphifls, nlat))
        Tic_rebin = np.zeros((nrfls, nphifls, nlat))
        
        # Mapping: 2D index = 360*i + j maps to 3D[i,j,:]
        for i in range(501):
            for j in range(360):
                idx_2d = 360 * i + j
                x_rebin[i, j, :] = x[idx_2d, :]
                y_rebin[i, j, :] = y[idx_2d, :]
                z_rebin[i, j, :] = z[idx_2d, :]
                r_rebin[i, j, :] = r[idx_2d, :]
                rho_rebin[i, j, :] = rho[idx_2d, :]
                phi_rebin[i, j, :] = phi[idx_2d, :]
                lat_rebin[i, j, :] = lat[idx_2d, :]
                ne_rebin[i, j, :] = nel[idx_2d, :]
                nsp_rebin[i, j, :] = nsp[idx_2d, :]
                ns2p_rebin[i, j, :] = ns2p[idx_2d, :]
                ns3p_rebin[i, j, :] = ns3p[idx_2d, :]
                nop_rebin[i, j, :] = nop[idx_2d, :]
                no2p_rebin[i, j, :] = no2p[idx_2d, :]
                Tec_rebin[i, j, :] = tec[idx_2d, :]
                Tic_rebin[i, j, :] = tico[idx_2d, :]
        
        # Assume S4+ is 10% of S3+ if not provided
        ns4p_rebin = 0.1 * ns3p_rebin
        
        model = {
            'x': x_rebin, 'y': y_rebin, 'z': z_rebin,
            'r': r_rebin, 'rho': rho_rebin, 'phi': phi_rebin, 'lat': lat_rebin,
            'nel': ne_rebin, 'nsp': nsp_rebin, 'ns2p': ns2p_rebin, 'ns3p': ns3p_rebin,
            'ns4p': ns4p_rebin, 'nop': nop_rebin, 'no2p': no2p_rebin,
            'Te': Tec_rebin, 'Ti': Tic_rebin
        }
        
        print('Model loaded and reshaped: 501 x 360 x 1001 grid')
        
        # Save for faster loading next time
        print('Saving processed model to variables_rebin3D.pkl...')
        with open('variables_rebin3D.pkl', 'wb') as f:
            pickle.dump(model, f)
        
        return model


def interpolate_field_model(x_point: float, y_point: float, z_point: float, 
                           model: Dict[str, np.ndarray]) -> Dict[str, Any]:
    """
    Interpolate plasma parameters from field-aligned model to arbitrary point.
    
    Args:
        x_point, y_point, z_point: Position to interpolate to [R_J]
        model: Dictionary with 3D field-aligned model
    
    Returns:
        Dictionary with interpolated plasma parameters
    
    Method:
        1. Find nearest field line in (ρ, φ) space
        2. Interpolate along that field line in latitude
        3. Use bilinear interpolation between adjacent field lines
    """
    
    # Calculate cylindrical coordinates of query point
    rho_point = np.sqrt(x_point**2 + y_point**2)
    phi_point = np.degrees(np.arctan2(y_point, x_point))
    if phi_point < 0:
        phi_point += 360.0
    r_point = np.sqrt(x_point**2 + y_point**2 + z_point**2)
    if r_point > 0:
        lat_point = np.degrees(np.arcsin(z_point/r_point))
    else:
        lat_point = 0.0
    
    # Check if point is in model domain
    if rho_point < 5.0 or rho_point > 10.0:
        # Outside radial range
        return {'nel': 0.0, 'nsp': 0.0, 'ns2p': 0.0, 'ns3p': 0.0, 'ns4p': 0.0,
                'nop': 0.0, 'no2p': 0.0, 'Te': 0.0, 'valid': False}
    
    if np.abs(lat_point) > 50.0:
        # Outside latitude range
        return {'nel': 0.0, 'nsp': 0.0, 'ns2p': 0.0, 'ns3p': 0.0, 'ns4p': 0.0,
                'nop': 0.0, 'no2p': 0.0, 'Te': 0.0, 'valid': False}
    
    # Find indices for interpolation
    # Radial index (0-500 for 5.00-10.00 RJ)
    rho_idx = (rho_point - 5.0) / 0.01
    i_rho = int(np.floor(rho_idx))
    if i_rho >= 500:
        i_rho = 499
    w_rho = rho_idx - i_rho
    
    # Azimuthal index (0-359 for 0-359 degrees)
    phi_idx = phi_point
    i_phi = int(np.floor(phi_idx))
    if i_phi >= 360:
        i_phi = 359
    w_phi = phi_idx - i_phi
    
    # Handle azimuthal wrap-around
    i_phi_next = i_phi + 1
    if i_phi_next >= 360:
        i_phi_next = 0
    
    # Latitude index (0-1000 for -50 to +50 degrees)
    lat_idx = (lat_point + 50.0) / 0.1
    i_lat = int(np.floor(lat_idx))
    if i_lat >= 1000:
        i_lat = 999
    if i_lat < 0:
        i_lat = 0
    w_lat = lat_idx - i_lat
    
    # Bilinear interpolation in (rho, phi) and linear in latitude
    ne_interp = 0.0
    nsp_interp = 0.0
    ns2p_interp = 0.0
    ns3p_interp = 0.0
    ns4p_interp = 0.0
    nop_interp = 0.0
    no2p_interp = 0.0
    Te_interp = 0.0
    
    # Loop over the 8 corners
    for di_rho in range(2):
        for di_phi in range(2):
            for di_lat in range(2):
                # Calculate indices
                i_r = min(i_rho + di_rho, 500)
                i_p = i_phi if di_phi == 0 else i_phi_next
                i_l = min(i_lat + di_lat, 1000)
                
                # Calculate weight
                w_r = (1 - w_rho) if di_rho == 0 else w_rho
                w_p = (1 - w_phi) if di_phi == 0 else w_phi
                w_l = (1 - w_lat) if di_lat == 0 else w_lat
                weight = w_r * w_p * w_l
                
                # Accumulate weighted values
                ne_interp += weight * model['nel'][i_r, i_p, i_l]
                nsp_interp += weight * model['nsp'][i_r, i_p, i_l]
                ns2p_interp += weight * model['ns2p'][i_r, i_p, i_l]
                ns3p_interp += weight * model['ns3p'][i_r, i_p, i_l]
                ns4p_interp += weight * model['ns4p'][i_r, i_p, i_l]
                nop_interp += weight * model['nop'][i_r, i_p, i_l]
                no2p_interp += weight * model['no2p'][i_r, i_p, i_l]
                Te_interp += weight * model['Te'][i_r, i_p, i_l]
    
    return {'nel': ne_interp, 'nsp': nsp_interp, 'ns2p': ns2p_interp,
            'ns3p': ns3p_interp, 'ns4p': ns4p_interp,
            'nop': nop_interp, 'no2p': no2p_interp, 'Te': Te_interp, 'valid': True}


# =============================================================================
# INTERPOLATION AND CONVOLUTION FUNCTIONS
# =============================================================================

def interpolate_emissivity_2D(temp_eV: float, dens_cm3: float,
                              temp_arr: np.ndarray, dens_arr: np.ndarray,
                              emiss_table: np.ndarray) -> np.ndarray:
    """
    Interpolate emissivity from 2D temperature-density table.
    
    Args:
        temp_eV: Electron temperature in eV
        dens_cm3: Electron density in cm^-3
        temp_arr: Temperature grid points
        dens_arr: Density grid points
        emiss_table: 3D emissivity table (n_temp x n_dens x n_lines)
    
    Returns:
        Interpolated emissivities for all spectral lines
    """
    
    log_temp = np.log10(temp_eV)
    log_dens = np.log10(dens_cm3)
    log_temp_arr = np.log10(temp_arr)
    log_dens_arr = np.log10(dens_arr)
    
    n_temp = len(temp_arr)
    n_dens = len(dens_arr)
    n_lines = emiss_table.shape[2]
    
    # Find indices for interpolation
    it = np.searchsorted(log_temp_arr, log_temp) - 1
    id = np.searchsorted(log_dens_arr, log_dens) - 1
    
    # Bound checking
    it = max(0, min(it, n_temp - 2))
    id = max(0, min(id, n_dens - 2))
    
    # Calculate weights
    wt = (log_temp - log_temp_arr[it]) / (log_temp_arr[it+1] - log_temp_arr[it])
    wd = (log_dens - log_dens_arr[id]) / (log_dens_arr[id+1] - log_dens_arr[id])
    
    # Bilinear interpolation
    emiss_interp = ((1-wt) * (1-wd) * emiss_table[it, id, :] +
                    (1-wt) * wd * emiss_table[it, id+1, :] +
                    wt * (1-wd) * emiss_table[it+1, id, :] +
                    wt * wd * emiss_table[it+1, id+1, :])
    
    return emiss_interp


def simulate_IPT_spectrum_Rayleighs_ERF_form(x: np.ndarray, spec_binsize: float,
                                            xwavi: np.ndarray, yptsi: np.ndarray,
                                            fwhm: float) -> np.ndarray:
    """
    Convolve line emission with instrument response using error function.
    
    Args:
        x: Wavelength grid for output spectrum
        spec_binsize: Spectral bin size
        xwavi: Wavelengths of emission lines
        yptsi: Intensities of emission lines
        fwhm: Full width at half maximum of instrument
    
    Returns:
        Convolved spectrum
    """
    
    rootc = 2.0 * np.sqrt(np.log(2.0)) / fwhm
    ypts = np.zeros_like(x)
    
    for i in range(len(xwavi)):
        if yptsi[i] <= 0:
            continue
        ypts += yptsi[i] * 0.5 * (
            erf((x - xwavi[i] + spec_binsize/2.0) * rootc) -
            erf((x - xwavi[i] - spec_binsize/2.0) * rootc))
    
    ypts /= spec_binsize
    return ypts


# =============================================================================
# RAY TRACING FUNCTIONS
# =============================================================================

def ray_box_intersection(ray_origin: np.ndarray, ray_dir: np.ndarray,
                        box_min: np.ndarray, box_max: np.ndarray) -> Tuple[float, float]:
    """
    Calculate ray intersection with axis-aligned bounding box.
    
    Args:
        ray_origin: Starting point of ray
        ray_dir: Direction vector of ray (normalized)
        box_min: Minimum corner of box
        box_max: Maximum corner of box
    
    Returns:
        Tuple of (t_min, t_max) parameters, or (-1, -1) if no intersection
    """
    
    t_min = -1e30
    t_max = 1e30
    
    for i in range(3):
        if np.abs(ray_dir[i]) > 1e-10:
            t1 = (box_min[i] - ray_origin[i]) / ray_dir[i]
            t2 = (box_max[i] - ray_origin[i]) / ray_dir[i]
            
            if t1 > t2:
                t1, t2 = t2, t1
            
            t_min = max(t_min, t1)
            t_max = min(t_max, t2)
            
            if t_min > t_max:
                return -1.0, -1.0
        else:
            if ray_origin[i] < box_min[i] or ray_origin[i] > box_max[i]:
                return -1.0, -1.0
    
    if t_max < 0:
        return -1.0, -1.0
    if t_min < 0:
        t_min = 0.0
    
    return t_min, t_max


def sample_ray_through_torus(ray_origin: np.ndarray, ray_dir: np.ndarray,
                            ds: float) -> Dict[str, Any]:
    """
    Generate sample points along ray through torus volume.
    For field-aligned model, use 5-10 RJ radial range, ±3 RJ vertical.
    
    Args:
        ray_origin: Starting point of ray [R_J]
        ray_dir: Direction vector of ray (normalized)
        ds: Step size along ray [R_J]
    
    Returns:
        Dictionary with sample points and validity flags
    """
    
    torus_rmin = 5.0
    torus_rmax = 10.0
    torus_zmax = 3.0
    
    box_min = np.array([-torus_rmax, -torus_rmax, -torus_zmax])
    box_max = np.array([torus_rmax, torus_rmax, torus_zmax])
    
    t_range = ray_box_intersection(ray_origin, ray_dir, box_min, box_max)
    
    if t_range[0] == -1:
        return {'s': np.array([0.0]), 'x': np.array([ray_origin[0]]),
                'y': np.array([ray_origin[1]]), 'z': np.array([ray_origin[2]]),
                'valid': np.array([False]), 'n_valid': 0}
    
    n_samples = int(np.ceil((t_range[1] - t_range[0]) / ds)) + 1
    s_array = t_range[0] + ds * np.arange(n_samples)
    
    if s_array[-1] < t_range[1]:
        s_array = np.append(s_array, t_range[1])
        n_samples += 1
    
    x_array = ray_origin[0] + s_array * ray_dir[0]
    y_array = ray_origin[1] + s_array * ray_dir[1]
    z_array = ray_origin[2] + s_array * ray_dir[2]
    
    # For field model, valid points are those within model domain
    rho_array = np.sqrt(x_array**2 + y_array**2)
    r_array = np.sqrt(x_array**2 + y_array**2 + z_array**2)
    with np.errstate(invalid='ignore'):
        lat_array = np.degrees(np.arcsin(z_array/r_array))
    lat_array = np.nan_to_num(lat_array)
    
    valid = ((rho_array >= torus_rmin) & (rho_array <= torus_rmax) &
             (np.abs(lat_array) <= 50.0))  # Within ±50° latitude coverage
    n_valid = np.sum(valid)
    
    return {'s': s_array, 'x': x_array, 'y': y_array, 'z': z_array,
            'valid': valid, 'n_valid': n_valid}


def integrate_emission_along_los_field_model(ray_samples: Dict[str, Any],
                                            field_model: Dict[str, np.ndarray],
                                            xwav: np.ndarray, spec_binsize: float,
                                            fwhm: float, temp_arr: np.ndarray,
                                            dens_arr: np.ndarray,
                                            xwavi_struct: Dict[str, np.ndarray],
                                            yptsi_struct: Dict[str, np.ndarray]) -> np.ndarray:
    """
    Integrate emission along line of sight using field-aligned model.
    
    Args:
        ray_samples: Dictionary with sample points along ray
        field_model: Dictionary with 3D field-aligned model
        xwav: Output wavelength grid
        spec_binsize: Spectral bin size
        fwhm: Instrument FWHM
        temp_arr: Temperature grid for emissivity tables
        dens_arr: Density grid for emissivity tables
        xwavi_struct: Dictionary with wavelengths for each ion
        yptsi_struct: Dictionary with emissivities for each ion
    
    Returns:
        Spectrum in Rayleighs/Angstrom
    """
    
    valid_idx = np.where(ray_samples['valid'])[0]
    if len(valid_idx) == 0:
        return np.zeros(len(xwav))
    
    n_valid = len(valid_idx)
    
    # Get plasma parameters at each point by interpolating field model
    nel_array = np.zeros(n_valid)
    nsp_array = np.zeros(n_valid)
    ns2p_array = np.zeros(n_valid)
    ns3p_array = np.zeros(n_valid)
    ns4p_array = np.zeros(n_valid)
    nop_array = np.zeros(n_valid)
    no2p_array = np.zeros(n_valid)
    Tel_array = np.zeros(n_valid)
    
    for i in range(n_valid):
        idx = valid_idx[i]
        
        # Interpolate from field model
        params = interpolate_field_model(ray_samples['x'][idx],
                                        ray_samples['y'][idx],
                                        ray_samples['z'][idx],
                                        field_model)
        
        if params['valid']:
            nel_array[i] = params['nel']
            nsp_array[i] = params['nsp']
            ns2p_array[i] = params['ns2p']
            ns3p_array[i] = params['ns3p']
            ns4p_array[i] = params['ns4p']
            nop_array[i] = params['nop']
            no2p_array[i] = params['no2p']
            Tel_array[i] = params['Te']
    
    # Filter out points with no plasma
    good_idx = np.where((nel_array > 1.0) & (Tel_array > 0.1))[0]
    if len(good_idx) == 0:
        return np.zeros(len(xwav))
    
    n_good = len(good_idx)
    
    # Calculate path elements
    ds_array = np.zeros(n_good)
    
    for i in range(n_good):
        if i == 0:
            if n_good > 1:
                ds_array[i] = (ray_samples['s'][valid_idx[good_idx[1]]] -
                              ray_samples['s'][valid_idx[good_idx[0]]])
            else:
                ds_array[i] = 0.1
        elif i == n_good - 1:
            ds_array[i] = (ray_samples['s'][valid_idx[good_idx[i]]] -
                          ray_samples['s'][valid_idx[good_idx[i-1]]])
        else:
            ds_array[i] = ((ray_samples['s'][valid_idx[good_idx[i+1]]] -
                           ray_samples['s'][valid_idx[good_idx[i-1]]]) / 2.0)
    
    ds_cm = ds_array * RJ_CM
    
    # Combine all wavelengths and intensities
    xwavi_all = np.concatenate([xwavi_struct['sp'], xwavi_struct['s2p'],
                                xwavi_struct['s3p'], xwavi_struct['s4p'],
                                xwavi_struct['op'], xwavi_struct['o2p']])
    n_lines = len(xwavi_all)
    yptsi_all = np.zeros(n_lines)
    
    # Integrate along LOS
    for i in range(n_good):
        idx = good_idx[i]
        
        if nel_array[idx] < 1.0 or Tel_array[idx] < 0.1:
            continue
        
        # Get emissivities
        emiss_sp = interpolate_emissivity_2D(Tel_array[idx], nel_array[idx],
                                            temp_arr, dens_arr, yptsi_struct['sp'])
        emiss_s2p = interpolate_emissivity_2D(Tel_array[idx], nel_array[idx],
                                             temp_arr, dens_arr, yptsi_struct['s2p'])
        emiss_s3p = interpolate_emissivity_2D(Tel_array[idx], nel_array[idx],
                                             temp_arr, dens_arr, yptsi_struct['s3p'])
        emiss_s4p = interpolate_emissivity_2D(Tel_array[idx], nel_array[idx],
                                             temp_arr, dens_arr, yptsi_struct['s4p'])
        emiss_op = interpolate_emissivity_2D(Tel_array[idx], nel_array[idx],
                                            temp_arr, dens_arr, yptsi_struct['op'])
        emiss_o2p = interpolate_emissivity_2D(Tel_array[idx], nel_array[idx],
                                             temp_arr, dens_arr, yptsi_struct['o2p'])
        
        # Accumulate emission
        yptsi_step = np.concatenate([
            emiss_sp * nsp_array[idx],
            emiss_s2p * ns2p_array[idx],
            emiss_s3p * ns3p_array[idx],
            emiss_s4p * ns4p_array[idx],
            emiss_op * nop_array[idx],
            emiss_o2p * no2p_array[idx]
        ])
        
        yptsi_all += 1e-6 * yptsi_step * ds_cm[i]
    
    # Sort by wavelength
    wsort = np.argsort(xwavi_all)
    xwavi_sorted = xwavi_all[wsort]
    yptsi_sorted = yptsi_all[wsort]
    
    # Apply instrument response
    spectrum = simulate_IPT_spectrum_Rayleighs_ERF_form(xwav, spec_binsize,
                                                       xwavi_sorted, yptsi_sorted, fwhm)
    
    return spectrum


# =============================================================================
# MAIN FUNCTION
# =============================================================================

def main():
    """
    Main procedure for calculating IPT UV emission using field-aligned
    diffusive equilibrium model.
    """
    
    # =========================================================================
    # SETUP
    # =========================================================================
    print('================================================================')
    print('CITEP - Community Io Torus Emission Package')
    print('Field-Aligned Diffusive Equilibrium Model Version')
    print('================================================================')
    
    # Wavelength grid
    xmin = 550.0
    xmax = 2100.0
    binsize = 1.0
    nx = int(np.round((xmax - xmin)/binsize)) + 1
    xwav = np.arange(nx) * binsize + xmin
    
    # Instrument resolution
    fwhm = 4.47  # [Å]
    
    # =========================================================================
    # LOAD MODELS
    # =========================================================================
    print('\nLoading 3D field-aligned torus model...')
    field_model = load_3d_torus_model()
    
    print('\nLoading CHIANTI emission tables...')
    # Load from numpy compressed file
    chianti_data = np.load('CHIANTI_11.0.2_emiss_arrays_all_species_all_wavelengths_50x50_logspaced.npz')
    
    temp_arr = chianti_data['temp_arr']
    dens_arr = chianti_data['dens_arr']
    
    # Load wavelengths for each species
    xwavi = {
        'sp': chianti_data['xwavi_sp'],
        's2p': chianti_data['xwavi_s2p'],
        's3p': chianti_data['xwavi_s3p'],
        's4p': chianti_data['xwavi_s4p'],
        'op': chianti_data['xwavi_op'],
        'o2p': chianti_data['xwavi_o2p']
    }
    
    # Load emissivity tables for each species
    yptsi = {
        'sp': chianti_data['yptsi_sp'],
        's2p': chianti_data['yptsi_s2p'],
        's3p': chianti_data['yptsi_s3p'],
        's4p': chianti_data['yptsi_s4p'],
        'op': chianti_data['yptsi_op'],
        'o2p': chianti_data['yptsi_o2p']
    }
    
    print('Tables loaded:')
    print(f'  Temperature range: {temp_arr.min():.1f} - {temp_arr.max():.1f} eV')
    print(f'  Density range: {dens_arr.min():.1e} - {dens_arr.max():.1e} cm^-3')
    
    # =========================================================================
    # DISPLAY MODEL STRUCTURE
    # =========================================================================
    print('\n================================================================')
    print('FIELD-ALIGNED MODEL STRUCTURE')
    print('================================================================')
    print('Grid dimensions: 501 x 360 x 1001')
    print('  501 field lines: ρ_eq = 5.00 to 10.00 R_J (0.01 R_J steps)')
    print('  360 azimuths: φ = 0° to 359° (1° steps)')
    print('  1001 latitudes: λ_III = -50° to +50° (0.1° steps)')
    print('')
    print('Magnetic field model: JRM33 + CON2020')
    print('Plasma transport: Diffusive equilibrium along field lines')
    
    # Show sample densities at equator
    print('\nSample equatorial densities at φ = 180°:')
    print('  ρ [R_J]   n_e [cm^-3]   n_S++ [cm^-3]   n_O+ [cm^-3]')
    for i in range(0, 501, 50):
        rho_val = 5.0 + i * 0.01
        ne_eq = field_model['nel'][i, 180, 500]  # Latitude index 500 = equator
        ns2p_eq = field_model['ns2p'][i, 180, 500]
        nop_eq = field_model['nop'][i, 180, 500]
        if ne_eq > 0:
            print(f'{rho_val:7.2f} {ne_eq:14.3e} {ns2p_eq:14.3e} {nop_eq:14.3e}')
    
    # =========================================================================
    # CALCULATE SPECTRA
    # =========================================================================
    print('\n================================================================')
    print('CALCULATING UV EMISSION SPECTRA')
    print('================================================================')
    
    ds = 0.05  # Integration step size [R_J] - finer for better accuracy
    
    # Example 1: Equatorial ansa observation
    print('\nExample 1: Equatorial ansa observation')
    ray_origin = np.array([6.0, -20.0, 0.0])
    ray_dir = np.array([0.0, 1.0, 0.0])
    
    ray_samples = sample_ray_through_torus(ray_origin, ray_dir, ds)
    print(f'  Ray samples: {ray_samples["n_valid"]} points through torus')
    
    spectrum_eq = integrate_emission_along_los_field_model(
        ray_samples, field_model, xwav, binsize, fwhm,
        temp_arr, dens_arr, xwavi, yptsi)
    
    # Example 2: Off-equatorial observation
    print('\nExample 2: Off-equatorial observation (30° latitude)')
    ray_origin = np.array([6.0, -20.0, 3.46])  # 30° latitude at 6 R_J
    ray_dir = np.array([0.0, 1.0, 0.0])
    
    ray_samples = sample_ray_through_torus(ray_origin, ray_dir, ds)
    print(f'  Ray samples: {ray_samples["n_valid"]} points through torus')
    
    spectrum_off = integrate_emission_along_los_field_model(
        ray_samples, field_model, xwav, binsize, fwhm,
        temp_arr, dens_arr, xwavi, yptsi)
    
    # Example 3: Dawn vs Dusk asymmetry
    print('\nExample 3: Dawn observation (φ = 90°)')
    ray_origin = np.array([-20.0, 6.0, 0.0])  # Dawn side
    ray_dir = np.array([1.0, 0.0, 0.0])
    
    ray_samples = sample_ray_through_torus(ray_origin, ray_dir, ds)
    spectrum_dawn = integrate_emission_along_los_field_model(
        ray_samples, field_model, xwav, binsize, fwhm,
        temp_arr, dens_arr, xwavi, yptsi)
    
    print('Example 4: Dusk observation (φ = 270°)')
    ray_origin = np.array([20.0, -6.0, 0.0])  # Dusk side
    ray_dir = np.array([-1.0, 0.0, 0.0])
    
    ray_samples = sample_ray_through_torus(ray_origin, ray_dir, ds)
    spectrum_dusk = integrate_emission_along_los_field_model(
        ray_samples, field_model, xwav, binsize, fwhm,
        temp_arr, dens_arr, xwavi, yptsi)
    
    # =========================================================================
    # PLOTS
    # =========================================================================
    
    # Equatorial vs off-equatorial
    fig1, ax1 = plt.subplots(figsize=(10, 6))
    ax1.plot(xwav, spectrum_eq, 'k-', linewidth=2, label='Equatorial')
    ax1.plot(xwav, spectrum_off, 'b-', linewidth=2, alpha=0.7, label='30° Latitude')
    ax1.set_xlabel('Wavelength (Å)', fontsize=12)
    ax1.set_ylabel('Intensity (R/Å)', fontsize=12)
    ax1.set_title('Field-Aligned Model: Latitude Effects', fontsize=14)
    ax1.set_xlim(550, 2100)
    ax1.legend(loc='upper right')
    ax1.grid(True, alpha=0.3)
    plt.tight_layout()
    
    # Dawn-Dusk asymmetry
    fig2, ax2 = plt.subplots(figsize=(10, 6))
    ax2.plot(xwav, spectrum_dawn, 'r-', linewidth=2, label='Dawn (90°)')
    ax2.plot(xwav, spectrum_dusk, 'orange', linewidth=2, alpha=0.7, label='Dusk (270°)')
    ax2.set_xlabel('Wavelength (Å)', fontsize=12)
    ax2.set_ylabel('Intensity (R/Å)', fontsize=12)
    ax2.set_title('Field-Aligned Model: Dawn-Dusk Asymmetry', fontsize=14)
    ax2.set_xlim(550, 2100)
    ax2.legend(loc='upper right')
    ax2.grid(True, alpha=0.3)
    plt.tight_layout()
    
    # =========================================================================
    # 2D EMISSION MAP
    # =========================================================================
    print('\n================================================================')
    print('GENERATING 2D EMISSION MAP')
    print('================================================================')
    
    n_map = 31  # Reduced for field model due to computational cost
    x_positions = np.linspace(-10, 10, n_map)
    z_positions = np.linspace(-10, 10, n_map)
    
    total_emission = np.zeros((n_map, n_map))
    
    ray_dir = np.array([0.0, 1.0, 0.0])
    
    print(f'Calculating {n_map}x{n_map} pixel map...')
    
    for i in range(n_map):
        for j in range(n_map):
            ray_origin = np.array([x_positions[i], -20.0, z_positions[j]])
            
            ray_samples = sample_ray_through_torus(ray_origin, ray_dir, ds)
            
            if ray_samples['n_valid'] > 0:
                spectrum = integrate_emission_along_los_field_model(
                    ray_samples, field_model, xwav, binsize, fwhm,
                    temp_arr, dens_arr, xwavi, yptsi)
                total_emission[j, i] = np.sum(spectrum)
        
        if i % 5 == 0:
            print(f'  Row {i+1} of {n_map} complete')
    
    # Create contour plot
    fig3, ax3 = plt.subplots(figsize=(10, 8))
    im = ax3.contourf(x_positions, z_positions, total_emission, levels=20, cmap='viridis')
    ax3.set_xlabel('X Position (R_J)', fontsize=12)
    ax3.set_ylabel('Z Position (R_J)', fontsize=12)
    ax3.set_title('Total UV Emission: Field-Aligned Model', fontsize=14)
    ax3.set_aspect('equal')
    cbar = plt.colorbar(im, ax=ax3, label='Integrated Intensity (R)')
    plt.tight_layout()
    
    # =========================================================================
    # SUMMARY
    # =========================================================================
    print('\n================================================================')
    print('CALCULATION COMPLETE')
    print('================================================================')
    print('Model features captured:')
    print('  - Non-dipolar magnetic field (JRM33+CON2020)')
    print('  - Field-aligned plasma distribution')
    print('  - Dawn-dusk asymmetries')
    print('  - Realistic vertical stratification')
    print('')
    print('Total integrated emission (550-2100 Å):')
    print(f'  Equatorial: {np.sum(spectrum_eq):.1f} Rayleighs')
    print(f'  30° Latitude: {np.sum(spectrum_off):.1f} Rayleighs')
    print(f'  Dawn: {np.sum(spectrum_dawn):.1f} Rayleighs')
    print(f'  Dusk: {np.sum(spectrum_dusk):.1f} Rayleighs')
    print('')
    print(f'Dawn/Dusk ratio: {np.sum(spectrum_dawn)/np.sum(spectrum_dusk):.3f}')
    
    # Show all plots
    plt.show()


if __name__ == '__main__':
    main()
