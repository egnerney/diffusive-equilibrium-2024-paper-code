# -*- coding: utf-8 -*-
"""
Complete code to compute and plot parallel and perpendicular temperatures for all species
and cases, including the requested adjustments for Cases C-F, and plotting.

@author: Assistant
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker  # Import ticker module
from scipy.optimize import bisect
from scipy.interpolate import interp1d
import sys
import mpmath as mp_mp
from scipy.special import hyp2f1, expn

# Set mpmath precision
mp_mp.mp.dps = 15  # 15 decimal places

# Define the Species class
class Species:
    def __init__(self, name, m, q, T, A, kappa, lam, n, dist_type):
        self.name = name
        self.m = m
        self.q = q
        self.T = T
        self.A = A
        self.lam = lam
        self.kappa = kappa
        self.n = n
        self.type = dist_type

# Function to define planet-specific variables
def define_planet(planet):
    if planet.lower() == 'jupiter':
        RP = 7.1492e7       # Radius of Jupiter in meters
        Omega = 1.7583586898863878e-4  # Rotation rate in rad/s
        GM = 1.2668689059e17  # Gravitational constant * mass (m^3/s^2)
        return RP, Omega, GM
    elif planet.lower() == 'saturn':
        RP = 6.0268e7       # Radius of Saturn in meters
        Omega = 1.63788e-4  # Rotation rate in rad/s
        GM = 3.7931187e16   # Gravitational constant * mass (m^3/s^2)
        return RP, Omega, GM
    else:
        print('Planet ' + planet + ' is not supported')
        sys.exit(1)

# Isotropic Maxwellian density function
def maxwellian_density(species, deltaU, Phi, Bratio):
    n = species.n * np.exp(-(species.m * deltaU + species.q * Phi) / species.T)
    return n

# Anisotropic Maxwellian density function
def aniso_maxwellian_density(species, deltaU, Phi, Bratio):
    B0_over_B = 1.0 / Bratio  # B0 / B

    denom = species.A + (1.0 - species.A) * B0_over_B
    if denom <= 0:
        print("Warning: Non-Physical Density (Denominator negative or zero), setting denom to a small positive number.")
        denom = np.finfo(float).eps  # Avoid division by zero

    n = (species.n / denom) * np.exp(-(species.m * deltaU + species.q * Phi) / species.T)
    return n

# Anisotropic Kappa density function, Correlated Case
def aniso_kappa_density(species, deltaU, Phi, Bratio):
    B0_over_B = 1.0 / Bratio  # B0 / B

    denom = species.A + (1.0 - species.A) * B0_over_B
    if denom <= 0.:
        print("Warning: Non-Physical Density (Denominator negative or zero), setting denom to small positive number.")
        denom = np.finfo(float).eps  # Avoid division by zero

    PEratio = (species.m * deltaU + species.q * Phi) / (species.kappa * species.T)
    if PEratio <= -1.:
        PEratio  = -1. + np.finfo(float).eps  # setting to -1 + epsilon to prevent 0 or negative density while bisecting to find phi.
    n = (species.n / denom) * (1. + PEratio ) ** (0.5 - species.kappa)
    return n

# Anisotropic Product Kappa density function, unCorrelated Case
def aniso_product_kappa_density(species, deltaU, Phi, Bratio):
    if Bratio < 1.:
        print('integral does not converge if B/B0 < 1 setting to 1 + epsilon to proceed')
        Bratio = 1. + np.finfo(float).eps
    B0_over_B = 1.0 / Bratio  # B0 / B

    denom = species.A + (1.0 - species.A) * B0_over_B
    if denom <= 0.:
        print("Warning: Non-Physical Density (Denominator negative or zero), setting denom to small positive number.")
        denom = np.finfo(float).eps  # Avoid division by zero

    PEratio = (species.m * deltaU + species.q * Phi) / (species.kappa * species.T)
    if PEratio <= -1.:
        PEratio  = -1. + np.finfo(float).eps  # setting to -1 + epsilon to prevent 0 or negative density while bisecting to find phi.

    kappa_par = species.kappa
    kappa_perp = species.kappa * species.lam
    if ((kappa_par < 17.) and (kappa_perp < 17.)):
        delta = (species.lam* species.A)*((Bratio -1.) / (1. + PEratio))

        result  = kappa_perp * hyp2f1(1., kappa_par + 0.5,  kappa_perp + kappa_par + 1.5, 1. - delta) / (kappa_par + kappa_perp + 0.5)
        n =(species.n *Bratio )*( (1. + PEratio ) ** (0.5 - species.kappa))*result
    else:
        n = (species.n / denom) * (1. + PEratio ) ** (0.5 - species.kappa)
    return n

# Fried-Egg distribution density function
def fried_egg_density(species, deltaU, Phi, Bratio):
    if Bratio < 1.0:
        print("Warning: Density integral doesn't converge if Bratio = B/B0 < 1, setting Bratio to 1.0")
        Bratio = 1.0

    x = species.A * (Bratio - 1.0)
    exponent = -(species.m * deltaU + species.q * Phi) / species.T  # - Delta PE / parallel temp

    if Bratio != 1.0:
        if species.kappa > 15.:
            # Use approximation 1 / (1 + z)
            n = (species.n / (species.A + (1. - species.A) * (1. / Bratio))) * np.exp(exponent)
        else:
            # Use exact computation using mpmath
            func_value = float(species.kappa * mp_mp.expint(1.0 + species.kappa, species.kappa * x) * mp_mp.exp(species.kappa * x))
            n = species.n * Bratio * func_value * np.exp(exponent)
    else:
        n = species.n * np.exp(exponent)  # Isotropic Maxwellian solution when Bratio = 1

    return n

# Dictionary to hold different density functions
density_functions = {
    'Maxwellian': maxwellian_density,
    'Aniso_Maxwellian': aniso_maxwellian_density,
    'Aniso_kappa': aniso_kappa_density,
    'Aniso_product_kappa': aniso_product_kappa_density,
    'Fried_Egg': fried_egg_density
}

# Function to compute densities
def diff_eq_eddie(Pos_in, Pos0_in, Bratio, species0_in, phi, posc, nchoke, phic, planet, fcor, if_psh):
    # Set planet-specific variables
    RP, Omega, GM = define_planet(planet)
    Omega *= fcor

    # Convert species parameters to SI units
    species0 = []
    for sp in species0_in:
        sp_new = Species(
            name=sp.name,
            m=sp.m * 1.67e-27,    # Convert mass from AMU to kg
            q=sp.q * 1.6e-19,     # Convert charge from elementary charges to Coulombs
            T=sp.T * 1.6e-19,     # Convert temperature from eV to Joules
            A=sp.A,
            kappa=sp.kappa,
            lam=sp.lam,
            n=sp.n,
            dist_type=sp.type
        )
        species0.append(sp_new)
    nspec = len(species0)

    # Ensure the plasma is charge neutral at Pos0
    # Adjust electron density
    # Assuming species0[nspec-1] is electrons and species0[nspec-2] is hot electrons
    total_ion_charge_density = np.sum([species0[i].q * species0[i].n for i in range(0, nspec - 2)])
    species0[nspec - 1].n = -total_ion_charge_density / species0[nspec - 1].q  # Adjust cold electron density

    # Calculate gravitational-centrifugal potentials at Pos0 and Pos_in
    Pos0_r = Pos0_in[0] * RP  # Convert radius to meters
    Pos_r = Pos_in[0] * RP    # Convert radius to meters

    U0 = -GM / Pos0_r - 0.5 * Pos0_r ** 2. * np.cos(np.deg2rad(Pos0_in[1])) ** 2. * Omega ** 2.
    U = -GM / Pos_r - 0.5 * Pos_r ** 2. * np.cos(np.deg2rad(Pos_in[1])) ** 2. * Omega ** 2.
    deltaU = U - U0

    def net_charge_density(Phi):
        n_local = np.zeros(nspec)
        for i in range(nspec):
            density_func = density_functions.get(species0[i].type, maxwellian_density)
            n_local[i] = density_func(species0[i], deltaU, Phi, Bratio)
        # Exclude hot electrons from charge neutrality calculation
        nq_temp = [species0[i].q * n_local[i] for i in range(0, nspec - 2)]
        nq_temp.append(species0[nspec - 1].q * n_local[nspec - 1])  # Include cold electrons
        nq = np.sum(nq_temp)
        return nq

    # Initialize variables for root-finding
    Phi = 0.0
    Phi2 = 0.0
    dPhi = 1.0

    # Compute initial net charge densities
    nq = net_charge_density(Phi)
    nq2 = net_charge_density(Phi2)

    # Adjust Phi and Phi2 to bracket the root
    max_iterations = 1000
    iteration = 0
    while nq * nq2 > 0 and iteration < max_iterations:
        Phi -= dPhi
        Phi2 += dPhi
        nq = net_charge_density(Phi)
        nq2 = net_charge_density(Phi2)
        iteration += 1
    if iteration >= max_iterations:
        print("Failed to bracket the root.")
        return None

    # Use bisection method to find Phi where net charge density is zero
    Phi_root = bisect(net_charge_density, Phi, Phi2, xtol=1e-5)

    # Compute densities at Phi_root
    n = np.zeros(nspec)
    for i in range(nspec):
        density_func = density_functions.get(species0[i].type, maxwellian_density)
        n[i] = density_func(species0[i], deltaU, Phi_root, Bratio)
    phi = Phi_root
    # Return densities, phi, and deltaU
    return n, phi, deltaU

# Function to process a single field line
def process_field_line(case_name, field_line_index, species0_in, planet, fcor, if_psh, x_fl, y_fl, z_fl, rho_fl, lat_fl, phi_fl_wlong, r_fl, Btot_fl, r_00,
         nop0, no2p0, nsp0, ns2p0, ns3p0, nhp0, nnap0, noph0, nec0, neh0,
         ti0, thp0, tec0, teh0, toph0, kappa_Te, kappa_Top, kappa_e, kappa_op, kappa_Thp, pkappa_Te, pkappa_Top, pkappa_e, pkappa_op, pkappa_Thp):
    j = field_line_index
    # Set planet-specific variables
    RP, Omega, GM = define_planet(planet)
    Omega *= fcor

    nspec = len(species0_in)

    # Create double precision arrays
    dlatfl = 0.1
    lat_int_fl_max = 70.0
    lat_int_fl_min = -lat_int_fl_max
    npointsfieldline = round((lat_int_fl_max - lat_int_fl_min) / dlatfl) + 1
    lat_int_fl = lat_int_fl_min + dlatfl * np.arange(npointsfieldline)

    npoints = npointsfieldline

    densities_keys = ['op', 'o2p', 'sp', 's2p', 's3p', 'hp', 'nap', 'oph', 'eh', 'elec']
    n_out_line = {key: np.zeros(npoints) for key in densities_keys}
    T_out_line = {key: np.zeros(npoints) for key in densities_keys}

    # Initialize temperature dictionaries
    T_para = {key: np.zeros(npoints) for key in densities_keys}
    T_perp = {key: np.zeros(npoints) for key in densities_keys}

    rho_out_line = np.zeros(npoints, dtype=np.float64)
    z_out_line = np.zeros(npoints, dtype=np.float64)
    x_out_line = np.zeros(npoints, dtype=np.float64)
    y_out_line = np.zeros(npoints, dtype=np.float64)
    B_out_line = np.zeros(npoints, dtype=np.float64)

    r_fl_ceq = None  # Equatorial radial distance

    # Start processing field line 'j'
    # Interpolate field line data
    idx_finite = np.where(np.isfinite(x_fl[:, j]))[0]

    lat_min_bound = lat_fl[idx_finite, j].min()
    lat_max_bound = lat_fl[idx_finite, j].max()

    x_int_fl = interp1d(lat_fl[idx_finite, j], x_fl[idx_finite, j], fill_value="extrapolate")(lat_int_fl)
    y_int_fl = interp1d(lat_fl[idx_finite, j], y_fl[idx_finite, j], fill_value="extrapolate")(lat_int_fl)
    z_int_fl = interp1d(lat_fl[idx_finite, j], z_fl[idx_finite, j], fill_value="extrapolate")(lat_int_fl)
    rho_int_fl = interp1d(lat_fl[idx_finite, j], rho_fl[idx_finite, j], fill_value="extrapolate")(lat_int_fl)
    B_int_fl = interp1d(lat_fl[idx_finite, j], Btot_fl[idx_finite, j], fill_value="extrapolate")(lat_int_fl)
    r_int_fl = interp1d(lat_fl[idx_finite, j], r_fl[idx_finite, j], fill_value="extrapolate")(lat_int_fl)
    phi_int_fl_wlong = interp1d(lat_fl[idx_finite, j], phi_fl_wlong[idx_finite, j], fill_value="extrapolate")(lat_int_fl)
    # Create positions
    Pos = np.array([r_int_fl, lat_int_fl, phi_int_fl_wlong, B_int_fl]).T

    idx_ceq = np.argmax(rho_int_fl)
    r_fl_ceq = r_int_fl[idx_ceq]
    lat_fl_ceq = lat_int_fl[idx_ceq]
    longw_fl_ceq = phi_int_fl_wlong[idx_ceq]
    B_fl_ceq = B_int_fl[idx_ceq]

    if (r_fl_ceq >= 4.00) and (r_fl_ceq <= 10.0):
        Pos0 = [r_fl_ceq, lat_fl_ceq, longw_fl_ceq, B_fl_ceq]

        # Now, set species parameters according to the case
        species_list = setup_species_for_case(case_name, species0_in, r_fl_ceq, r_00, nop0, no2p0, nsp0, ns2p0, ns3p0,
                                              nhp0, nnap0, noph0, nec0, neh0, kappa_e, kappa_op, kappa_Te, kappa_Top, kappa_Thp, pkappa_e, pkappa_op, pkappa_Te, pkappa_Top, pkappa_Thp, ti0, thp0, teh0, toph0)

        # Prepare arrays for computation
        n_ions = np.zeros((nspec, npoints))
        phis = np.zeros(npoints)
        deltaUs = np.zeros(npoints)

        # Magnetic field ratio
        Bratio = B_int_fl / B_fl_ceq
        B0_over_B = 1.0 / Bratio
        phi = None
        posc = None
        nchoke = None
        phic = None

        # Compute deltaU along the field line
        deltaUs_array = np.zeros(npoints)
        for i in range(npoints):
            # Calculate deltaU at each point
            Pos_r = Pos[i, 0] * RP
            U = -GM / Pos_r - 0.5 * Pos_r ** 2. * np.cos(np.deg2rad(Pos[i, 1])) ** 2. * Omega ** 2.
            deltaUs_array[i] = U - (-GM / (Pos0[0]*RP) - 0.5 * (Pos0[0]*RP) ** 2. * np.cos(np.deg2rad(Pos0[1])) ** 2. * Omega ** 2.)

        # Loop over each point to compute densities along field line
        for i in range(npoints):
            if ((Pos[i, 0] > 1.0) and (lat_min_bound <= lat_int_fl[i] <= lat_max_bound) and
                    (rho_int_fl[i] > 0.) and (np.abs(z_int_fl[i]) <= 5.) and (rho_int_fl[i] <= 10.)):
                result = diff_eq_eddie(
                    Pos[i, 0:3],
                    Pos0[0:3],
                    Bratio[i],
                    species_list,
                    phi,
                    posc,
                    nchoke,
                    phic,
                    planet,
                    fcor,
                    if_psh
                )
                if result is not None:
                    n_ions[:, i], phis[i], deltaUs[i] = result
                else:
                    n_ions[:, i] = np.zeros(nspec)
                    phis[i] = np.nan
                    deltaUs[i] = np.nan
            else:
                n_ions[:, i] = np.zeros(nspec)
                phis[i] = np.nan
                deltaUs[i] = np.nan

        # Assign computed densities to n_out_line
        n_out_line['op'][:] = n_ions[species_names.index('O+'), :]
        n_out_line['o2p'][:] = n_ions[species_names.index('O++'), :]
        n_out_line['sp'][:] = n_ions[species_names.index('S+'), :]
        n_out_line['s2p'][:] = n_ions[species_names.index('S2+'), :]
        n_out_line['s3p'][:] = n_ions[species_names.index('S3+'), :]
        n_out_line['hp'][:] = n_ions[species_names.index('H+'), :]
        n_out_line['nap'][:] = n_ions[species_names.index('Na+'), :]
        n_out_line['oph'][:] = n_ions[species_names.index('O+(hot)'), :]
        n_out_line['eh'][:] = n_ions[species_names.index('eh-'), :]
        n_out_line['elec'][:] = n_ions[species_names.index('e-'), :]

        T_out_line['elec'][:] = interp1d(r_00, tec0, fill_value="extrapolate")(r_fl_ceq)
        T_out_line['eh'][:] = interp1d(r_00, teh0, fill_value="extrapolate")(r_fl_ceq)
        T_out_line['sp'][:] = interp1d(r_00, ti0, fill_value="extrapolate")(r_fl_ceq)
        T_out_line['oph'][:] = interp1d(r_00, toph0, fill_value="extrapolate")(r_fl_ceq)
        T_out_line['hp'][:] = interp1d(r_00, thp0, fill_value="extrapolate")(r_fl_ceq)

        x_out_line[:] = x_int_fl
        y_out_line[:] = y_int_fl
        z_out_line[:] = z_int_fl
        B_out_line[:] = B_int_fl

        # Now compute temperatures for each species
        for idx_sp, sp in enumerate(species_list):
            key = species_name_to_key[sp.name]
            # Get necessary parameters
            A0 = sp.A  # A_alpha0
            kappa = sp.kappa
            T0 = sp.T  # Characteristic temperature at s0 (eV)
            B0 = B_fl_ceq  # Magnetic field at s0
            B = B_int_fl  # Magnetic field along the field line
            BBratio=np.zeros(len(B))
            for iii in range(len(B)):
                if B[iii]/B0 < 1. :
                    BBratio[iii] = 1. 
                else:
                    BBratio[iii] = B[iii]/B0
                
            
            delta_PE = sp.m * 1.67e-27 * deltaUs_array + sp.q * 1.6e-19 * phis  # Delta PE in Joules

            if (sp.name in ['O+(hot)', 'eh-']) and case_name in ['C', 'D', 'E', 'F']:
                # Set temperatures to zero for hot species in Cases C-F
                T_para[key][:] = 0.0
                T_perp[key][:] = 0.0
                continue  # Skip to next species

            if case_name == 'A':
                # Isotropic Maxwellian
                T_para[key][:] = T0  # Constant with s
                T_perp[key][:] = T0  # Constant with s
            elif case_name == 'B':
                # Anisotropic Maxwellian with A=2
                T_para[key][:] = T0  # Constant with s   
                denom = A0 + (1. - A0) * (1./BBratio)
                T_perp[key][:] = (A0*T0) / denom
            elif case_name == 'C':
                # Isotropic Kappa
                # Characteristic temperature
                delta_PE_eV = delta_PE / 1.6e-19  # Convert to eV
                T_char = T0 * (1. + delta_PE_eV / (kappa * T0))
                T_mean = T_char # / (1. - 3. / (2. * kappa))
                T_para[key][:] = T_mean
                T_perp[key][:] = T_mean
            elif case_name == 'D':
                # Anisotropic Kappa
                # Parallel characteristic temperature
                delta_PE_eV = delta_PE / 1.6e-19  # Convert to eV
                T_para_char = T0 * (1. + delta_PE_eV / (kappa * T0))
                T_para_mean = T_para_char #/ (1. - 3. / (2. * kappa))
                T_para[key][:] = T_para_mean

                # Perpendicular characteristic temperature
                denom = A0 + (1. - A0) * (1./BBratio)
                T_perp_char = ((A0*T0) / denom) * (1. + delta_PE_eV / (kappa * T0))
                T_perp_mean = T_perp_char# / (1. - 3. / (2. * kappa))
                T_perp[key][:] = T_perp_mean
            elif case_name == 'E':
                # Anisotropic Product Kappa
                # Implement the complex expressions using hyp2f1

                # Initialize constants
                kappa_par0 = kappa
                kappa_perp0 = kappa * sp.lam  # Assuming sp.lam = 1.0
                delta_PE_eV = delta_PE / 1.6e-19  # Convert to eV
                delta_PE_ratio = delta_PE_eV / (kappa_par0 * T0)
                #delta_PE_ratio = np.clip(delta_PE_ratio, -1 + 1e-10, None)  # Avoid division by zero or negative sqrt

                delta_alpha = sp.lam * A0 * ((BBratio - 1.) / (1. + delta_PE_ratio))

                # Compute F(nu) and H(q)
                def F(nu):
                    return hyp2f1(nu, kappa_par0 + 0.5, kappa_par0 + kappa_perp0 + 1.5, 1. - delta_alpha)

                def H(q):
                    return hyp2f1(1, kappa_par0 + (1. - q) / 2., kappa_par0 + kappa_perp0 + (1. - q) / 2. + 1., 1 - delta_alpha)

                # Constants
                C1 = -2. * kappa_par0 - 2. * kappa_perp0 + 3.
                C2 = 2. * kappa_par0 + 2. * kappa_perp0 - 1.
                #C3 = -4. * kappa_par0 - 4 * kappa_perp0 + 6.
                C3 = 4. * (kappa_par0 + kappa_perp0) ** 2. - 1.
                C4 = 2. * kappa_par0 - 1.
                C5 = 2. * kappa_par0 - 3.
                C6 = 2. * kappa_par0 + 2. * kappa_perp0 + 1.

                # Compute T_perp_c and T_para_c
                F1 = F(1.)
                F2 = F(2.)
                F3 = F(3.)
                H2 = H(2.)
                H4 = H(4.)
                
                Tperp0 = A0*T0
                Tpar0 = T0

                numerator_perp = 2.*(BBratio) * kappa_perp0 * Tperp0 * F2 * F3
                denominator_perp = C1 *( F2 ** 2.) + 2.* C2 * F1 * F3
                T_perp_char = numerator_perp / denominator_perp

                numerator_para = 4. * C3 * kappa_par0 * Tpar0 * (1. + delta_PE_ratio) * H2 * H4
                denominator_para = 3. * C4 * (C2 ** 2.) * F1 * H4 + C5 * C1 * C6 * (H2 ** 2.)
                T_para_char = numerator_para / denominator_para

                # Mean temperatures
                T_para_mean = T_para_char# / (1. - 1. / (2. * kappa_par0))
                T_perp_mean = T_perp_char #/ (1. - 1. / kappa_perp0)

                T_para[key][:] = T_para_mean
                T_perp[key][:] = T_perp_mean

            elif case_name == 'F':
                # Fried-Egg distribution
                # Parallel temperature remains constant
                T_para[key][:] = T0

                # Perpendicular temperature using mpmath's arbitrary precision
                # Set higher precision for mpmath calculations
                mp_mp.mp.dps = 25  # You can adjust this value as needed

                # Convert scalar parameters to mpmath's mpf type
                A0_mpf = mp_mp.mpf(A0)
                kappa0_mpf = mp_mp.mpf(kappa)
                T0_mpf = mp_mp.mpf(A0*T0) #Tperp0

                # Convert arrays to lists and then to mpmath's mpf type
                Bratio_list = BBratio.tolist()
                Bratio_mpf = [mp_mp.mpf(bi) for bi in Bratio_list]

                # Compute x = A0 * (Bratio - 1)
                x_mpf = [A0_mpf * (bi - mp_mp.mpf(1.0)) for bi in Bratio_mpf]

                # Compute z = kappa0 * x
                z_mpf = [kappa0_mpf * xi for xi in x_mpf]

                # Compute E_kappa and E_kappa_minus1 using mpmath's expint
                E_kappa_mpf = [mp_mp.expint(kappa0_mpf, zi) for zi in z_mpf]
                E_kappa_minus1_mpf = [mp_mp.expint(kappa0_mpf - mp_mp.mpf(1.0), zi) for zi in z_mpf]

                numerator_mpf = []
                denominator_mpf = []

                # Loop over each element to compute numerator and denominator with high precision
                for i in range(len(x_mpf)):
                    xi = x_mpf[i]
                    Bi_over_B0 = Bratio_mpf[i]
                    #zi = z_mpf[i]
                    E_kappa_i = E_kappa_mpf[i]
                    E_kappa_minus1_i = E_kappa_minus1_mpf[i]

                    exp_kappa_x = mp_mp.exp(kappa0_mpf * xi)
                    x_plus_1 = xi + mp_mp.mpf(1.0)

                    # Compute terms step by step
                    term1 = Bi_over_B0 * kappa0_mpf * T0_mpf * xi

                    term2 = (
                        kappa0_mpf * exp_kappa_x * (kappa0_mpf * (x_plus_1) ** 2. - mp_mp.mpf(1.0)) * E_kappa_minus1_i
                        - (kappa0_mpf * x_plus_1 + mp_mp.mpf(1.0))
                    )

                    term3 = mp_mp.mpf(1.0) - kappa0_mpf * x_plus_1 * exp_kappa_x * E_kappa_i

                    numerator_i = term1 * term2 * term3

                    denom_term1 = mp_mp.mpf(2.0) * (kappa0_mpf - mp_mp.mpf(1.0))

                    denom_term2 = (
                        exp_kappa_x * E_kappa_i * (kappa0_mpf * xi * exp_kappa_x * E_kappa_i + kappa0_mpf * xi + kappa0_mpf - mp_mp.mpf(1.0))
                        - mp_mp.mpf(1.0)
                    )

                    denominator_i = denom_term1 * denom_term2

                    numerator_mpf.append(numerator_i)
                    denominator_mpf.append(denominator_i)

                # Compute T_perp for each element, converting from Joules to eV
                T_perp_list = []
                for i in range(len(numerator_mpf)):
                    # Avoid division by zero
                    if denominator_mpf[i] == 0:
                        T_perp_i = mp_mp.mpf('nan')  # Assign NaN if denominator is zero
                    else:
                        T_perp_i = (numerator_mpf[i] / denominator_mpf[i])#/(mp_mp.mpf(1.) - mp_mp.mpf(1.)/(mp_mp.mpf(1.)*kappa0_mpf))
                    T_perp_list.append(T_perp_i)

                # Convert the list of mpmath values to a NumPy array of floats
                T_perp[key][:] = np.array([float(tp) for tp in T_perp_list])
            else:
                print(f"Unknown case {case_name}")
                sys.exit(1)
    else:
        # Assign zeros if outside range
        for key in densities_keys:
            n_out_line[key][:] = 0.
            T_out_line[key][:] = 0.
            T_para[key][:] = 0.
            T_perp[key][:] = 0.

        x_out_line[:] = x_int_fl
        y_out_line[:] = y_int_fl
        z_out_line[:] = z_int_fl
        B_out_line[:] = B_int_fl
        phis[:] = np.nan
        deltaUs[:] = np.nan

    # Return the results including phi and deltaU and temperatures
    return (j, n_out_line, T_out_line, x_out_line, y_out_line, z_out_line, B_out_line, phis, deltaUs, species_list, T_para, T_perp)

# Function to set up species parameters for each case
def setup_species_for_case(case_name, species0_in, r_fl_ceq, r_00, nop0, no2p0, nsp0, ns2p0, ns3p0,
                           nhp0, nnap0, noph0, nec0, neh0, kappa_e, kappa_op, kappa_Te, kappa_Top, kappa_Thp, pkappa_e, pkappa_op, pkappa_Te, pkappa_Top, pkappa_Thp,  ti0, thp0, teh0, toph0):
    # Create a deep copy of species0_in
    import copy
    species_list = copy.deepcopy(species0_in)
    nspec = len(species_list)

    # Interpolate densities and kappa values at equator
    nopnew = interp1d(r_00, nop0, fill_value="extrapolate")(r_fl_ceq).item()
    no2pnew = interp1d(r_00, no2p0, fill_value="extrapolate")(r_fl_ceq).item()
    nspnew = interp1d(r_00, nsp0, fill_value="extrapolate")(r_fl_ceq).item()
    ns2pnew = interp1d(r_00, ns2p0, fill_value="extrapolate")(r_fl_ceq).item()
    ns3pnew = interp1d(r_00, ns3p0, fill_value="extrapolate")(r_fl_ceq).item()
    nhpnew = interp1d(r_00, nhp0, fill_value="extrapolate")(r_fl_ceq).item()
    nnapnew = interp1d(r_00, nnap0, fill_value="extrapolate")(r_fl_ceq).item()
    nophnew = interp1d(r_00, noph0, fill_value="extrapolate")(r_fl_ceq).item()
    necnew = interp1d(r_00, nec0, fill_value="extrapolate")(r_fl_ceq).item()
    nehnew = interp1d(r_00, neh0, fill_value="extrapolate")(r_fl_ceq).item()
    netotnew = necnew + nehnew

    kappa_Te_new = interp1d(r_00, kappa_Te, fill_value="extrapolate")(r_fl_ceq).item()
    kappa_Top_new = interp1d(r_00, kappa_Top, fill_value="extrapolate")(r_fl_ceq).item()
    kappa_e_new = interp1d(r_00, kappa_e, fill_value="extrapolate")(r_fl_ceq).item()
    kappa_op_new = interp1d(r_00, kappa_op, fill_value="extrapolate")(r_fl_ceq).item()
    kappa_Thp_new = interp1d(r_00, kappa_Thp, fill_value="extrapolate")(r_fl_ceq).item()
    
    pkappa_Te_new = interp1d(r_00, pkappa_Te, fill_value="extrapolate")(r_fl_ceq).item()
    pkappa_Top_new = interp1d(r_00, pkappa_Top, fill_value="extrapolate")(r_fl_ceq).item()
    pkappa_e_new = interp1d(r_00, pkappa_e, fill_value="extrapolate")(r_fl_ceq).item()
    pkappa_op_new = interp1d(r_00, pkappa_op, fill_value="extrapolate")(r_fl_ceq).item()
    pkappa_Thp_new = interp1d(r_00, pkappa_Thp, fill_value="extrapolate")(r_fl_ceq).item()
    

    #Ticnew = interp1d(r_00, kappa_Top, fill_value="extrapolate")(r_fl_ceq)
    Ticnew = interp1d(r_00, ti0, fill_value="extrapolate")(r_fl_ceq).item()
    #Thpnew = interp1d(r_00, kappa_Thp, fill_value="extrapolate")(r_fl_ceq)
    Thpnew = interp1d(r_00, thp0, fill_value="extrapolate")(r_fl_ceq).item()
    #Tecnew = interp1d(r_00, kappa_Te, fill_value="extrapolate")(r_fl_ceq)
    Tecnew = interp1d(r_00, tec0, fill_value="extrapolate")(r_fl_ceq).item()
    Tehnew = interp1d(r_00, teh0, fill_value="extrapolate")(r_fl_ceq).item()
    Tophnew = interp1d(r_00, toph0, fill_value="extrapolate")(r_fl_ceq).item()
    
    nop = nopnew
    ns2p = ns2pnew
    noph = nophnew

    if ns2p > 0.:
        ratio = nop / ns2p
        nophot = (noph * ratio) / (ratio + 1.)      # Hot O⁺ density
        ns2phot = noph / (ratio + 1.)               # Hot S²⁺ density
    else:
        nophot = noph  # Hot O⁺ density
        ns2phot = 0.0  # Hot S²⁺ density

    noptotnew = nop + nophot
    ns2ptotnew = ns2p + ns2phot

    # Update species parameters based on the case
    for sp in species_list:
        if sp.name == 'O+':
            sp.n = nopnew
            sp.T = Ticnew
            sp.kappa = kappa_op_new
            sp.lam = 1.0
        elif sp.name == 'O++':
            sp.n = no2pnew
            sp.T = Ticnew
            sp.kappa = kappa_op_new
            sp.lam = 1.0
        elif sp.name == 'S+':
            sp.n = nspnew
            sp.T = Ticnew
            sp.kappa = kappa_op_new
            sp.lam = 1.0
        elif sp.name == 'S2+':
            sp.n = ns2pnew
            sp.T = Ticnew
            sp.kappa = kappa_op_new
            sp.lam = 1.0
        elif sp.name == 'S3+':
            sp.n = ns3pnew
            sp.T = Ticnew
            sp.kappa = kappa_op_new
            sp.lam = 1.0
        elif sp.name == 'H+':
            sp.n = nhpnew
            sp.T = Thpnew
            sp.kappa = kappa_op_new
            sp.lam = 1.0
        elif sp.name == 'Na+':
            sp.n = nnapnew
            sp.T = Ticnew
            sp.kappa = kappa_op_new
            sp.lam = 1.0
        elif sp.name == 'O+(hot)':
            if case_name in ['C', 'D', 'E', 'F']:
                sp.n = 0.
                sp.T = 1.
            else:
                sp.n = nophnew
                sp.T = Tophnew  
            sp.kappa = kappa_op_new
            sp.lam = 1.0
        elif sp.name == 'eh-':
            if case_name in ['C', 'D', 'E', 'F']:
                sp.n = 0.
                sp.T = 1.
            else:
                sp.n = nehnew
                sp.T = Tehnew 
            sp.kappa = kappa_e_new
            sp.lam = 1.0
        elif sp.name == 'e-':
            sp.n = necnew
            sp.T = Tecnew
            sp.kappa = kappa_e_new
            sp.lam = 1.0

    # Now set the distribution types and anisotropy A based on the case
    if case_name == 'A':
        for sp in species_list:
            sp.type = 'Maxwellian'
            sp.A = 1.0
    elif case_name == 'B':
        for sp in species_list:
            sp.type = 'Aniso_Maxwellian'
            sp.A = 2.0
            sp.T = sp.T / 2.0  # T_parallel = T_isotropic / 2
    elif case_name == 'C':
        for sp in species_list:
            sp.type = 'Aniso_kappa'
            sp.A = 1.0
            if sp.name == 'O+':
                sp.n = noptotnew
                sp.T = kappa_Top_new
                sp.kappa = kappa_op_new
                sp.lam = 1.0
            elif sp.name == 'O++':
                sp.n = no2pnew
                sp.T = kappa_Top_new
                sp.kappa = kappa_op_new
                sp.lam = 1.0
            elif sp.name == 'S+':
                sp.n = nspnew
                sp.T = kappa_Top_new
                sp.kappa = kappa_op_new
                sp.lam = 1.0
            elif sp.name == 'S2+':
                sp.n = ns2ptotnew
                sp.T = kappa_Top_new
                sp.kappa = kappa_op_new
                sp.lam = 1.0
            elif sp.name == 'S3+':
                sp.n = ns3pnew
                sp.T = kappa_Top_new
                sp.kappa = kappa_op_new
                sp.lam = 1.0
            elif sp.name == 'H+':
                sp.n = nhpnew
                sp.T = kappa_Thp_new
                sp.kappa = kappa_op_new
                sp.lam = 1.0
            elif sp.name == 'Na+':
                sp.n = nnapnew
                sp.T = kappa_Top_new
                sp.kappa = kappa_op_new
                sp.lam = 1.0
            elif sp.name == 'O+(hot)':
                sp.n = 0.
                sp.T = 1.  
                sp.kappa = kappa_op_new
                sp.lam = 1.0
            elif sp.name == 'eh-':
                sp.n = 0.
                sp.T = 1. 
                sp.kappa = kappa_e_new
                sp.lam = 1.0
            elif sp.name == 'e-':
                sp.n = netotnew
                sp.T = kappa_Te_new
                sp.kappa = kappa_e_new
                sp.lam = 1.0
    elif case_name == 'D':
        for sp in species_list:
            sp.type = 'Aniso_kappa'
            sp.A = 2.0
            if sp.name == 'O+':
                sp.n = noptotnew
                sp.T = kappa_Top_new
                sp.kappa = kappa_op_new
                sp.lam = 1.0
            elif sp.name == 'O++':
                sp.n = no2pnew
                sp.T = kappa_Top_new
                sp.kappa = kappa_op_new
                sp.lam = 1.0
            elif sp.name == 'S+':
                sp.n = nspnew
                sp.T = kappa_Top_new
                sp.kappa = kappa_op_new
                sp.lam = 1.0
            elif sp.name == 'S2+':
                sp.n = ns2ptotnew
                sp.T = kappa_Top_new
                sp.kappa = kappa_op_new
                sp.lam = 1.0
            elif sp.name == 'S3+':
                sp.n = ns3pnew
                sp.T = kappa_Top_new
                sp.kappa = kappa_op_new
                sp.lam = 1.0
            elif sp.name == 'H+':
                sp.n = nhpnew
                sp.T =  kappa_Thp_new
                sp.kappa = kappa_op_new
                sp.lam = 1.0
            elif sp.name == 'Na+':
                sp.n = nnapnew
                sp.T = kappa_Top_new
                sp.kappa = kappa_op_new
                sp.lam = 1.0
            elif sp.name == 'O+(hot)':
                sp.n = 0.
                sp.T = 1.  
                sp.kappa = kappa_op_new
                sp.lam = 1.0
            elif sp.name == 'eh-':
                sp.n = 0.
                sp.T = 1. 
                sp.kappa = kappa_e_new
                sp.lam = 1.0
            elif sp.name == 'e-':
                sp.n = netotnew
                sp.T = kappa_Te_new
                sp.kappa = kappa_e_new
                sp.lam = 1.0
            sp.T = sp.T / 2.0  # T_parallel = T_from_C / 2
    elif case_name == 'E':
        for sp in species_list:
            sp.type = 'Aniso_product_kappa'
            sp.A = 1.0
            sp.lam = 1.0
            if sp.name == 'O+':
                sp.n = noptotnew
                sp.T = pkappa_Top_new
                sp.kappa = pkappa_op_new
            elif sp.name == 'O++':
                sp.n = no2pnew
                sp.T = pkappa_Top_new
                sp.kappa = pkappa_op_new
            elif sp.name == 'S+':
                sp.n = nspnew
                sp.T = pkappa_Top_new
                sp.kappa = pkappa_op_new
            elif sp.name == 'S2+':
                sp.n = ns2ptotnew
                sp.T = pkappa_Top_new
                sp.kappa = pkappa_op_new
            elif sp.name == 'S3+':
                sp.n = ns3pnew
                sp.T = pkappa_Top_new
                sp.kappa = pkappa_op_new
            elif sp.name == 'H+':
                sp.n = nhpnew
                sp.T = pkappa_Thp_new
                sp.kappa = pkappa_op_new
            elif sp.name == 'Na+':
                sp.n = nnapnew
                sp.T = pkappa_Top_new
                sp.kappa = pkappa_op_new
            elif sp.name == 'O+(hot)':
                sp.n = 0.
                sp.T = 1.
                sp.kappa = kappa_op_new
            elif sp.name == 'eh-':
                sp.n = 0.
                sp.T = 1.
                sp.kappa = kappa_e_new
            elif sp.name == 'e-':
                sp.n = netotnew
                sp.T = pkappa_Te_new
                sp.kappa = pkappa_e_new
    elif case_name == 'F':
        for sp in species_list:
            sp.type = 'Fried_Egg'
            if sp.name == 'O+':
                sp.n = noptotnew
                sp.T = Ticnew
                sp.A = pkappa_Top_new/Ticnew
                sp.kappa = pkappa_op_new
            elif sp.name == 'O++':
                sp.n = no2pnew
                sp.T = Ticnew
                sp.A = pkappa_Top_new/Ticnew
                sp.kappa = pkappa_op_new
            elif sp.name == 'S+':
                sp.n = nspnew
                sp.T = Ticnew
                sp.A = pkappa_Top_new/Ticnew
                sp.kappa = pkappa_op_new
            elif sp.name == 'S2+':
                sp.n = ns2ptotnew
                sp.T = Ticnew
                sp.A = pkappa_Top_new/Ticnew
                sp.kappa = pkappa_op_new
            elif sp.name == 'S3+':
                sp.n = ns3pnew
                sp.T = Ticnew
                sp.A = pkappa_Top_new/Ticnew
                sp.kappa = pkappa_op_new
            elif sp.name == 'H+':
                sp.n = nhpnew
                sp.T = Thpnew
                sp.A = pkappa_Thp_new/Thpnew
                sp.kappa = pkappa_op_new
            elif sp.name == 'Na+':
                sp.n = nnapnew
                sp.T = Ticnew
                sp.A = pkappa_Top_new/Ticnew
                sp.kappa = pkappa_op_new
            elif sp.name == 'O+(hot)':
                sp.n = 0.
                sp.T = 1.
                sp.kappa = kappa_op_new
            elif sp.name == 'eh-':
                sp.n = 0.
                sp.T = 1.
                sp.kappa = kappa_e_new
            elif sp.name == 'e-':
                sp.n = netotnew
                sp.T = Tecnew
                sp.A = pkappa_Te_new/Tecnew
                sp.kappa = pkappa_e_new
    else:
        print(f"Unknown case {case_name}")
        sys.exit(1)

    return species_list

if __name__ == '__main__':
    # Define species
    planet = 'Jupiter'
    fcor = 1.0  # azimuthal velocity relative to corotation
    if_psh = 0
    RP, Omega, GM = define_planet(planet)
    Omega *= fcor

    # Use your specified species_names
    species_names = ['O+', 'O++', 'S+', 'S2+', 'S3+', 'H+', 'Na+', 'O+(hot)', 'eh-', 'e-']
    species_keys = ['op', 'o2p', 'sp', 's2p', 's3p', 'hp', 'nap', 'oph', 'eh', 'elec']
    species_name_to_key = dict(zip(species_names, species_keys))
    nspec = len(species_names)
    species_m = [16.0, 16.0, 32.0, 64.0, 96.0, 1.0, 23.0, 16.0, 1.0 / 1837.0, 1.0 / 1837.0]  # AMU
    species_q = [1.0, 2.0, 1.0, 2.0, 3.0, 1.0, 1.0, 1.0, -1.0, -1.0]  # Elementary charges
    species_T = [79.3, 79.3, 79.3, 79.3, 79.3, 79.3, 94.1, 362.0, 46.0, 4.6]  # eV
    species_n = [592.0, 76.3, 163.0, 538.0, 90.7, 50.6, 97.2, 134.0, 2.5375, 2537.5]  # Density units
    species_A = [1.0] * nspec
    species_kappa = [100.0] * nspec
    species_lam = [1.0] * nspec
    species_type = ['Aniso_product_kappa'] * nspec

    # Create species list
    species0_in = []
    for i in range(nspec):
        sp = Species(
            name=species_names[i],
            m=species_m[i],
            q=species_q[i],
            T=species_T[i],
            A=species_A[i],
            kappa=species_kappa[i],
            lam=species_lam[i],
            n=species_n[i],
            dist_type=species_type[i]
        )
        species0_in.append(sp)

    # Load field line data (ensure these files exist)
    x_fl = np.transpose(np.loadtxt('int_aligned_x_jrm33+con20integral_601_4-10_properly_aligned_phi0=75.85_perfectaligned_for_4.70-9.71.csv', delimiter=','))
    y_fl = np.transpose(np.loadtxt('int_aligned_y_jrm33+con20integral_601_4-10_properly_aligned_phi0=75.85_perfectaligned_for_4.70-9.71.csv', delimiter=','))
    z_fl = np.transpose(np.loadtxt('int_aligned_z_jrm33+con20integral_601_4-10_properly_aligned_phi0=75.85_perfectaligned_for_4.70-9.71.csv', delimiter=','))
    rho_fl = np.sqrt(x_fl ** 2. + y_fl ** 2.)  # rho systemIII RJ
    r_fl = np.sqrt(x_fl ** 2. + y_fl ** 2. + z_fl ** 2.)  # spherical r systemIII RJ
    lat_fl = np.degrees(np.arcsin(z_fl / r_fl))  # sysIII lat degrees
    phi_fl = np.degrees(np.arctan2(y_fl, x_fl))  # elong degrees
    phi_fl_wlong = 360. - phi_fl

    r_00 = 0.01 * np.arange(601) + 4.0

    # Define field line indices for Io and Europa
    field_line_indices = [189, 538]  # 190 and 539 in 0-based indexing
    field_line_labels = {189: 'Io Field Line', 538: 'Europa Field Line'}

    # Load plasma parameters (ensure these files exist)
    nop0 = np.loadtxt('nop_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv', delimiter=',')
    no2p0 = np.loadtxt('no2p_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv', delimiter=',')
    nsp0 = np.loadtxt('nsp_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv', delimiter=',')
    ns2p0 = np.loadtxt('ns2p_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv', delimiter=',')
    ns3p0 = np.loadtxt('ns3p_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv', delimiter=',')
    nhp0 = np.loadtxt('nhp_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv', delimiter=',')
    nnap0 = np.loadtxt('nnap_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv', delimiter=',')
    noph0 = np.loadtxt('noph_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv', delimiter=',')
    nec0 = np.loadtxt('new_nominal_model_nec_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv', delimiter=',')
    neh0 = np.loadtxt('new_nominal_model_neh_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv', delimiter=',')
    # Temperatures
    ti0 = np.loadtxt('Tic_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv', delimiter=',')
    thp0 = np.loadtxt('Thp_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv', delimiter=',')
    tec0 = np.loadtxt('Tec_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv', delimiter=',')
    teh0 = np.loadtxt('new_nominal_model_Teh_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv', delimiter=',')
    toph0 = np.loadtxt('Toph_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv', delimiter=',')
    # Btot_fl
    Btot_fl = np.transpose(np.loadtxt('Btot_int_aligned_jrm33+con20integral_601_4-10_properly_aligned_phi0=75.85_perfectaligned_for_4.70-9.71.csv', delimiter=','))
    # Kappa parameters
    kappa_Te = np.loadtxt('new_nominal_model_Te_electrons_012_momentsequal_bimax_to_kappa.csv', delimiter=',')
    kappa_e = np.loadtxt('new_nominal_model_kappa_electrons_012_momentsequal_bimax_to_kappa.csv', delimiter=',')
    kappa_Top = np.loadtxt('Te_O_012_momentsequal_bimax_to_kappa.csv', delimiter=',')
    kappa_op = np.loadtxt('kappa_O_012_momentsequal_bimax_to_kappa.csv', delimiter=',')
    kappa_Thp = kappa_Top * (thp0 / ti0)
    
    pkappa_Te = kappa_Te #np.loadtxt('new_nominal_model_Te_electrons_012_momentsequal_bimax_to_product_kappa.csv', delimiter=',')
    pkappa_e = kappa_e #np.loadtxt('new_nominal_model_kappa_electrons_012_momentsequal_bimax_to_product_kappa.csv', delimiter=',')
    pkappa_Top = kappa_Top #np.loadtxt('Te_O_012_momentsequal_bimax_to_product_kappa.csv', delimiter=',')
    pkappa_op = kappa_op #np.loadtxt('kappa_O_012_momentsequal_bimax_to_product_kappa.csv', delimiter=',')
    pkappa_Thp = kappa_Thp #kappa_Top * (thp0 / ti0)

    # Define case names
    case_names = ['A', 'B', 'C', 'D', 'E', 'F']

    # Dictionary to store results
    results = {}

    for case_name in case_names:
        print(f"Processing case {case_name}")
        results[case_name] = {}
        for field_line_index in field_line_indices:
            # Process field line
            j, n_out_line, T_out_line, x_out_line, y_out_line, z_out_line, B_out_line, phis, deltaUs, species_list, T_para, T_perp = process_field_line(
                case_name, field_line_index, species0_in, planet, fcor, if_psh, x_fl, y_fl, z_fl, rho_fl, lat_fl,
                phi_fl_wlong, r_fl, Btot_fl, r_00,
                nop0, no2p0, nsp0, ns2p0, ns3p0, nhp0, nnap0, noph0, nec0, neh0,
                ti0, thp0, tec0, teh0, toph0, kappa_Te, kappa_Top, kappa_e, kappa_op, kappa_Thp, pkappa_Te, pkappa_Top, pkappa_e, pkappa_op, pkappa_Thp)

            # Compute total potential energy difference for each species
            total_PE = {}
            for i, sp in enumerate(species_list):
                total_PE[sp.name] = sp.m * 1.67e-27 * deltaUs + sp.q * 1.6e-19 * phis  # in Joules

            # Save results
            results[case_name][field_line_index] = {
                'phis': phis,  # in Volts
                'deltaUs': deltaUs,  # in Joules per kg
                'total_PE': total_PE,  # in Joules
                's': None,  # Will compute s below
                'species_list': species_list,
                'x': x_out_line,
                'y': y_out_line,
                'z': z_out_line,
                'n_out_line': n_out_line,
                'T_para': T_para,
                'T_perp': T_perp,
                'B_out_line': B_out_line
            }

            # Compute s along the field line
            x_line = x_out_line
            y_line = y_out_line
            z_line = z_out_line
            dx = np.diff(x_line)
            dy = np.diff(y_line)
            dz = np.diff(z_line)
            ds = np.sqrt(dx**2 + dy**2 + dz**2)
            s = np.concatenate(([0], np.cumsum(ds)))
            equator_index = np.argmin(np.abs(z_line))
            s = s - s[equator_index]

            results[case_name][field_line_index]['s'] = s

    # Now, generate the plots
    field_line_labels = {189: 'Io Field Line', 538: 'Europa Field Line'}

    # Define species labels
    species_keys = ['elec', 'eh', 'op', 'sp', 'o2p', 's2p', 'oph', 's3p', 'hp', 'nap']
    species_labels = ['e$^{-}$', 'e$^{-}$(hot)', 'O$^{+}$', 'S$^{+}$', 'O$^{++}$', 'S$^{2+}$', 'O$^{+}$(hot)', 'S$^{3+}$', 'H$^{+}$', 'Na$^{+}$']

    # Plotting temperatures for all species, in multiplot panels
    for field_line_index in field_line_indices:
        field_label = field_line_labels[field_line_index]
        fig, axs = plt.subplots(5, 2, figsize=(14, 20), sharex=True)
        fig.suptitle(f'Parallel Temperatures along {field_label} (Northern Hemisphere)', fontsize=20, fontweight='bold')
        for idx_sp, (key, label) in enumerate(zip(species_keys, species_labels)):
            row = idx_sp // 2
            col = idx_sp % 2
            ax = axs[row, col]
            for case_name in case_names:
                s = results[case_name][field_line_index]['s']
                T_para = results[case_name][field_line_index]['T_para'][key]
                # Only consider s >= 0 (northern hemisphere)
                positive_s_indices = s >= 0
                s_positive = s[positive_s_indices]
                T_para_positive = T_para[positive_s_indices]
                # For Europa, limit s ≤ 11.39
                if field_line_index == 538:
                    s_limit_indices = s_positive <= 11.39
                    s_positive = s_positive[s_limit_indices]
                    T_para_positive = T_para_positive[s_limit_indices]
                else:
                    s_limit_indices = s_positive <= 7.2
                    s_positive = s_positive[s_limit_indices]
                    T_para_positive = T_para_positive[s_limit_indices]
                    
                    
                # Exclude zero temperatures (for hot species in Cases C-F)
                if np.any(T_para_positive > 0):
                    ax.plot(s_positive, T_para_positive, label=f"Case {case_name}")
            #ax.set_ylabel('Temperature (eV)', fontsize=14, fontweight='bold')
            
            if (key == 'eh') or (key == 'oph'):
                ax.text(0.5, 0.75, f'{label}', transform=ax.transAxes, fontsize=22, fontweight='bold', ha='center')
            else:
                ax.text(0.5, 0.9, f'{label}', transform=ax.transAxes, fontsize=22, fontweight='bold', ha='center')
                
            
            ax.grid(True, which='both', linestyle='--')
            ax.tick_params(axis='both', which='both', labelsize=17, width=2)
            for tick in ax.get_xticklabels():
                tick.set_fontweight('bold')
            for tick in ax.get_yticklabels():
                tick.set_fontweight('bold')
            if row == 4:
                ax.set_xlabel(r's ($\bf{R_J}$)', fontsize=20, fontweight='bold')
        fig.text(0.04, 0.5, 'Parallel Temperature (eV)', va='center', rotation='vertical', fontsize=20, fontweight='bold')
        axs[3,1].legend(loc='upper left', frameon=False, prop={'weight': 'bold', 'size':16})
        plt.tight_layout(rect=[0.05, 0.03, 1, 0.99])
        plt.savefig(f'mean_Parallel_Temperatures_vs_s_Species_FieldLine_{field_label.replace(" ", "_")}_NorthernHemisphere.png', dpi=600)
        plt.show()

        # Similar plot for perpendicular temperatures
        fig, axs = plt.subplots(5, 2, figsize=(14, 20), sharex=True)
        fig.suptitle(f'Perpendicular Temperatures along {field_label} (Northern Hemisphere)', fontsize=20, fontweight='bold')
        for idx_sp, (key, label) in enumerate(zip(species_keys, species_labels)):
            row = idx_sp // 2
            col = idx_sp % 2
            ax = axs[row, col]
            for case_name in case_names:
                s = results[case_name][field_line_index]['s']
                T_perp = results[case_name][field_line_index]['T_perp'][key]
                # Only consider s >= 0 (northern hemisphere)
                positive_s_indices = s >= 0
                s_positive = s[positive_s_indices]
                T_perp_positive = T_perp[positive_s_indices]
                # For Europa, limit s ≤ 11.39
                if field_line_index == 538:
                    s_limit_indices = s_positive <= 11.39
                    s_positive = s_positive[s_limit_indices]
                    T_perp_positive = T_perp_positive[s_limit_indices]
                else:
                    s_limit_indices = s_positive <= 7.2
                    s_positive = s_positive[s_limit_indices]
                    T_perp_positive = T_perp_positive[s_limit_indices]
                # Exclude zero temperatures (for hot species in Cases C-F)
                if np.any(T_perp_positive > 0):
                    ax.plot(s_positive, T_perp_positive, label=f"Case {case_name}")
            #ax.set_ylabel('Temperature (eV)', fontsize=14, fontweight='bold')
            #ax.text(0.5, 0.9, f'{label}', transform=ax.transAxes, fontsize=22, fontweight='bold', ha='center')
            if (key == 'eh') or (key == 'oph'):
                ax.text(0.5, 0.75, f'{label}', transform=ax.transAxes, fontsize=22, fontweight='bold', ha='center')
            elif (key == 'hp'):
                ax.text(0.5, 0.8, f'{label}', transform=ax.transAxes, fontsize=22, fontweight='bold', ha='center')
            else:
                ax.text(0.5, 0.9, f'{label}', transform=ax.transAxes, fontsize=22, fontweight='bold', ha='center')
                
                
            ax.grid(True, which='both', linestyle='--')
            ax.tick_params(axis='both', which='both', labelsize=17, width=2)
            for tick in ax.get_xticklabels():
                tick.set_fontweight('bold')
            for tick in ax.get_yticklabels():
                tick.set_fontweight('bold')
            if row == 4:
                ax.set_xlabel(r's ($\bf{R_J}$)', fontsize=20, fontweight='bold')
        fig.text(0.04, 0.5, 'Perpendicular Temperature (eV)', va='center', rotation='vertical', fontsize=20, fontweight='bold')
        axs[3,1].legend(loc='upper left',frameon=False, prop={'weight': 'bold', 'size':16})
        plt.tight_layout(rect=[0.05, 0.03, 1, 0.99])
        plt.savefig(f'mean_Perpendicular_Temperatures_vs_s_Species_FieldLine_{field_label.replace(" ", "_")}_NorthernHemisphere.png', dpi=600)
        plt.show()
        
    # ... [Rest of your code above]

    # Now, generate the plots
    #field_line_labels = {189: 'Io Field Line', 538: 'Europa Field Line'}

    # Define species labels
    #species_keys = ['elec', 'eh', 'op', 'sp', 'o2p', 's2p', 'oph', 's3p']
    #species_labels = ['e$^{-}$', 'e$^{-}$(hot)', 'O$^{+}$', 'S$^{+}$', 'O$^{++}$', 'S$^{2+}$', 'O$^{+}$(hot)', 'S$^{3+}$']

    # Plotting densities for all species, in multiplot panels
    for field_line_index in field_line_indices:
        field_label = field_line_labels[field_line_index]
        fig, axs = plt.subplots(5, 2, figsize=(14, 16), sharex=True)
        fig.suptitle(f'Densities along {field_label} (Northern Hemisphere)', fontsize=20, fontweight='bold')
        for idx_sp, (key, label) in enumerate(zip(species_keys, species_labels)):
            row = idx_sp // 2
            col = idx_sp % 2
            ax = axs[row, col]
            for case_name in case_names:
                s = results[case_name][field_line_index]['s']
                n_out_line = results[case_name][field_line_index]['n_out_line'][key]
                # Only consider s >= 0 (northern hemisphere)
                positive_s_indices = s >= 0
                s_positive = s[positive_s_indices]
                n_positive = n_out_line[positive_s_indices]
                # For Europa, limit s ≤ 11.39
                if field_line_index == 538:
                    s_limit_indices = s_positive <= 11.39
                    s_positive = s_positive[s_limit_indices]
                    n_positive = n_positive[s_limit_indices]
                else:
                    s_limit_indices = s_positive <= 7.2
                    s_positive = s_positive[s_limit_indices]
                    n_positive = n_positive[s_limit_indices]
                # Exclude zero densities (for hot species in Cases C-F)
                if np.any(n_positive > 0):
                    ax.plot(s_positive, n_positive, label=f"Case {case_name}")
            #ax.set_ylabel('Density', fontsize=14, fontweight='bold')
            #ax.text(0.5, 0.9, f'{label}', transform=ax.transAxes, fontsize=22, fontweight='bold', ha='center')
            if (key == 'eh'):
                ax.text(0.5, 0.7, f'{label}', transform=ax.transAxes, fontsize=22, fontweight='bold', ha='center')
            else:
                ax.text(0.5, 0.8, f'{label}', transform=ax.transAxes, fontsize=22, fontweight='bold', ha='center')
            ax.grid(True, which='both', linestyle='--')
            ax.tick_params(axis='both', which='both', labelsize=17, width=2)
            for tick in ax.get_xticklabels():
                tick.set_fontweight('bold')
            for tick in ax.get_yticklabels():
                tick.set_fontweight('bold')
            if row == 3:
                ax.set_xlabel(r's ($\bf{R_J}$)', fontsize=20, fontweight='bold')
        fig.text(0.04, 0.5, 'Density (cm$^{-3}$)', va='center', rotation='vertical', fontsize=20, fontweight='bold')
        axs[3,1].legend(loc='center right', frameon=False, prop={'weight': 'bold', 'size':16})
        plt.tight_layout(rect=[0.05, 0.03, 1, 0.95])
        plt.savefig(f'Densities_vs_s_Species_FieldLine_{field_label.replace(" ", "_")}_NorthernHemisphere.png', dpi=600)
        plt.show()
    
    # Plotting Phi in volts over s, for all cases
    for field_line_index in field_line_indices:
        field_label = field_line_labels[field_line_index]
        fig, ax = plt.subplots(figsize=(10, 6))
        for case_name in case_names:
            s = results[case_name][field_line_index]['s']
            phis = results[case_name][field_line_index]['phis']
            # Only consider s >= 0 (northern hemisphere)
            positive_s_indices = s >= 0
            s_positive = s[positive_s_indices]
            phis_positive = phis[positive_s_indices]
            # For Europa, limit s ≤ 11.39
            if field_line_index == 538:
                s_limit_indices = s_positive <= 11.39
                s_positive = s_positive[s_limit_indices]
                phis_positive = phis_positive[s_limit_indices]
            else:
                s_limit_indices = s_positive <= 7.2
                s_positive = s_positive[s_limit_indices]
                phis_positive = phis_positive[s_limit_indices]
            ax.plot(s_positive, phis_positive, label=f"Case {case_name}")
        ax.set_xlabel(r's ($\bf{R_J}$)', fontsize=20, fontweight='bold')
        ax.set_ylabel(r'$\Phi$ (V)', fontsize=20, fontweight='bold')
        ax.set_title(f'$\Phi$ along {field_label} (Northern Hemisphere)', fontsize=20, fontweight='bold')
        ax.grid(True, which='both', linestyle='--')
        ax.tick_params(axis='both', which='both', labelsize=17, width=2)
        for tick in ax.get_xticklabels():
            tick.set_fontweight('bold')
        for tick in ax.get_yticklabels():
            tick.set_fontweight('bold')
        ax.legend(frameon=False, prop={'weight': 'bold', 'size':16})
        plt.tight_layout()
        plt.savefig(f'Phi_vs_s_FieldLine_{field_label.replace(" ", "_")}_NorthernHemisphere.png', dpi=600)
        plt.show()
        
    # ... [Rest of your code above]
    
    # Save data to CSV files for further use in Python or Mathematica
    
    # Species keys and labels
    species_keys = ['elec', 'eh', 'op', 'sp', 'o2p', 's2p', 'oph', 's3p', 'hp', 'nap']
    species_labels = ['e-', 'e-(hot)', 'O+', 'S+', 'O++', 'S$^{2+}$', 'O+(hot)', 'S$^{3+}$', 'H+', 'Na+']
    
    for field_line_index in field_line_indices:
        field_label = field_line_labels[field_line_index].replace(' ', '_')
    
        # Set s_max based on field line index
        if field_line_index == 189:  # Io Field Line
            s_max = 7.2
        elif field_line_index == 538:  # Europa Field Line
            s_max = 11.39
        else:
            s_max = np.inf  # Default value if other field lines are added
    
        for case_name in case_names:
            s = results[case_name][field_line_index]['s']  # Array of length N
            phis = results[case_name][field_line_index]['phis']
            x = results[case_name][field_line_index]['x']
            y = results[case_name][field_line_index]['y']
            z = results[case_name][field_line_index]['z']
            B = results[case_name][field_line_index]['B_out_line']
    
            # Create a boolean mask to limit s to the specified range
            s_mask = (s >= 0) & (s < s_max)
    
            # Apply the mask to all data arrays
            s = s[s_mask]
            phis = phis[s_mask]
            x = x[s_mask]
            y = y[s_mask]
            z = z[s_mask]
            B = B[s_mask]
    
            # Initialize list of data columns
            data_columns = [s, phis, x, y, z, B]
    
            # Initialize list of headers
            headers = ['s_RJ', 'phi_V', 'x_RJ', 'y_RJ', 'z_RJ', 'B_T']
    
            # Densities
            for key, label in zip(species_keys, species_labels):
                n_out_line = results[case_name][field_line_index]['n_out_line'][key]
                n_out_line = n_out_line[s_mask]  # Apply the mask
                data_columns.append(n_out_line)
                headers.append(f'density_{label}_cm3')
    
            # T_parallel
            for key, label in zip(species_keys, species_labels):
                T_para = results[case_name][field_line_index]['T_para'][key]
                T_para = T_para[s_mask]
                data_columns.append(T_para)
                headers.append(f'T_para_{label}_eV')
    
            # T_perp
            for key, label in zip(species_keys, species_labels):
                T_perp = results[case_name][field_line_index]['T_perp'][key]
                T_perp = T_perp[s_mask]
                data_columns.append(T_perp)
                headers.append(f'T_perp_{label}_eV')
    
            # Stack all data columns horizontally to form a 2D array
            data_array = np.column_stack(data_columns)  # Shape: (M, num_columns), where M <= N
    
            # Create header string
            header_str = ','.join(headers)
    
            # Save to CSV using np.savetxt
            filename = f'Data_FieldLine_{field_label}_Case_{case_name}.csv'
            np.savetxt(filename, data_array, delimiter=',', header=header_str, comments='', fmt='%g')
            print(f"Data saved to {filename}")

        
   

