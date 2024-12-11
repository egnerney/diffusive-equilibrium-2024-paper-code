# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 22:46:59 2024

@author: Owner
"""

# -*- coding: utf-8 -*-
"""
Created on Sun Oct  6 22:38:05 2024

@author: Owner
"""

import numpy as np
from scipy.optimize import bisect
from scipy.interpolate import interp1d
import sys
import multiprocessing as mp
# import matplotlib.pyplot as plt
# import matplotlib.colors as colors
# from scipy.ndimage import uniform_filter1d  # For smoothing
# import scipy.special
# import matplotlib.gridspec as gridspec
# from scipy.special import gamma, gammaincc
# from scipy.integrate import quad
import mpmath as mp_mp

# Set mpmath precision
mp_mp.mp.dps = 15  # 15 decimal places

# Define the Species class
class Species:
    def __init__(self, name, m, q, T, A, kappa, n, dist_type):
        self.name = name
        self.m = m
        self.q = q
        self.T = T
        self.A = A
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
    """
    Computes the density for a Maxwellian distribution. f ∝ Exp[-(m*vperp^2 + m*vpar^2)/(2*T)]

    Parameters:
    - species: The species object containing properties like n, m, q, T
    - deltaU: The difference in gravitational-centrifugal potential energy per unit mass
    - Phi: Electrostatic potential difference
    - Bratio: The ratio B / B0 (not used in this function but kept for consistency)
    """
    n = species.n * np.exp(-(species.m * deltaU + species.q * Phi) / species.T)
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
    'Fried_Egg': fried_egg_density
}

# Main function to compute densities for a single field line
def process_field_line(args):
    (j, species0_in, planet, fcor, if_psh, x_fl, y_fl, z_fl, rho_fl, lat_fl, phi_fl_wlong, r_fl, Btot_fl, r_00,
     nop0, no2p0, nsp0, ns2p0, ns3p0, nhp0, nnap0, noph0, nec0, neh0,
     ti0, thp0, tec0, kappa_Te, kappa_Top, kappa_e, kappa_op, kappa_Thp) = args

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

        # Interpolate parameters at equator
        Ticnew = interp1d(r_00, kappa_Top, fill_value="extrapolate")(r_fl_ceq)
        Ticold = interp1d(r_00, ti0, fill_value="extrapolate")(r_fl_ceq)
        Thpnew = interp1d(r_00, kappa_Thp, fill_value="extrapolate")(r_fl_ceq)
        Thpold = interp1d(r_00, thp0, fill_value="extrapolate")(r_fl_ceq)
        Tecnew = interp1d(r_00, kappa_Te, fill_value="extrapolate")(r_fl_ceq)
        Tecold = interp1d(r_00, tec0, fill_value="extrapolate")(r_fl_ceq)

        Aic = Ticnew / Ticold  # perp/par = kappa/max for fried egg
        Ahp = Thpnew / Thpold
        Ae = Tecnew / Tecold
        Ticnew = Ticold.item()
        Thpnew = Thpold.item()
        Tecnew = Tecold.item()

        T_values = [Ticnew,  # O+
                    Ticnew,  # O2+
                    Ticnew,  # S+
                    Ticnew,  # S2+
                    Ticnew,  # S3+
                    Thpnew,  # H+
                    Ticnew,  # Na+
                    10.,     # O+(hot)
                    10.,     # eh-
                    Tecnew   # e-
                    ]

        A_values = [Aic,   # O+
                    Aic,   # O2+
                    Aic,   # S+
                    Aic,   # S2+
                    Aic,   # S3+
                    Ahp,   # H+
                    Aic,   # Na+
                    1.,    # O+(hot)
                    1.,    # eh-
                    Ae     # e-
                    ]

        kappainew = interp1d(r_00, kappa_op, fill_value="extrapolate")(r_fl_ceq)
        kappaenew = interp1d(r_00, kappa_e, fill_value="extrapolate")(r_fl_ceq)
        kappainew = kappainew.item()
        kappaenew = kappaenew.item()
        kappa_values = [kappainew,  # O+
                        kappainew,  # O2+
                        kappainew,  # S+
                        kappainew,  # S2+
                        kappainew,  # S3+
                        kappainew,  # H+
                        kappainew,  # Na+
                        kappainew,  # O+(hot)
                        kappainew,  # eh-
                        kappaenew   # e-
                        ]

        nopnew = interp1d(r_00, nop0, fill_value="extrapolate")(r_fl_ceq)
        no2pnew = interp1d(r_00, no2p0, fill_value="extrapolate")(r_fl_ceq)
        nspnew = interp1d(r_00, nsp0, fill_value="extrapolate")(r_fl_ceq)
        ns2pnew = interp1d(r_00, ns2p0, fill_value="extrapolate")(r_fl_ceq)
        ns3pnew = interp1d(r_00, ns3p0, fill_value="extrapolate")(r_fl_ceq)
        nnapnew = interp1d(r_00, nnap0, fill_value="extrapolate")(r_fl_ceq)
        nhpnew = interp1d(r_00, nhp0, fill_value="extrapolate")(r_fl_ceq)
        nophnew = interp1d(r_00, noph0, fill_value="extrapolate")(r_fl_ceq)
        necnew = interp1d(r_00, nec0, fill_value="extrapolate")(r_fl_ceq)
        nehnew = interp1d(r_00, neh0, fill_value="extrapolate")(r_fl_ceq)

        nenew = necnew + nehnew

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

        nspnew = nspnew.item()
        ns3pnew = ns3pnew.item()
        no2pnew = no2pnew.item()
        nhpnew = nhpnew.item()
        nnapnew = nnapnew.item()

        # Create the array of densities (n values) for each species
        n_values = [noptotnew,  # O+
                    no2pnew,    # O2+
                    nspnew,     # S+
                    ns2ptotnew, # S2+
                    ns3pnew,    # S3+
                    nhpnew,     # H+
                    nnapnew,    # Na+
                    0.0,        # O+(hot)
                    0.0,        # eh-
                    nenew       # e-
                    ]

        # Update the temperatures (T), densities (n), anisotropy (A), and kappa for all species
        species_list = []
        for ii, sp in enumerate(species0_in):
            sp_new = Species(
                name=sp.name,
                m=sp.m,
                q=sp.q,
                T=T_values[ii],
                A=A_values[ii],
                kappa=kappa_values[ii],
                n=n_values[ii],
                dist_type=sp.type
            )
            species_list.append(sp_new)

        # Prepare arrays for computation
        n_ions = np.zeros((nspec, npoints))
        phis = np.zeros(npoints)
        deltaUs = np.zeros(npoints)

        # Magnetic field ratio
        Bratio = B_int_fl / B_fl_ceq
        phi = None
        posc = None
        nchoke = None
        phic = None

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
            else:
                n_ions[:, i] = np.zeros(nspec)
                phis[i] = np.nan

        # Assign computed densities to n_out_line
        n_out_line['op'][:] = n_ions[0, :]
        n_out_line['o2p'][:] = n_ions[1, :]
        n_out_line['sp'][:] = n_ions[2, :]
        n_out_line['s2p'][:] = n_ions[3, :]
        n_out_line['s3p'][:] = n_ions[4, :]
        n_out_line['hp'][:] = n_ions[5, :]
        n_out_line['nap'][:] = n_ions[6, :]
        n_out_line['oph'][:] = n_ions[7, :]
        n_out_line['eh'][:] = n_ions[8, :]
        n_out_line['elec'][:] = n_ions[9, :]

        T_out_line['elec'][:] = interp1d(r_00, tec0, fill_value="extrapolate")(r_fl_ceq)
        T_out_line['eh'][:] = interp1d(r_00, neh0, fill_value="extrapolate")(r_fl_ceq)
        T_out_line['sp'][:] = interp1d(r_00, ti0, fill_value="extrapolate")(r_fl_ceq)
        T_out_line['oph'][:] = interp1d(r_00, tec0, fill_value="extrapolate")(r_fl_ceq)
        T_out_line['hp'][:] = interp1d(r_00, thp0, fill_value="extrapolate")(r_fl_ceq)

        x_out_line[:] = x_int_fl
        y_out_line[:] = y_int_fl
        z_out_line[:] = z_int_fl
        B_out_line[:] = B_int_fl

    else:
        # Assign zeros if outside range
        for key in densities_keys:
            n_out_line[key][:] = 0.
            T_out_line[key][:] = 0.

        x_out_line[:] = x_int_fl
        y_out_line[:] = y_int_fl
        z_out_line[:] = z_int_fl
        B_out_line[:] = B_int_fl

    return (j, n_out_line, T_out_line, x_out_line, y_out_line, z_out_line, B_out_line)

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
    # Return densities and phi
    return n, phi, deltaU

if __name__ == '__main__':
    # Define species
    planet = 'Jupiter'
    fcor = 1.0  # azimuthal velocity relative to corotation
    if_psh = 0
    RP, Omega, GM = define_planet(planet)
    Omega *= fcor

    species_names = ['O+', 'O++', 'S+', 'S++', 'S+++', 'H+', 'Na+', 'O+(hot)', 'eh-', 'e-']
    nspec = len(species_names)
    species_m = [16.0, 16.0, 32.0, 32.0, 32.0, 1.0, 23.0, 16.0, 1.0 / 1837.0, 1.0 / 1837.0]  # AMU
    species_q = [1.0, 2.0, 1.0, 2.0, 3.0, 1.0, 1.0, 1.0, -1.0, -1.0]  # Elementary charges
    species_T = [79.3, 79.3, 79.3, 79.3, 79.3, 79.3, 94.1, 362.0, 46.0, 4.6]  # eV
    species_n = [592.0, 76.3, 163.0, 538.0, 90.7, 50.6, 97.2, 134.0, 2.5375, 2537.5]  # Density units
    species_A = [1.0] * nspec
    species_kappa = [100.0] * nspec
    species_type = ['Fried_Egg'] * nspec

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
            n=species_n[i],
            dist_type=species_type[i]
        )
        species0_in.append(sp)

    nrfls = 601
    nphifls = 1

    # Load data files (ensure these files exist)
    x_fl = np.transpose(np.loadtxt('int_aligned_x_jrm33+con20integral_601_4-10_properly_aligned_phi0=75.85_perfectaligned_for_4.70-9.71.csv', delimiter=','))
    y_fl = np.transpose(np.loadtxt('int_aligned_y_jrm33+con20integral_601_4-10_properly_aligned_phi0=75.85_perfectaligned_for_4.70-9.71.csv', delimiter=','))
    z_fl = np.transpose(np.loadtxt('int_aligned_z_jrm33+con20integral_601_4-10_properly_aligned_phi0=75.85_perfectaligned_for_4.70-9.71.csv', delimiter=','))
    rho_fl = np.sqrt(x_fl ** 2. + y_fl ** 2.)  # rho systemIII RJ
    r_fl = np.sqrt(x_fl ** 2. + y_fl ** 2. + z_fl ** 2.)  # spherical r systemIII RJ
    lat_fl = np.degrees(np.arcsin(z_fl / r_fl))  # sysIII lat degrees
    phi_fl = np.degrees(np.arctan2(y_fl, x_fl))  # elong degrees
    phi_fl_wlong = 360. - phi_fl

    r_00 = 0.01 * np.arange(601) + 4.0

    dlatfl = 0.1
    lat_int_fl_max = 70.0
    lat_int_fl_min = -lat_int_fl_max
    npointsfieldline = round((lat_int_fl_max - lat_int_fl_min) / dlatfl) + 1
    lat_int_fl = lat_int_fl_min + dlatfl * np.arange(npointsfieldline)

    # Define the 'n_out' and 'T_out' as dictionaries of numpy arrays
    densities_keys = ['op', 'o2p', 'sp', 's2p', 's3p', 'hp', 'nap', 'oph', 'eh', 'elec']
    n_out = {key: np.zeros((nrfls * nphifls, npointsfieldline)) for key in densities_keys}
    T_out = {key: np.zeros((nrfls * nphifls, npointsfieldline)) for key in densities_keys}

    # Create double precision arrays
    rho_out = np.zeros((nrfls * nphifls, npointsfieldline), dtype=np.float64)
    z_out = np.zeros((nrfls * nphifls, npointsfieldline), dtype=np.float64)
    x_out = np.zeros((nrfls * nphifls, npointsfieldline), dtype=np.float64)
    y_out = np.zeros((nrfls * nphifls, npointsfieldline), dtype=np.float64)
    B_out = np.zeros((nrfls * nphifls, npointsfieldline), dtype=np.float64)

    # Load more data files (ensure these files exist)
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

    ti0 = np.loadtxt('Tic_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv', delimiter=',')
    thp0 = np.loadtxt('Thp_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv', delimiter=',')
    tec0 = np.loadtxt('Tec_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv', delimiter=',')

    Btot_fl = np.transpose(np.loadtxt('Btot_int_aligned_jrm33+con20integral_601_4-10_properly_aligned_phi0=75.85_perfectaligned_for_4.70-9.71.csv', delimiter=','))

    kappa_Te = np.loadtxt('new_nominal_model_Te_electrons_012_momentsequal_bimax_to_kappa.csv', delimiter=',')
    kappa_e = np.loadtxt('new_nominal_model_kappa_electrons_012_momentsequal_bimax_to_kappa.csv', delimiter=',')
    kappa_Top = np.loadtxt('Te_O_012_momentsequal_bimax_to_kappa.csv', delimiter=',')
    kappa_op = np.loadtxt('kappa_O_012_momentsequal_bimax_to_kappa.csv', delimiter=',')
    kappa_Thp = kappa_Top * (thp0 / ti0)

    # Prepare arguments for multiprocessing
    args_list = []
    for j in range(len(r_00)):
        args = (j, species0_in, planet, fcor, if_psh, x_fl, y_fl, z_fl, rho_fl, lat_fl, phi_fl_wlong, r_fl, Btot_fl, r_00,
                nop0, no2p0, nsp0, ns2p0, ns3p0, nhp0, nnap0, noph0, nec0, neh0,
                ti0, thp0, tec0, kappa_Te, kappa_Top, kappa_e, kappa_op, kappa_Thp)
        args_list.append(args)

    # Use multiprocessing Pool
    num_processes = mp.cpu_count()  # Use all available CPUs
    with mp.Pool(processes=num_processes) as pool:
        results = pool.map(process_field_line, args_list)

    # Collect results
    for result in results:
        j, n_out_line, T_out_line, x_out_line, y_out_line, z_out_line, B_out_line = result
        # Assign to main arrays
        for key in densities_keys:
            n_out[key][j, :] = n_out_line[key]
            T_out[key][j, :] = T_out_line[key]
        x_out[j, :] = x_out_line
        y_out[j, :] = y_out_line
        z_out[j, :] = z_out_line
        B_out[j, :] = B_out_line

    # Compute additional outputs
    r_out = np.sqrt(x_out ** 2. + y_out ** 2. + z_out ** 2.)
    rho_out = np.sqrt(x_out ** 2. + y_out ** 2.)
    lat_out_deg = np.degrees(np.arcsin(z_out / r_out))

    # Save n_out (dictionary of species densities)
    np.savez_compressed('new_nominal_model_n_out_nominal_model_4-10_aniso_fried_egg_all_nominal_model.npz', **n_out)

    # Save T_out (dictionary of species temperatures)
    np.savez_compressed('new_nominal_model_T_out_nominal_model_4-10_aniso_fried_egg_all_nominal_model.npz', **T_out)

    # Save x_out, y_out, z_out, B_out (NumPy arrays)
    np.savez_compressed('new_nominal_model_field_data_nominal_model_4-10_aniso_fried_egg_all_nominal_model.npz',
                        x_out=x_out, y_out=y_out, z_out=z_out, B_out=B_out)
