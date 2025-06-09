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
interp=np.interp
from scipy.optimize import bisect
from scipy.interpolate import CubicSpline, interp1d, PchipInterpolator
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

from scipy.special import hyp2f1, gamma, hyperu
# Set mpmath precision
mp_mp.mp.dps = 30  # 15 decimal places
# Set mpmath precision for special functions if needed

# CODATA RECOMMENDED VALUES OF THE FUNDAMENTAL PHYSICAL CONSTANTS: 2022, physics.nist.gov/constants
# Unnecessary precision
kg_per_u = 1.66053906892e-27  # Atomic mass unit (u or daltons, defined as 1/12 carbon 12 mass) in kg/u
u_per_kg = 1./kg_per_u        # u/kg
m_e = 5.485799090441e-4       # electron mass in u ~ 1.0/1822.9 u 

#Standard Periodic table values taking into account isotopic weighting and differences due to nuclear binding energy
m_H = 1.0078                  # Hydrogen mass in u
m_O = 15.999                  # Oxygen mass in u
m_S = 32.065                  # Sulfur mass in u
m_Na = 22.990                 # Sodium mass in u
# Even more Unwarranted Unnecessary precision (ignoring electron binding energy  decrease too)
m_Hp = m_H - m_e              # Sinlgy (Fully) Ionized Hydrogen or H^{+} or proton mass in u
m_Op = m_O - m_e              # Sinlgy Ionized Oxygen or O^{+} mass in u
m_O2p = m_O - 2.*m_e          # Doubly Ionized Oxygen or O^{++} mass in u
m_Sp = m_S - m_e              # Sinlgy Ionized Sulfur or S^{+} mass in u
m_S2p = m_S - 2.*m_e          # Doubly Ionized Oxygen or S^{++} mass in u
m_S3p = m_S - 3.*m_e          # Triply Ionized Oxygen or S^{+++} mass in u
m_Nap = m_Na - m_e            # Sinlgy Ionized Sodium or Na^{+} mass in u

ELEMENTARY_CHARGE = 1.602176634e-19  # Elementary charge in Coulombs
EV_TO_JOULE = ELEMENTARY_CHARGE   # Conversion from eV to Joules
JOULE_TO_EV = 1./ELEMENTARY_CHARGE 


###############################################################################
#                               CLASS DEFINITIONS
###############################################################################
class Species:
    """
    Container for plasma species parameters.

    Attributes:
    -----------
    name : str
        Species name (e.g., 'O+', 'e-', 'S++', etc.).
    m : float
        Mass in AMU (if prior to conversion) or in kg (if SI).
    q : float
        Charge in elementary charges (if prior to conversion) or in Coulombs (if SI).
    T : float
        Characteristic temperature in eV (if prior to conversion) or Joules (if SI).
    A : float
        Anisotropy factor A0 = T_perp / T_par at reference location, if relevant.
    kappa : float
        Kappa parameter. For Maxwellian, kappa -> infinity.
    lam : float
        λ = kappa_perp / kappa_par (only needed in product-kappa or more exotic models).
    n : float
        Number density at reference location (in cm^-3 or other consistent units).
    type : str
        Distribution type. Must match a key in density_functions dict (e.g., 'Maxwellian',
        'Aniso_Maxwellian', 'Aniso_kappa', 'Fried_Egg', etc.).
    """
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


###############################################################################
#                            HELPER FUNCTIONS
###############################################################################
def define_planet(planet):
    """
    Returns planet-specific parameters: radius (m), rotation rate (rad/s),
    and GM (m^3/s^2).

    Parameters
    ----------
    planet : str
        Name of the planet

    Returns
    -------
    RP : float
        Planetary Equatorial radius (1 bar level) in meters.
    Omega : float
        Rotation rate in rad/s. 
    GM : float
        Gravitational parameter G*M_planet in m^3/s^2.
    """

    
    if planet.lower() == 'jupiter':    # From NASA Jupiter Fact Sheet
        RP = 7.1492e7                  # Equatorial Radius (1 bar level) in meters
        Omega = 1.7585e-4              # Rotation rate in rad/s 2.*np.pi/(9.9250*3600) using sidereal rotation period
        GM = 1.26687e17                # G*M (m^3/s^2) or 126.687*10^6 km^3/s^2
        return RP, Omega, GM

    elif planet.lower() == 'saturn':   # From NASA Saturn Fact Sheet
        RP = 6.0268e7                  # Equatorial Radius (1 bar level) in meters
        Omega = 1.6379e-4              # Rotation rate in rad/s 2.*np.pi/(10.656*3600) using sidereal rotation period
        GM = 3.7931e16                 # G*M (m^3/s^2) or 37.931*10^6 km^3/s^2
        return RP, Omega, GM
    
    elif planet.lower() == 'earth':    # From NASA Earth Fact Sheet
        RP = 3.3781e6                  # Equatorial Radius (1 bar level) in meters
        Omega = 7.29210e-5             # Rotation rate in rad/s 2.*np.pi/(23.9345*3600) using sidereal rotation period
        GM = 3.9860e14                 # G*M (m^3/s^2) or 0.39860*10^6 km^3/s^2 
        return RP, Omega, GM
    
    elif planet.lower() == 'uranus':   # From NASA Uranus Fact Sheet
        RP = 2.5559e7                  # Equatorial Radius (1 bar level) in meters
        Omega = 1.012e-4               # Rotation rate in rad/s 2.*np.pi/(17.24*3600) using sidereal rotation period
        GM = 5.7940e15                 # G*M (m^3/s^2) or 5.7940*10^6 km^3/s^2
        return RP, Omega, GM
    
    elif planet.lower() == 'neptune':  # From NASA Neptune Fact Sheet
        RP = 2.4764e7                  # Equatorial Radius (1 bar level) in meters
        Omega = 1.083e-4               # Rotation rate in rad/s 2.*np.pi/(16.11*3600) using sidereal rotation period
        GM = 6.8351e15                 # G*M (m^3/s^2) or 6.8351*10^6 km^3/s^2
        return RP, Omega, GM
    
    else:
        print(f"Planet {planet} is not supported in define_planet(). Exiting.")
        sys.exit(1)


###############################################################################
#                  DENSITY FUNCTIONS FOR DIFFERENT DISTRIBUTIONS
###############################################################################
def maxwellian_density(species, deltaU, Phi, Bratio):
    """
    Isotropic Maxwellian distribution: n(s) ~ exp[-(m*deltaU + q*Phi)/T ].

    B-field ratio (Bratio) is not used in the isotropic Maxwellian formula,
    but is included for uniform interface.

    References:
    -----------
    Bagenal & Sullivan (1981), eq. for isotropic Maxwellian

    Parameters:
    -----------
    species : Species
        Contains n, m, q, T (SI or consistent units).
    deltaU : float
        Gravitational + centrifugal potential energy difference per unit mass.
    Phi : float
        Electrostatic potential difference (voltage).
    Bratio : float
        B(s) / B0. (Unused here but included for interface consistency.)

    Returns:
    --------
    n_local : float
        Number density at location s.
    """
    n_local = species.n * np.exp(-(species.m * deltaU + species.q * Phi) / species.T)
    return n_local


def aniso_maxwellian_density(species, deltaU, Phi, Bratio):
    """
    Anisotropic Maxwellian distribution (standard form), where T_par is constant
    along the field line, but T_perp(s) ~ 1 / [ A + (1 - A)*B0/B ].

    n(s) = n0 / [ A0 + (1 - A0)*(B0/B) ] * exp[ - (m*deltaU + q*Phi)/T_par ]

    References:
    -----------
    Huang & Birmingham (1992); Bagenal & Sullivan (1981)

    Parameters:
    -----------
    species : Species
        Must have .A = T_perp0 / T_par0 at reference location.
    deltaU : float
        Potential energy difference per unit mass (grav + cent).
    Phi : float
        Electrostatic potential difference.
    Bratio : float
        B(s)/B0.

    Returns:
    --------
    n_local : float
        Computed number density at s.
    """
    B0_over_B = 1.0 / Bratio
    denom = species.A + (1.0 - species.A) * B0_over_B

    if denom <= 0:
        print("Warning: Non-physical denom in aniso_maxwellian_density. Setting denom -> eps.")
        denom = np.finfo(float).eps

    factor_exp = np.exp(-(species.m * deltaU + species.q * Phi) / species.T)
    n_local = (species.n / denom) * factor_exp
    return n_local




def aniso_maxwellian_const_Tperp_density(species, deltaU, Phi, Bratio):
    """
    Alternate anisotropic Maxwellian variant where T_perp is also assumed constant
    (along with T_par). This leads to a different expression:

    n(s) ~ (B(s)/B0)^(1 - A0) * exp[-(m*deltaU + q*Phi)/T_par ].

    References:
    -----------
    Bagenal & Sullivan (1981), Bagenal (1994), etc.

    Parameters:
    -----------
    species : Species
        .A is T_perp/T_par (constant).
    deltaU : float
        Potential energy difference per unit mass.
    Phi : float
        Electrostatic potential difference.
    Bratio : float
        B(s)/B0.

    Returns:
    --------
    n_local : float
        Number density at s.
    """
    factor = Bratio**(1.0 - species.A)
    exponent_term = np.exp(-(species.m * deltaU + species.q * Phi) / species.T)
    n_local = species.n * factor * exponent_term
    return n_local


def aniso_kappa_density(species, deltaU, Phi, Bratio):
    """
    Standard anisotropic kappa distribution with single kappa (same parallel & perp):
      n(s) = [ n0 / ( A0 + (1-A0)(B0/B) ) ] * [ 1 + (m deltaU + q Phi) / (kappa T_par) ]^(0.5 - kappa)

    References:
    -----------
    Meyer-Vernet et al. (1995); Moncuquet et al. (2002)

    Parameters:
    -----------
    species : Species
        Has .kappa, .A, .T, .m, .q
    deltaU : float
        Potential difference per unit mass (grav + cent).
    Phi : float
        Electrostatic potential difference.
    Bratio : float
        B(s)/B0.

    Returns:
    --------
    n_local : float
        Number density at s.
    """
    B0_over_B = 1.0 / Bratio
    denom = species.A + (1.0 - species.A) * B0_over_B
    if denom <= 0.:
        print("Warning: Negative denom in aniso_kappa_density. Setting denom -> eps.")
        denom = np.finfo(float).eps

    PEratio = (species.m * deltaU + species.q * Phi) / (species.kappa * species.T)
    # For physical consistency in root-finding, clamp if PEratio <= -1
    if PEratio <= -1.:
        PEratio = -1.0 + np.finfo(float).eps
    

    
    exponent_factor = (1.0 + PEratio)**(0.5 - species.kappa)
    n_local = (species.n / denom) * exponent_factor
    
    return n_local


def aniso_product_kappa_density(species, deltaU, Phi, Bratio):
    """
    Product kappa distribution with possibly different kappa_par and kappa_perp:
       f ~ [1 + (m v_par^2)/(2 kappa_par T_par)]^(-kappa_par - a) *
            [1 + (m v_perp^2)/(2 kappa_perp T_perp)]^(-kappa_perp - b)
    The integral solution typically involves hypergeometric functions.

    This function includes a rough check for B/B0 >= 1. If B < B0, we set Bratio = 1 + eps
    as integral diverges for B/B0<1

    References:
    -----------
    Nerney (2025) 

    Parameters:
    -----------
    species : Species
        Has .kappa (interpreted as kappa_par), .lam ( = kappa_perp / kappa_par ), etc.
    deltaU : float
    Phi : float
    Bratio : float
        B(s)/B0.

    Returns:
    --------
    n_local : float
        Number density at s (via hypergeometric for small kappa; a simpler expression for large).
    """
    
    k0 = species.kappa
    k_perp0 = k0 * species.lam
    k_tot0 = k0 + k_perp0
    n0 = species.n
    T_par0 = species.T
    A0 = species.A
    #T_perp0 = T_par0*A0
    l0 = species.lam
    
    
    # For convergence of integral and positive density (physical)
    if Bratio < 1.0:
        Bratio = 1.0 + np.finfo(float).eps

    
    # For convergence of integral and positive density (physical)
    beta = (species.m * deltaU + species.q * Phi) / (k0 * T_par0 )
    if beta <= -1.:
        beta = -1.0 + np.finfo(float).eps


    
    #alpha = beta*k0 
    x = A0*(Bratio - 1.0)

    if Bratio != 1.0:
        
        # For smaller kappa (< ~17), attempt hypergeometric approach. Otherwise fallback.
        if (k0 < 17.) and (k_perp0 < 17.):
            delta = (l0 * A0) * ((Bratio - 1.) / (1. + beta))
            # Using 2F1 from scipy.special
 
            hf = hyp2f1(1., k_perp0 + 1.0, k_tot0  + 1.5, 1. - 1./delta)
    
            factor = k0 / (k_tot0  + 0.5)
            outer = (n0/ x) * ((1.0 + beta)**(0.5 - k0))
            n_local = Bratio *outer * factor * hf
    
        else:
            # Large kappa fallback: simpler expression ignoring the hypergeometric detail
            # Use approximation for large kappa_perp +k_par
            # k0*hf/(k0*(1+lambda0)+1/2) ~ 1 / (x + 1 + beta)
            # and good to a couple % for kappa_perp and kappa_par >17 (better for larger kappa)
            exponent_factor = Bratio * ((1. + beta)**(0.5 - k0))
            hypfac = 1. / (x + 1. + beta)
            n_local = n0 * exponent_factor * hypfac
    else:
        n_local = n0*((1. + beta)**(0.5 - k0))

    return n_local




def fried_egg_density(species, deltaU, Phi, Bratio):
    """
    "Fried-Egg" distribution: Maxwellian in parallel (kappa_par -> inf),
    kappa in the perpendicular direction. The velocity-space integral leads
    to exponential integrals E_n(x).

    n(s) ~ [some prefactor]* exp[-(m*deltaU + q*Phi)/T_par + ... ] * E_{kappa+1}(kappa * x)
    with x = A0*(B/B0 - 1).

    For physically consistent solutions in the derived formula, we
    require B >= B0.

    References:
    -----------
    Nerney (2025)

    Parameters:
    -----------
    species : Species
        Has .kappa (the perp kappa), .A = T_perp0/T_par0, .T = T_par in eV or J, etc.
    deltaU : float
        Potential difference per unit mass
    Phi : float
        Electrostatic potential difference
    Bratio : float
        B(s)/B0

    Returns:
    --------
    n_local : float
        Number density at s, using the "fried-egg" limit.
    """
    if Bratio < 1.0:
        print("Warning: B/B0 < 1 for fried_egg_density => non-convergent integral. Setting Bratio=1.")
        Bratio = 1.0 + np.finfo(float).eps
        
    k0 = species.kappa # k0 = kappa_0 = kappa_perp0 as there is no kappa_par0 for Fried Egg
    A0 = species.A # A0 = Tperp0/Tpar0
    x = A0 * (Bratio - 1.0)
    z = k0*x
    
    # exponent = \Delta PE / Tpar0
    exponent = -(species.m * deltaU + species.q * Phi) / species.T
    if np.abs(exponent) > 700:
        exponent = np.sign(exponent)*700.
    exp_factor = np.exp(exponent)
    
    # T assumed parallel T here
    # Though we can use mpmath approximations are highly accurate and faster
    # In addition scipy only allows positive integer n for E_n(x) 
    # and scipy special incomplete upper gamma function  only allows positive arguments so that respresntation won't work
    # Scipy special allows for real valued a, b, and x for hyperu(a,b,x) and checks against mpmath at dps=100 show good aggremment for kappa0 and x<30 
    # e^z * E_{k0+1}[z] = U[1, 1 - k0, z] 
    
    
    if Bratio != 1.0:
        # Small x or Bratio expansion to second order in x 
        # good to  0.9%. Much better away from k0 = 1
        # but near k0 = 1 this expansion doesn't do well so keep to within x< 0.001
        # Prevents nans from small Bratio>1 and x>0 in confluent hypergeometric scipy special call
        if x< 0.0001:
            factor = 1.
            n_local = species.n * Bratio * factor * exp_factor
        # if kappa=kappa_perp>30 and x>30. 0th order Approximation good to 0.1% that is same as Aniso-Maxwellian
        # But easier and better approximation to just use first order approximation for k0>15.
        elif k0 > 15.:
            # First order in 1/k0 expansion
            # k0 e^(k0*x)*E_{k0+1}[k0*x] ~1 / (1 + x) - x/(k0*(1+ x)^3) 
            # Good to 0.035% for k0>15. and all x values >0 => For Bratio>1
            factor = 1./(1. + x) - x / (k0*((x + 1.)**3.)) 
            n_local = species.n * Bratio * factor * exp_factor 
        elif (k0 < 15.) and (15.<x<30.):
            # Third order in 1/k0 expansion
            # k0 e^(k0*x)*E_{k0+1}[k0*x] ~1 / (1 + x) - x/(k0*(1+ x)^3) + x(2x-1)/(k0^2*(1+ x)^5)-  (6x^3 - 8x^2 +x)/(k0^3*(1+ x)^7)
            # Didn't expand in x but 1/k0 but still works but for large x that x^7 gets huge so use below for x>30
            # if x> 15. then for any k0 (though refer to above approximation if k0>15 for any x) 
            # Good to 0.02% for k0=1.001, better for larger k0 and x
            factor = 1./(1. + x) - x / (k0*((x + 1.)**3.)) + x*(2.*x - 1.)/((k0**2.)*((1.+ x)**5.)) -  (6.*(x**3.) - 8.*(x**2.) + x)/((k0**3.)*((1. + x)**7.))
            n_local = species.n * Bratio * factor * exp_factor
        elif (k0 < 15.) and (x>30.):
            # fourth order in 1/x expansion
            # For x>30 use this so that we don't have to calculate the x^7 above which gets large
            # if x> 30. then for any k0 (though refer to above approximation if k0>15) 
            # Good to 0.01% 
            numer = -6. + k0 *(-11. + 2.*x - (6. + (-3. + x)*x) * k0 + (-1. + x)*(1 + x**2.)*(k0**2.))
            denom = (x**4.)*(k0**3.)
            factor = numer/denom
            n_local = species.n * Bratio * factor * exp_factor
        else:
            func_value = k0 * hyperu(1., 1. - k0, z)
            n_local = species.n * Bratio * func_value * exp_factor
    else:
        # If Bratio = 1, x= 0, and we are at reference location s0
        # but just in case \Delta PE  does not equal 0 where Bratio = 1 
        # then we use limiting form for non zero \delta PE same as Isotropic Maxwellian
        n_local = species.n*exp_factor
        
 
        
    return n_local

###############################################################################
#                  TEMPERATURE FUNCTIONS FOR DIFFERENT DISTRIBUTIONS
###############################################################################

def maxwellian_temperature(species, deltaU, Phi, Bratio):
    """
    Temperature function for isotropic Maxwellian distribution.
    Temperature is constant along field line.
    
    Parameters:
    -----------
    species : Species
        Contains T (SI or consistent units).
    deltaU, Phi, Bratio : float
        Not used but included for interface consistency.
        
    Returns:
    --------
    T_par_local : float
        Parallel temperature at location s.
    T_perp_local : float
        Perpendicular temperature at location s.
    """
    T_par_local = species.T
    T_perp_local = species.T
    return T_par_local, T_perp_local


def aniso_maxwellian_temperature(species, deltaU, Phi, Bratio):
    """
    Temperature function for anisotropic Maxwellian distribution.
    T_parallel is constant, T_perpendicular varies with B.
    
    Parameters:
    -----------
    species : Species
        Contains T and A.
    deltaU, Phi, Bratio : float
        Bratio = B(s)/B0.
        
    Returns:
    --------
    T_par_local : float
        Parallel temperature at location s.
    T_perp_local : float
        Perpendicular temperature at location s.
    """
    T_par_local = species.T
    
    B0_over_B = 1.0 / Bratio
    denom = species.A + (1.0 - species.A) * B0_over_B
    
    if denom <= 0:
        print("Warning: Non-physical denom in aniso_maxwellian_temperature. Setting denom -> eps.")
        denom = np.finfo(float).eps
    
    T_perp_local = species.A * species.T / denom
    return T_par_local, T_perp_local

def aniso_maxwellian_const_Tperp_temperature(species, deltaU, Phi, Bratio):
    """
    Temperature function for anisotropic Maxwellian with constant T_perp.
    Both T_parallel and T_perpendicular are constant.
    
    Parameters:
    -----------
    species : Species
        Contains T and A.
    deltaU, Phi, Bratio : float
        Not used but included for interface consistency.
        
    Returns:
    --------
    T_par_local : float
        Parallel temperature at location s.
    T_perp_local : float
        Perpendicular temperature at location s.
    """
    T_par_local = species.T
    T_perp_local = species.A * species.T
    return T_par_local, T_perp_local


def aniso_kappa_temperature(species, deltaU, Phi, Bratio):
    """
    Temperature function for standard anisotropic kappa distribution.
    
    Parameters:
    -----------
    species : Species
        Contains T, A, kappa, m, q.
    deltaU : float
        Gravitational + centrifugal potential energy difference per unit mass.
    Phi : float
        Electrostatic potential difference (voltage).
    Bratio : float
        B(s) / B0.
        
    Returns:
    --------
    T_par_local : float
        Parallel temperature at location s.
    T_perp_local : float
        Perpendicular temperature at location s.
    """
    PEratio = (species.m * deltaU + species.q * Phi) / (species.kappa * species.T)
    
    # For physical consistency, clamp if PEratio <= -1
    if PEratio <= -1.:
        PEratio = -1.0 + np.finfo(float).eps
    
    beta = PEratio
    T_par_local = species.T * (1.0 + beta)
    
    B0_over_B = 1.0 / Bratio
    denom = species.A + (1.0 - species.A) * B0_over_B
    
    if denom <= 0:
        print("Warning: Non-physical denom in aniso_kappa_temperature. Setting denom -> eps.")
        denom = np.finfo(float).eps
    
    T_perp_local = species.A * species.T * (1.0 + beta) / denom
    return T_par_local, T_perp_local


def aniso_product_kappa_temperature(species, deltaU, Phi, Bratio):
    """
    Temperature function for anisotropic product kappa distribution.
    
    Parameters:
    -----------
    species : Species
        Contains T, A, kappa (interpreted as kappa_par), lam (kappa_perp/kappa_par).
    deltaU : float
        Gravitational + centrifugal potential energy difference per unit mass.
    Phi : float
        Electrostatic potential difference (voltage).
    Bratio : float
        B(s) / B0.
        
    Returns:
    --------
    T_par_local : float
        Parallel temperature at location s.
    T_perp_local : float
        Perpendicular temperature at location s.
    """
    # Check if B/B0 >= 1
    if Bratio < 1.0:
        print("Warning: B/B0 < 1 in aniso_product_kappa_temperature. Setting Bratio -> 1 + eps.")
        Bratio = 1.0 + np.finfo(float).eps
    
    beta = (species.m * deltaU + species.q * Phi) / (species.kappa * species.T)
    if beta <= -1.:
        beta = -1.0 + np.finfo(float).eps
    
    # Reference Values at s0
    kappa_par = species.kappa
    kappa_perp = kappa_par * species.lam
    A0 = species.A 
    T_par =  species.T
    T_perp =  A0*T_par
    x = A0 *(Bratio - 1.)
    g = (1.0 + beta)
    if Bratio != 1.0:
        # For larger kappa values, use simpler approximations
        if (kappa_par >= 17.) and (kappa_perp >= 17.):
            # Approximation for large kappa
            T_par_local = T_par * g 
            T_perp_local = Bratio* T_perp * g /(x + g )
        else:
            # Full calculation with hypergeometric functions
            x = species.A*(Bratio - 1.)
            delta_val = (species.lam * species.A) * ((Bratio - 1.) / (1. + beta))
            kappa_tot = kappa_par + kappa_perp
            
            # Calculate F_q values for Tperp
            # F_q = 2_F_1(q,kappa_perp + 0.5, kappa_tot + 1.5, 1.0 - 1./delta_val)
            F1 = hyp2f1(1.0, kappa_perp + 1.0, kappa_tot + 1.5, 1.0 - 1./delta_val)
            F2 = hyp2f1(2.0, kappa_perp + 1.0, kappa_tot + 1.5, 1.0 - 1./delta_val)
            F3q = hyp2f1(0.75, kappa_perp + 1.0, kappa_tot + 1.5, 1.0 - 1./delta_val)
            F7q = hyp2f1(1.75, kappa_perp + 1.0, kappa_tot + 1.5, 1.0 - 1./delta_val)
            
            
            # Calculate regularized 2_F_1 values for Tpar
            # H_q = 2_F_1(1.0, kappa_perp + 1.0, kappa_tot + 1.5, 1.0 - 1./delta_val)
            H1h = hyp2f1(1.0, kappa_perp + 1.0, kappa_tot + 0.5, 1.0 - 1./delta_val)/gamma(kappa_tot + 0.5)
            H3q = hyp2f1(1.0, kappa_perp + 1.0, kappa_tot + 0.75, 1.0 - 1./delta_val)/gamma(kappa_tot + 0.75)
            H3h = F1/gamma(kappa_tot + 1.5)
            H7q = hyp2f1(1.0, kappa_perp + 1.0, kappa_tot + 1.75, 1.0 - 1./delta_val)/gamma(kappa_tot + 1.75)
            
            # Calculate constants
    
            C1 = 4.0 * kappa_tot - 1.0
            C2 = 2.0 * kappa_tot - 1.0
            C3 = 4.0 * kappa_par - 1.0
            C4 = 2.0 * kappa_par - 1.0
    
            # Calculate temperatures
            T_perp_num = (beta + 1.0)*kappa_par*T_perp*(A0 + x)/(3. * A0 * x)
            
            Dperp1 =  C1*F3q / (3.*F7q)
            Dperp2 = C2*F1 / (2.*F2)
            
            T_perp_denom = Dperp1 - Dperp2
            
            # Value of Perp Temperature at s 
            T_perp_local = T_perp_num / T_perp_denom
    
    
            T_par_num = 8.0*(beta + 1.0)*kappa_par*T_par
            
            Dpar1 =  C3*C1*H7q / H3q
            Dpar2 = 2.*C4*C2*H3h / H1h
            T_par_denom = Dpar1 - Dpar2
            
            # Value of Parallel Temperature at s 
            T_par_local = T_par_num / T_par_denom
    else:
        # If B/B0 = 1 then we are at reference location s0
        T_perp_local = T_perp 
        T_par_local = T_par
    
    return T_par_local, T_perp_local



def fried_egg_temperature(species, deltaU, Phi, Bratio):
    """
    Temperature function for "Fried-Egg" distribution.
    
    Parameters:
    -----------
    species : Species
        Contains T, A, kappa (perp kappa).
    deltaU : float
        Gravitational + centrifugal potential energy difference per unit mass.
    Phi : float
        Electrostatic potential difference (voltage).
    Bratio : float
        B(s) / B0.
        
    Returns:
    --------
    T_par_local : float
        Parallel temperature at location s.
    T_perp_local : float
        Perpendicular temperature at location s.
    """
    # Check if B/B0 >= 1
    if Bratio <1.0:
        print("Warning: B/B0 < 1 in fried_egg_temperature. Setting Bratio -> 1 + eps.")
        Bratio = 1.0 + np.finfo(float).eps
    
    # T_parallel is constant for Fried-Egg
    T_par_local = species.T
    T_perp0 =  species.T *species.A
    
    k0 = species.kappa
    # For T_perpendicular, we use the derived formula with x and z
    x = species.A * (Bratio - 1.0)
    z = k0 * x
    
    
    
    if Bratio != 1.0:
        if k0 > 100.:
            # 0th order approximation is just
            # Tperp result is same as Aniso Maxwellian (Huang & Birmingham)
            # factor = 1./(1. + x) as Bratio/(1. + x) = 1/(A0+(1-A0)B0/B)
            # Good to 0.25 % for k0>100 for all x>0
            # Better for larger k0 and larger x 
            factor = 1./(1. + x) 
            T_perp_local = Bratio*T_perp0*factor   
        elif k0 > 15.:
            # Good to 0.05% for all x>0 
            # Better for larger kappa and larger x 
            factor = 1./(1. + x) - x/(((1. + x)**3.)*k0)
            T_perp_local =  Bratio*T_perp0*factor
        elif (x>7.) and (k0<15.):
            # For k0> 1, Bratio>1, and x>7. this 5th order approximation
            # Good to better than 0.22% for all all k0>1 
            # better for larger x or k0
            C1 = (3. + k0)*(4. + k0 * (7. + k0))
            C2 = -x*((1. + k0) ** 2.)*(4. + k0)
            C3 = -( x**3.) * (k0**2.) *(1. + k0) + ( x**2.) * k0 *(1. + k0) *(2. + k0)
            C4 = ( x**4.) * (k0**3.) + C3 + C2 + C1
            numerator =  -4. + k0 *C4
            denominator = ( x**5.) * (k0**4.)
            factor = numerator/denominator
            T_perp_local = Bratio*T_perp0*factor
        else:
            #for x<7 and k0 <15
            U5q = hyperu(k0 + 1., k0 + 5./4., z)
            U1q = hyperu(k0 + 1., k0 + 1./4., z)
            U1 = hyperu(k0 + 1., k0 + 1., z)
            U0 = hyperu(k0 + 1., k0, z)
            denominator = 4.* (U5q/U1q) - 3. *(U1/U0)
            factor = (1./(x*denominator))
            T_perp_local = Bratio * T_perp0 * factor
 
    else:
        # If B/B0 = 1 then we are at reference location s0
        T_perp_local = T_perp0 
    
    return T_par_local, T_perp_local




###############################################################################
#                  KAPPA FUNCTIONS FOR DIFFERENT DISTRIBUTIONS
###############################################################################

def maxwellian_kappa(species, deltaU, Phi, Bratio):
    """
    Kappa function for Maxwellian distribution.
    Not applicable but included for interface consistency.
    
    Returns:
    --------
    None, None : No kappa parameters for Maxwellian.
    Though actually effectivly Infinite as Maxwellian is infinite kappa limit
    """
    return None, None


def aniso_maxwellian_kappa(species, deltaU, Phi, Bratio):
    """
    Kappa function for anisotropic Maxwellian distribution.
    Not applicable but included for interface consistency.
    
    Returns:
    --------
    None, None : No kappa parameters for anisotropic Maxwellian.
    Though actually effectivly Infinite as Maxwellian is infinite kappa limit
    """
    return None, None


def aniso_maxwellian_const_Tperp_kappa(species, deltaU, Phi, Bratio):
    """
    Kappa function for anisotropic Maxwellian with constant T_perp.
    Not applicable but included for interface consistency.
    
    Returns:
    --------
    None, None : No kappa parameters for this distribution.
    Though actually effectivly Infinite as Maxwellian is infinite kappa limit
    """
    return None, None


def aniso_kappa_kappa(species, deltaU, Phi, Bratio):
    """
    Kappa function for standard anisotropic kappa distribution.
    Kappa parameter is constant along field line.
    
    Parameters:
    -----------
    species : Species
        Contains kappa.
    deltaU, Phi, Bratio : float
        Not used but included for interface consistency.
        
    Returns:
    --------
    kappa_par_local : float
        Parallel kappa parameter at location s.
    kappa_perp_local : float
        Perpendicular kappa parameter at location s.
    there is no par or perp in standard kappa only kappa 
    but include both for interface consistency.
    """
    kappa_par_local = species.kappa
    kappa_perp_local = species.kappa  # Same kappa for parallel and perpendicular
    # but include both for interface consistency.
    return kappa_par_local, kappa_perp_local


def aniso_product_kappa_kappa(species, deltaU, Phi, Bratio):
    """
    Kappa function for anisotropic product kappa distribution.
    
    Parameters:
    -----------
    species : Species
        Contains kappa (kappa_par), lam (kappa_perp/kappa_par).
    deltaU : float
        Gravitational + centrifugal potential energy difference per unit mass.
    Phi : float
        Electrostatic potential difference (voltage).
    Bratio : float
        B(s) / B0.
        
    Returns:
    --------
    kappa_par_local : float
        Parallel kappa parameter at location s.
    kappa_perp_local : float
        Perpendicular kappa parameter at location s.
    """
    # Check if B/B0 >= 1
    if Bratio < 1.0:
        print("Warning: B/B0 < 1 in aniso_product_kappa_kappa. Setting Bratio -> 1 + eps.")
        Bratio = 1.0 + np.finfo(float).eps
    
    beta = (species.m * deltaU + species.q * Phi) / (species.kappa * species.T)
    if beta <= -1.:
        beta = -1.0 + np.finfo(float).eps
    
    # kappa values at reference location s0
    # kappa_par0
    # kappa_perp0
    kappa_par = species.kappa
    kappa_perp = kappa_par * species.lam
    
    if Bratio != 1.0:
        # For larger kappa values, return constant kappa 
        # (approximation for large kappa matches standard kappa behavior)
        if (kappa_par >= 30.) and (kappa_perp >= 30.):
            kappa_par_local = kappa_par
            kappa_perp_local = kappa_perp
        else:
            # Full calculation with hypergeometric functions
            delta_val = (species.lam * species.A) * ((Bratio - 1.) / (1. + beta))
            kappa_tot = kappa_par + kappa_perp
            
            # Calculate F_q values for Tperp
            # F_q = 2_F_1(q,kappa_perp + 0.5, kappa_tot + 1.5, 1.0 - 1./delta_val)
            F1 = hyp2f1(1.0, kappa_perp + 1.0, kappa_tot + 1.5, 1.0 - 1./delta_val)
            F2 = hyp2f1(2.0, kappa_perp + 1.0, kappa_tot + 1.5, 1.0 - 1./delta_val)
            F3q = hyp2f1(0.75, kappa_perp + 1.0, kappa_tot + 1.5, 1.0 - 1./delta_val)
            F7q = hyp2f1(1.75, kappa_perp + 1.0, kappa_tot + 1.5, 1.0 - 1./delta_val)
            
            
            # Calculate regularized 2_F_1 values for Tpar
            # H_q = 2_F_1(1.0, kappa_perp + 1.0, kappa_tot + 1.5, 1.0 - 1./delta_val)
            H1h = hyp2f1(1.0, kappa_perp + 1.0, kappa_tot + 0.5, 1.0 - 1./delta_val)/gamma(kappa_tot + 0.5)
            H3q = hyp2f1(1.0, kappa_perp + 1.0, kappa_tot + 0.75, 1.0 - 1./delta_val)/gamma(kappa_tot + 0.75)
            H3h = F1/gamma(kappa_tot + 1.5)
            H7q = hyp2f1(1.0, kappa_perp + 1.0, kappa_tot + 1.75, 1.0 - 1./delta_val)/gamma(kappa_tot + 1.75)
            
            # Calculate constants

            C1 = 4.0 * kappa_tot - 1.0
            C2 = 2.0 * kappa_tot - 1.0
            C3 = 4.0 * kappa_par - 1.0
            C4 = 2.0 * kappa_par - 1.0
             
            num_kappa_perp = 2.0 * C1* F3q * F2 / (C2* F1 * F7q) - 4.0 
            # Calculate kappa_perp at s
            kappa_perp_local = 1.0/num_kappa_perp + 1.0
            
            
            
            # Calculate kappa_par at s
            fac1 = C3*C1/(C4*C2)
            fac2 = fac1*H1h*H7q/(H3q*H3h)
            kappa_par_local = (fac2 -2.0) / (2*fac2 - 8.0 )
    else:
        # If Bratio = 1 then at reference location s0
        kappa_par_local, kappa_perp_local = kappa_par, kappa_perp
    
    return kappa_par_local, kappa_perp_local


def fried_egg_kappa(species, deltaU, Phi, Bratio):
    """
    Kappa function for "Fried-Egg" distribution.
    
    Parameters:
    -----------
    species : Species
        Contains kappa (perp kappa).
    deltaU : float
        Gravitational + centrifugal potential energy difference per unit mass.
    Phi : float
        Electrostatic potential difference (voltage).
    Bratio : float
        B(s) / B0.
        
    Returns:
    --------
    None : Since there's no parallel kappa 
    (Maxwellian in parallel direction so parallel kappa effectivly infinite).
    kappa_perp_local : float
        Perpendicular kappa parameter at location s.
    """
    # Check if B/B0 >= 1
    if Bratio <1.0:
        print("Warning: B/B0 < 1 in fried_egg_kappa. Setting Bratio -> 1 + eps.")
        Bratio = 1.0 + np.finfo(float).eps
    
    # No parallel kappa parameter (Maxwellian in parallel direction)
    kappa_par_local = None
    
    #k0 = kappa_perp0
    k0 = species.kappa 
    # For kappa_perpendicular, use the derived formula with x and z
    x = species.A * (Bratio - 1.0)
    z = k0 * x
 
    if Bratio != 1.0:
        if  (k0 > 50.0):
            kappa_perp_local = k0
        elif (x >15.0):
            #x>15.
            # Series expansion about x = infinity to 2nd order
            # Good to 0.4% for k0 >1 and x>15
            term1 = ((x**2.)* (k0**2.))/(1. + k0) 
            term2 =  (289.* (2. + k0))/(16. * x*k0 *(1. + k0))
            term3 = (x*k0 *(27. + 8.*k0))/(4.* (1. + k0))
            kappa_perp_local = term1 + term2 + term3
        else: 
            #x<15. and kappa_0 <15.0
            # or kappa_0 > 15
            # Via scipy special direct calculation using confluent hypergeometric function U(a,b,z)
            U0 = hyperu(k0 + 1., k0, z)
            U1 = hyperu(k0 + 1., k0 + 1., z)
            U1q = hyperu(k0 + 1., k0 + 1./4., z)
            U5q = hyperu(k0 + 1., k0 + 5./4., z)
            
            fac = U0*U5q/(U1q*U1)
            kappa_perp_local = (0.75 - fac)/(1.0 - fac)
    else:
        # If Bratio = 1 we are at reference location
        kappa_perp_local = species.kappa
    
    return kappa_par_local, kappa_perp_local


###############################################################################
#         DICTIONARIES OF DENSITY, TEMPERATURE, AND KAPPA FUNCTIONS
###############################################################################

density_functions = {
    'Maxwellian': maxwellian_density,
    'Aniso_Maxwellian': aniso_maxwellian_density,
    'Aniso_Maxwellian_const_Tperp': aniso_maxwellian_const_Tperp_density,
    'Aniso_kappa': aniso_kappa_density,
    'Aniso_product_kappa': aniso_product_kappa_density,
    'Fried_Egg': fried_egg_density
}

temperature_functions = {
    'Maxwellian': maxwellian_temperature, # Constant but included for consistent structures
    'Aniso_Maxwellian': aniso_maxwellian_temperature,
    'Aniso_Maxwellian_const_Tperp': aniso_maxwellian_const_Tperp_temperature, # Constant but included for consistent structures
    'Aniso_kappa': aniso_kappa_temperature,
    'Aniso_product_kappa': aniso_product_kappa_temperature,
    'Fried_Egg': fried_egg_temperature
}

kappa_functions = {
    'Maxwellian': maxwellian_kappa, #NA But included for consistency accross structures
    'Aniso_Maxwellian': aniso_maxwellian_kappa, #NA But included for consistency accross structures
    'Aniso_Maxwellian_const_Tperp': aniso_maxwellian_const_Tperp_kappa, #NA But included for consistency accross structures
    'Aniso_kappa': aniso_kappa_kappa, #Constant But included for consistency accross structures
    'Aniso_product_kappa': aniso_product_kappa_kappa, 
    'Fried_Egg': fried_egg_kappa
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

    # Prepare output dictionaries for temperatures & kappa 
    # (parallel & perpendicular) for each species key
    T_out_line = {}
    kappa_out_line = {}
    for key in densities_keys:
        T_out_line[key + "_par"] = np.zeros(npoints)
        T_out_line[key + "_perp"] = np.zeros(npoints)
        kappa_out_line[key + "_par"] = np.zeros(npoints)
        kappa_out_line[key + "_perp"] = np.zeros(npoints)

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


 
    latfl_temp = lat_fl[idx_finite, j]
    xfltemp = x_fl[idx_finite, j]
    yfltemp = y_fl[idx_finite, j]
    zfltemp = z_fl[idx_finite, j]
    Btotfltemp = Btot_fl[idx_finite, j]
    rhofltemp = rho_fl[idx_finite, j]
    rfltemp = r_fl[idx_finite, j]
    
    phiflwlongtemp = phi_fl_wlong[idx_finite, j]
    sorted_indices = np.argsort(latfl_temp)

    
    
    x_int_fl = interp1d(lat_fl[idx_finite, j], x_fl[idx_finite, j], fill_value="extrapolate")(lat_int_fl)
    y_int_fl = interp1d(lat_fl[idx_finite, j], y_fl[idx_finite, j], fill_value="extrapolate")(lat_int_fl)
    z_int_fl = interp1d(lat_fl[idx_finite, j], z_fl[idx_finite, j], fill_value="extrapolate")(lat_int_fl)
    rho_int_fl = interp1d(lat_fl[idx_finite, j], rho_fl[idx_finite, j], fill_value="extrapolate")(lat_int_fl)
    B_int_fl = interp1d(lat_fl[idx_finite, j], Btot_fl[idx_finite, j], fill_value="extrapolate")(lat_int_fl)
    r_int_fl = interp1d(lat_fl[idx_finite, j], r_fl[idx_finite, j], fill_value="extrapolate")(lat_int_fl)
    phi_int_fl_wlong = interp1d(lat_fl[idx_finite, j], phi_fl_wlong[idx_finite, j], fill_value="extrapolate")(lat_int_fl)
  
    """
    
    x_int_fl = PchipInterpolator(latfl_temp[sorted_indices], xfltemp[sorted_indices])(lat_int_fl)
    y_int_fl = PchipInterpolator(latfl_temp[sorted_indices], yfltemp[sorted_indices])(lat_int_fl)
    z_int_fl = PchipInterpolator(latfl_temp[sorted_indices], zfltemp[sorted_indices])(lat_int_fl)
    rho_int_fl = PchipInterpolator(latfl_temp[sorted_indices], rhofltemp[sorted_indices])(lat_int_fl)
    B_int_fl = PchipInterpolator(latfl_temp[sorted_indices], Btotfltemp[sorted_indices])(lat_int_fl)
    r_int_fl = PchipInterpolator(latfl_temp[sorted_indices], rfltemp[sorted_indices])(lat_int_fl)
    phi_int_fl_wlong = PchipInterpolator(latfl_temp[sorted_indices], phiflwlongtemp[sorted_indices])(lat_int_fl)

    
    x_int_fl = interp1d(lat_fl[idx_finite, j], x_fl[idx_finite, j], fill_value="extrapolate")(lat_int_fl)
    y_int_fl = interp1d(lat_fl[idx_finite, j], y_fl[idx_finite, j], fill_value="extrapolate")(lat_int_fl)
    z_int_fl = interp1d(lat_fl[idx_finite, j], z_fl[idx_finite, j], fill_value="extrapolate")(lat_int_fl)
    rho_int_fl = interp1d(lat_fl[idx_finite, j], rho_fl[idx_finite, j], fill_value="extrapolate")(lat_int_fl)
    B_int_fl = interp1d(lat_fl[idx_finite, j], Btot_fl[idx_finite, j], fill_value="extrapolate")(lat_int_fl)
    r_int_fl = interp1d(lat_fl[idx_finite, j], r_fl[idx_finite, j], fill_value="extrapolate")(lat_int_fl)
    phi_int_fl_wlong = interp1d(lat_fl[idx_finite, j], phi_fl_wlong[idx_finite, j], fill_value="extrapolate")(lat_int_fl)
   """
    
    # Create positions
    Pos = np.array([r_int_fl, lat_int_fl, phi_int_fl_wlong, B_int_fl]).T

    idx_ceq = np.argmax(rho_int_fl)
    r_fl_ceq = r_int_fl[idx_ceq]
    lat_fl_ceq = lat_int_fl[idx_ceq]
    longw_fl_ceq = phi_int_fl_wlong[idx_ceq]
    B_fl_ceq = B_int_fl[idx_ceq]
    
        
    phis = np.zeros(npoints)
    deltaUs = np.zeros(npoints)
    #Bratio = np.zeros(npoints)
    Bratio_vec = np.zeros(npoints)
    

    if (r_fl_ceq >= 4.00) and (r_fl_ceq <= 10.0):
        Pos0 = [r_fl_ceq, lat_fl_ceq, longw_fl_ceq, B_fl_ceq]

        # Interpolate parameters at equator
        Ticnew = interp1d(r_00, kappa_Top, fill_value="extrapolate")(r_fl_ceq)
        # Ticold = interp1d(r_00, ti0, fill_value="extrapolate")(r_fl_ceq)
        Thpnew = interp1d(r_00, kappa_Thp, fill_value="extrapolate")(r_fl_ceq)
        #Thpold = interp1d(r_00, thp0, fill_value="extrapolate")(r_fl_ceq)
        Tecnew = interp1d(r_00, kappa_Te, fill_value="extrapolate")(r_fl_ceq)
        # Tecold = interp1d(r_00, tec0, fill_value="extrapolate")(r_fl_ceq)

        Aic = 1. #Ticnew / Ticold  # perp/par = kappa/max for fried egg
        Ahp = Aic #Thpnew / Thpold
        Ae =  Aic #Tecnew / Tecold
        #Ticnew = Ticold.item()
        #Thpnew = Thpold.item()
        #Tecnew = Tecold.item()
        
        Ticnew = Ticnew / Aic
        
        Thpnew = Thpnew / Ahp
        
        Tecnew = Tecnew /  Ae


        lamic = 1.
        lame=1.
        lamhp = 1.
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
        
        lam_values =[lamic,   # O+
                    lamic,   # O2+
                    lamic,   # S+
                    lamic,   # S2+
                    lamic,   # S3+
                    lamhp,   # H+
                    lamic,   # Na+
                    1.,    # O+(hot)
                    1.,    # eh-
                    lame     # e-
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
                lam=lam_values[ii],
                n=n_values[ii],
                dist_type=sp.type
            )
            species_list.append(sp_new)

        # Prepare arrays for computation
        n_ions = np.zeros((nspec, npoints))


        # Magnetic field ratio
        #Bratio = B_int_fl / B_fl_ceq
        Bratio_vec = B_int_fl / B_fl_ceq 
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
                    Bratio_vec[i],
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
        # Assign computed densities to n_out_line
        for idx_sp, key in enumerate(densities_keys):
            n_out_line[key] = n_ions[idx_sp, :]
        
        kg_per_u = 1.66053906892e-27  # Atomic mass unit (u or daltons, defined as 1/12 carbon 12 mass) in kg/u
        ELEMENTARY_CHARGE = 1.602176634e-19  # Elementary charge in Coulombs
        species_list2 = []
        for sp in species_list:
            sp_new = Species(
                name=sp.name,
                m=sp.m * kg_per_u,    # Convert mass from AMU to kg
                q=sp.q * ELEMENTARY_CHARGE,     # Convert charge from elementary charges to Coulombs
                T=sp.T * ELEMENTARY_CHARGE,     # Convert temperature from eV to Joules
                A=sp.A,
                kappa=sp.kappa,
                lam=sp.lam,
                n=sp.n,
                dist_type=sp.type
            )
            species_list2.append(sp_new)
        # Now compute T_par, T_perp, kappa_par, kappa_perp at each point 
        # using the final potential (phis[i]) and deltaUs[i].
        for i in range(npoints):
            # If n_ions was zero, skip
            if np.all(n_ions[:, i] == 0.0):
                for key in densities_keys:
                    T_out_line[key + "_par"][i]  = 0.0
                    T_out_line[key + "_perp"][i] = 0.0
                    kappa_out_line[key + "_par"][i]  = 0.0
                    kappa_out_line[key + "_perp"][i] = 0.0
                continue
        
            # We found phi at this lat
            phi_local   = phis[i]
            deltaU_here = deltaUs[i]
            B_ratio     = Bratio_vec[i]
           
            # For each species, compute T_par, T_perp, kappa_par, kappa_perp
            for idx_sp, sp in enumerate(species_list2):
                # temperature_functions & kappa_functions 
                # call them in the same "units" the Species object is using
                T_par_loc, T_perp_loc = temperature_functions[sp.type](sp, deltaU_here, phi_local, B_ratio)
                kpar_loc, kperp_loc   = kappa_functions[sp.type](sp, deltaU_here, phi_local, B_ratio)

                # Store and put back in  eV
                species_key = densities_keys[idx_sp]
                T_out_line[species_key + "_par"][i]  = T_par_loc/ELEMENTARY_CHARGE
                T_out_line[species_key + "_perp"][i] = T_perp_loc/ELEMENTARY_CHARGE
                if kpar_loc is not None:
                    kappa_out_line[species_key + "_par"][i]  = kpar_loc
                else:
                    kappa_out_line[species_key + "_par"][i]  = 0.0
                if kperp_loc is not None:
                    kappa_out_line[species_key + "_perp"][i] = kperp_loc
                else:
                    kappa_out_line[species_key + "_perp"][i] = 0.0
        
        x_out_line[:] = x_int_fl
        y_out_line[:] = y_int_fl
        z_out_line[:] = z_int_fl
        B_out_line[:] = B_int_fl

    else:
        # Assign zeros if outside range
        for key in densities_keys:
            n_out_line[key][:] = 0.
            T_out_line[key + "_par"][:]  = 0.0
            T_out_line[key + "_perp"][:] = 0.0
            kappa_out_line[key + "_par"][:]  = 0.0
            kappa_out_line[key + "_perp"][:] = 0.0

        x_out_line[:] = x_int_fl
        y_out_line[:] = y_int_fl
        z_out_line[:] = z_int_fl
        B_out_line[:] = B_int_fl

    return (j, n_out_line, T_out_line, kappa_out_line, x_out_line, y_out_line, z_out_line, B_out_line, phis, deltaUs, Bratio_vec)


# Function to compute densities
def diff_eq_eddie(Pos_in, Pos0_in, Bratio, species0_in, phi, posc, nchoke, phic, planet, fcor, if_psh):
    # Set planet-specific variables
    RP, Omega, GM = define_planet(planet)
    Omega *= fcor
    kg_per_u = 1.66053906892e-27  # Atomic mass unit (u or daltons, defined as 1/12 carbon 12 mass) in kg/u
    ELEMENTARY_CHARGE = 1.602176634e-19  # Elementary charge in Coulombs
    

    # Convert species parameters to SI units
    species0 = []
    for sp in species0_in:
        sp_new = Species(
            name=sp.name,
            m=sp.m * kg_per_u,    # Convert mass from AMU to kg
            q=sp.q * ELEMENTARY_CHARGE,     # Convert charge from elementary charges to Coulombs
            T=sp.T * ELEMENTARY_CHARGE,     # Convert temperature from eV to Joules
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

    if (Pos0_in[0] == Pos_in[0] and
        Pos0_in[1] == Pos_in[1] and
        Pos0_in[2] == Pos_in[2]):
        # —— all three coordinates match exactly ——
        # If at reference location s=s0 so phi = 0.0 by definition
        # Consistent with density and temp variation along field line definitions
        Phi_root = 0.0
    else:
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
        # xtol of 1e-12 ~ neutrality to ~1 part in 10^11 for typical 6RJ paramters
        # Far more than needed but why not it's still fast
        Phi_root = bisect(net_charge_density, Phi, Phi2, xtol=1e-12, rtol=8.881784197001252e-16, maxiter=300)


    
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
        
    # Unnecessary precision
    kg_per_u = 1.66053906892e-27  # Atomic mass unit (u or daltons, defined as 1/12 carbon 12 mass) in kg/u
    u_per_kg = 1./kg_per_u        # u/kg
    m_e = 5.485799090441e-4       # electron mass in u ~ 1.0/1822.9 u 

    #Standard Periodic table values taking into account isotopic weighting and differences due to nuclear binding energy
    m_H = 1.0078                  # Hydrogen mass in u
    m_O = 15.999                  # Oxygen mass in u
    m_S = 32.065                  # Sulfur mass in u
    m_Na = 22.990                 # Sodium mass in u
    # Even more Unwarranted Unnecessary precision (ignoring electron binding energy  decrease too)
    m_Hp = m_H - m_e              # Sinlgy (Fully) Ionized Hydrogen or H^{+} or proton mass in u
    m_Op = m_O - m_e              # Sinlgy Ionized Oxygen or O^{+} mass in u
    m_O2p = m_O - 2.*m_e          # Doubly Ionized Oxygen or O^{++} mass in u
    m_Sp = m_S - m_e              # Sinlgy Ionized Sulfur or S^{+} mass in u
    m_S2p = m_S - 2.*m_e          # Doubly Ionized Oxygen or S^{++} mass in u
    m_S3p = m_S - 3.*m_e          # Triply Ionized Oxygen or S^{+++} mass in u
    m_Nap = m_Na - m_e            # Sinlgy Ionized Sodium or Na^{+} mass in u

    ELEMENTARY_CHARGE = 1.602176634e-19  # Elementary charge in Coulombs
    EV_TO_JOULE = ELEMENTARY_CHARGE   # Conversion from eV to Joules
    JOULE_TO_EV = 1./ELEMENTARY_CHARGE 

    #species_names = ['O+', 'O++', 'S+', 'S++', 'S+++', 'H+', 'Na+', 'O+(hot)', 'eh-', 'e-']
    # Making same as output for now consistency with output names but still 
    # for ease of use with already written code to use outputs
    # So change this instead of outputs
    species_names = ['op', 'o2p', 'sp', 's2p', 's3p', 'hp', 'nap', 'oph', 'eh', 'elec']
   
    
    nspec = len(species_names)
    
    #species_m = [16.0, 16.0, 32.0, 32.0, 32.0, 1.0, 23.0, 16.0, 1.0 / 1823., 1.0 / 1823.]  # AMU
    
    species_m =[m_Op, m_O2p, m_Sp, m_S2p, m_S3p, m_Hp, m_Nap, m_Op , m_e, m_e ]
    
    species_q = [1.0, 2.0, 1.0, 2.0, 3.0, 1.0, 1.0, 1.0, -1.0, -1.0]  # Elementary charges
    
    species_type = ['Fried_Egg'] * nspec
    #species_type = ['Aniso_kappa'] * nspec
    #species_type = ['Aniso_product_kappa'] * nspec
    
    #Pleaceholders will be filled later
    species_T = [79.3, 79.3, 79.3, 79.3, 79.3, 79.3, 94.1, 362.0, 46.0, 4.6]  # eV
    species_n = [592.0, 76.3, 163.0, 538.0, 90.7, 50.6, 97.2, 134.0, 2.5375, 2537.5]  # Density units
    species_A = [1.0] * nspec
    species_kappa = [100.0] * nspec
    species_lam = [1.] * nspec

    
    """
    density_functions = {
        'Maxwellian': maxwellian_density,
        'Aniso_Maxwellian': aniso_maxwellian_density,
        'Aniso_Maxwellian_const_Tperp': aniso_maxwellian_const_Tperp_density,
        'Aniso_kappa': aniso_kappa_density,
        'Aniso_product_kappa': aniso_product_kappa_density,
        'Fried_Egg': fried_egg_density
    }
    """

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

    nrfls = 601
    nphifls = 1
    

    # Load data files (ensure these files exist)
    x_fl = np.loadtxt('1000x601_shape_int_x_jrm33+con20integral_601_4-10_properly_aligned_whole_range_by_shifting_ends+-0.52_from_phi0=75.85_which_is_only_goodfor_4.70-9.71_without_shifts_at_ends.csv', delimiter=',')
    y_fl = np.loadtxt('1000x601_shape_int_y_jrm33+con20integral_601_4-10_properly_aligned_whole_range_by_shifting_ends+-0.52_from_phi0=75.85_which_is_only_goodfor_4.70-9.71_without_shifts_at_ends.csv', delimiter=',')
    z_fl = np.loadtxt('1000x601_shape_int_z_jrm33+con20integral_601_4-10_properly_aligned_whole_range_by_shifting_ends+-0.52_from_phi0=75.85_which_is_only_goodfor_4.70-9.71_without_shifts_at_ends.csv', delimiter=',')
    Btot_fl = np.loadtxt('1000x601_shape_int_Btot_jrm33+con20integral_601_4-10_properly_aligned_whole_range_by_shifting_ends+-0.52_from_phi0=75.85_which_is_only_goodfor_4.70-9.71_without_shifts_at_ends.csv', delimiter=',')
    

    
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
    #densities_keys = ['op', 'o2p', 'sp', 's2p', 's3p', 'hp', 'nap', 'oph', 'eh', 'elec']
    densities_keys = species_names
    n_out = {key: np.zeros((nrfls * nphifls, npointsfieldline)) for key in densities_keys}
    # Prepare output dictionaries for temperatures & kappa 
    # (parallel & perpendicular) for each species key
    # Create a dict with shape [nrfls*nphifls, npointsfieldline] for each 'species_component'
    T_out = {
        f"{species_key}_{component}": np.zeros((nrfls * nphifls, npointsfieldline)) 
        for species_key in densities_keys
        for component in ["par", "perp"]
    }
    
    kappa_out = {
        f"{species_key}_{component}": np.zeros((nrfls * nphifls, npointsfieldline)) 
        for species_key in densities_keys
        for component in ["par", "perp"]
    }

    # Create double precision arrays
    rho_out = np.zeros((nrfls * nphifls, npointsfieldline), dtype=np.float64)
    z_out = np.zeros((nrfls * nphifls, npointsfieldline), dtype=np.float64)
    x_out = np.zeros((nrfls * nphifls, npointsfieldline), dtype=np.float64)
    y_out = np.zeros((nrfls * nphifls, npointsfieldline), dtype=np.float64)
    B_out = np.zeros((nrfls * nphifls, npointsfieldline), dtype=np.float64)
    Bratio_out = np.zeros((nrfls * nphifls, npointsfieldline), dtype=np.float64)
    phi_out = np.zeros((nrfls * nphifls, npointsfieldline), dtype=np.float64)
    deltaU_out = np.zeros((nrfls * nphifls, npointsfieldline), dtype=np.float64)

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






    """
    #Standard Kappa Isotropic match numerically 0,1,2 moments to isotropic double Maxwellian Nominal Model
    kappa_Te = np.loadtxt('ceq_model_kappa_temp_elec_4-10RJ.csv', delimiter=',')
    kappa_e = np.loadtxt('ceq_model_kappa_elec_4-10RJ.csv', delimiter=',')
    kappa_Top = np.loadtxt('ceq_model_kappa_temp_ions_4-10RJ.csv', delimiter=',')
    kappa_op = np.loadtxt('ceq_model_kappa_ions_4-10RJ.csv', delimiter=',')

    
    """
    """
    #Product Kappa Isotropic match numerically 0,1,2 moments to isotropic double Maxwellian Nominal Model
    kappa_Te = np.loadtxt('ceq_model_product_kappa_temp_elec_4-10RJ.csv', delimiter=',')
    kappa_e = np.loadtxt('ceq_model_product_kappa_elec_4-10RJ.csv', delimiter=',')
    kappa_Top = np.loadtxt('ceq_model_product_kappa_temp_ions_4-10RJ.csv', delimiter=',')
    kappa_op = np.loadtxt('ceq_model_product_kappa_ions_4-10RJ.csv', delimiter=',')
    """

   

    # Fried-Egg Isotropic match numerically 0,1,2 moments to isotropic double Maxwellian Nominal Model
    kappa_Te = np.loadtxt('ceq_model_fried-egg_temp_elec_4-10RJ.csv', delimiter=',')
    kappa_e = np.loadtxt('ceq_model_fried-egg_kappa_elec_4-10RJ.csv', delimiter=',')
    kappa_Top = np.loadtxt('ceq_model_fried-egg_temp_ions_4-10RJ.csv', delimiter=',')
    kappa_op = np.loadtxt('ceq_model_fried-egg_kappa_ions_4-10RJ.csv', delimiter=',')

    
    
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

    #(j, n_out_line, T_out_line, kappa_out_line, x_out_line, y_out_line, z_out_line, B_out_line, phis, deltaUs, Bratio)
    # Collect results
    for result in results:
        j, n_out_line, T_out_line, kappa_out_line, x_out_line, y_out_line, z_out_line, B_out_line, phis, deltaUs, Bratio_vec = result
        # Assign to main arrays
        for key in densities_keys:
            n_out[key][j, :] = n_out_line[key]
            # Here we store the parallel and perpendicular temperatures 
            # into the pre-allocated arrays in T_out:
            T_out[f"{key}_par"][j, :] = T_out_line[f"{key}_par"]
            T_out[f"{key}_perp"][j, :] = T_out_line[f"{key}_perp"]
            kappa_out[f"{key}_par"][j, :] = kappa_out_line[f"{key}_par"]
            kappa_out[f"{key}_perp"][j, :] = kappa_out_line[f"{key}_perp"]
            
        x_out[j, :] = x_out_line
        y_out[j, :] = y_out_line
        z_out[j, :] = z_out_line
        B_out[j, :] = B_out_line
        Bratio_out[j, :] = Bratio_vec
        phi_out[j, :] = phis
        deltaU_out[j, :] = deltaUs

    # Compute additional outputs
    r_out = np.sqrt(x_out ** 2. + y_out ** 2. + z_out ** 2.)
    rho_out = np.sqrt(x_out ** 2. + y_out ** 2.)
    lat_out_deg = np.degrees(np.arcsin(z_out / r_out))

    # Save n_out (dictionary of species densities)
    np.savez_compressed('final_new_nominal_model_n_out_nominal_model_4-10_isoT_isokappa_A=1_fried_egg_all_nominal_model.npz', **n_out)

    # Save T_out (dictionary of species temperatures)
    np.savez_compressed('final_new_nominal_model_T_out_nominal_model_4-10_isoT_isokappa_A=1_fried_egg_all_nominal_model.npz', **T_out)

    # Save T_out (dictionary of species temperatures)
    np.savez_compressed('final_new_nominal_model_kappa_out_nominal_model_4-10_isoT_isokappa_A=1_fried_egg_all_nominal_model.npz', **kappa_out)


    # Save x_out, y_out, z_out, B_out (NumPy arrays)
    np.savez_compressed('final_new_nominal_model_field_data_nominal_model_4-10_isoT_isokappa_A=1_fried_egg_all_nominal_model.npz',
                        x_out=x_out, y_out=y_out, z_out=z_out, B_out=B_out, Bratio_out=Bratio_out, phi_out=phi_out, deltaU_out=deltaU_out)
