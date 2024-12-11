# -*- coding: utf-8 -*-
"""
Created on Sun Sep 29 16:59:39 2024

@author: Owner
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 27 08:59:59 2024

Modified to numerically solve for kappa and T for each species.

@author: edne8319
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gamma
from scipy.optimize import fsolve
from scipy.interpolate import interp1d

# Load data files (assuming they are in the same directory)
nop0 = np.loadtxt('nop_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv', delimiter=',')
no2p0 = np.loadtxt('no2p_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv', delimiter=',')
noph0 = np.loadtxt('noph_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv', delimiter=',')
nsp0 = np.loadtxt('nsp_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv', delimiter=',')
ns2p0 = np.loadtxt('ns2p_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv', delimiter=',')
ns3p0 = np.loadtxt('ns3p_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv', delimiter=',')
nhp0 = np.loadtxt('nhp_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv', delimiter=',')
nnap0 = np.loadtxt('nnap_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv', delimiter=',')
nec0 = np.loadtxt('nec_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv', delimiter=',')
neh0 = np.loadtxt('neh_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv', delimiter=',')

ti0 = np.loadtxt('Tic_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv', delimiter=',')
thp0 = np.loadtxt('Thp_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv', delimiter=',')
toph0 = np.loadtxt('Toph_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv', delimiter=',')
tec0 = np.loadtxt('Tec_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv', delimiter=',')
teh0 = np.loadtxt('Teh_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv', delimiter=',')

# Ensure all data arrays have the same length
num_points = len(tec0)

# Physical constants
m_e = 9.10938356e-31        # Electron mass in kg
e = 1.602176634e-19         # Elementary charge in C (J/eV)
pi = np.pi

# Ion masses in kg
m_O = 16. * 1.66053906660e-27        # Mass of O+ ion (16 u)
m_S = 32. * 1.66053906660e-27        # Mass of S++ ion (32 u)

# Arrays to store computed parameters
Te_electrons = np.zeros(num_points)
kappa_electrons = np.zeros(num_points)
Te_O = np.zeros(num_points)
kappa_O = np.zeros(num_points)
Te_S = np.zeros(num_points)
kappa_S = np.zeros(num_points)

# Radial positions from 4 to 10 RJ
R = np.linspace(4, 10, num_points)  # Radial distance array

alpha_ions = -0.92
beta_ions = -17.655

alpha_e = -7.221406#-4.07075
beta_e = -21.4182#-31.3772
gamma_e = -40.3772

feh_nom = neh0/(neh0 + nec0)
plt.figure(figsize=(10, 6))
plt.semilogy(R, neh0 + nec0, label='n_{e Total}')
plt.semilogy(R, neh0, label='n_{eh} Nominal Model')
plt.semilogy(R, 0.001*((R/6.)**6.)*(neh0 + nec0), label='n_{eh} for Feh0 = 0.001, fehpl = 6')
plt.semilogy(R, 0.001*((R/6.)**7.5)*(neh0 + nec0), label='n_{eh} for Feh0 = 0.001, fehpl = 7.5')
plt.semilogy(R, 0.002*((R/6.)**6.)*(neh0 + nec0), label='n_{eh} for Feh0 = 0.002, fehpl = 6')
plt.xlabel(r'Radial Distance ($R_J$)', fontsize=14)
plt.ylabel(r'n$_{e Total}$', fontsize=14)
plt.legend(fontsize=12)
plt.grid(True)
plt.tight_layout()
#plt.savefig('Teh='+f'{Teh_ev:.0f}'+', feh0=' + f'{feh0:.3f}' +'_fehpl='+f'{fehpl:.1f}'+'_'+'kappa_temperatures_vs_radial_distance_012_momentsequal_bimax_to_kappa.png', dpi=300)
plt.show()

plt.figure(figsize=(10, 6))
plt.semilogy(R, feh_nom, label='Nominal model')
plt.semilogy(R, 0.002*((R/6.)**6.), label='feh0=0.002, tecpl=6')
plt.semilogy(R, 0.001*((R/6.)**6.), label='feh0=0.001, tecpl=6')
plt.semilogy(R, 0.001*((R/6.)**7.), label='feh0=0.001, tecpl=7')
plt.semilogy(R, 0.001*((R/6.)**7.5), label='feh0=0.001, tecpl=7.5')
plt.semilogy(R, 0.001*((R/6.)**8.), label='feh0=0.001, tecpl=8')
plt.xlabel(r'Radial Distance ($R_J$)', fontsize=14)
plt.ylabel(r'F$_{eh}$', fontsize=14)
plt.legend(fontsize=12)
plt.grid(True)
plt.tight_layout()
#plt.savefig('Teh='+f'{Teh_ev:.0f}'+', feh0=' + f'{feh0:.3f}' +'_fehpl='+f'{fehpl:.1f}'+'_'+'kappa_temperatures_vs_radial_distance_012_momentsequal_bimax_to_kappa.png', dpi=300)
plt.show()

numnpoints = len(R)
nop_tot = np.zeros(num_points)
nop_h = np.zeros(num_points)
fop_h = np.zeros(num_points)

ns2p_tot = np.zeros(num_points)
ns2p_h = np.zeros(num_points)
fs2p_h = np.zeros(num_points)


feh_nom = np.zeros(num_points)
teh_nom = np.zeros(num_points)


# Loop over all radial positions
for i in range(num_points):
    # Now perform the same calculations for O⁺ ions
    # Calculate hot densities for O⁺ and S²⁺
    nop = nop0[i]
    ns2p = ns2p0[i]
    noph = noph0[i]

    if ns2p > 0.:
        ratio = nop / ns2p
        nophot = (noph * ratio) / (ratio + 1.)      # Hot O⁺ density
        ns2phot = noph  / (ratio + 1.)  # Hot S²⁺ density
    else:
        nophot = noph  # Hot O⁺ density
        ns2phot = 0.0  # Hot S²⁺ density
    nop_h[i] = nophot
    ns2p_h[i] = ns2phot
    nop_tot[i] = nop + nophot
    ns2p_tot[i] = ns2p + ns2phot
    fop_h[i] = nophot/(nop + nophot)
    fs2p_h[i] = ns2phot/(ns2p + ns2phot)
    
    if R[i] < 5.:
        fop_h[i] =0.
    if R[i] < 5.5:
        fs2p_h[i] =0.
        


plt.figure(figsize=(10, 6))
plt.semilogy(R, fop_h, label='Nominal model')
plt.xlabel(r'Radial Distance ($R_J$)', fontsize=14)
plt.ylabel(r'F$_{O^+}$', fontsize=14)
plt.legend(fontsize=12)
plt.grid(True)
plt.tight_layout()
#plt.savefig('Teh='+f'{Teh_ev:.0f}'+', feh0=' + f'{feh0:.3f}' +'_fehpl='+f'{fehpl:.1f}'+'_'+'kappa_temperatures_vs_radial_distance_012_momentsequal_bimax_to_kappa.png', dpi=300)
plt.show()

plt.figure(figsize=(10, 6))
plt.semilogy(R, fs2p_h, label='Nominal model')
plt.xlabel(r'Radial Distance ($R_J$)', fontsize=14)
plt.ylabel(r'F$_{S^{++}}$', fontsize=14)
plt.legend(fontsize=12)
plt.grid(True)
plt.tight_layout()
#plt.savefig('Teh='+f'{Teh_ev:.0f}'+', feh0=' + f'{feh0:.3f}' +'_fehpl='+f'{fehpl:.1f}'+'_'+'kappa_temperatures_vs_radial_distance_012_momentsequal_bimax_to_kappa.png', dpi=300)
plt.show()


# Loop over all radial positions
for i in range(num_points):
    if R[i] > 6.:
        kappa_guess_electrons = 3.#100.*((R[i]/6.) ** alpha_e)
        kappa_guess_s2p_and_Op = 4.*((R[i]/6.) ** alpha_ions)
    else:
        if R[i] > 5.8:
            kappa_guess_electrons = 3.#100.*((R[i]/6.) ** beta_e)
            kappa_guess_s2p_and_Op = 4.*((R[i]/6.) ** beta_ions)
        else:
            if R[i] > 5.5:
                kappa_guess_electrons = 3.#100.*((R[i]/6.) ** gamma_e)
                kappa_guess_s2p_and_Op = 4.*((R[i]/6.) ** beta_ions)
            else:
                kappa_guess_electrons = 3.# 500.
                kappa_guess_s2p_and_Op =4.

            
        
        
    
    # Electrons
    n_ec_cm3 = nec0[i]             # Core electron density in cm^-3
    n_eh_cm3 = neh0[i]             # Hot electron density in cm^-3
    T_ec_eV = tec0[i]              # Core electron temperature in eV
    #T_eh_eV = teh0[i]              # Hot electron temperature in eV

    # Convert densities to m^-3
    n_ec = n_ec_cm3 * 1e6       # Core electron density in m^-3
    n_eh = n_eh_cm3 * 1e6       # Hot electron density in m^-3
    n_e = n_ec + n_eh           # Total electron density in m^-3
    
    feh0 = 0.002
    fehpl = 6.
    teh0 = 270.
    tehpl = 1.5632

    
    if R[i] >= 5.7:
        feh = feh0*((R[i]/6.) ** fehpl) 
        T_eh_eV = teh0*((R[i]/6.) ** tehpl)
    else:
        feh = 0.0
        T_eh_eV = 1.0
        
    feh_nom[i] = feh
    teh_nom[i] = T_eh_eV


    
    #pppl = np.log(270./28.216762748831577)/np.log(6./5.7)
    #if R[i] >= 6.:
    #    T_eh_eV = teh0*((R[i]/6.) ** tehpl)
    #else:
    #    T_eh_eV = 28.216762748831577*((R[i]/5.7) ** pppl)
    
    #T_eh_eV = 500.
    Teh_ev = T_eh_eV
        
        
    n_eh = n_e*feh
    n_ec = n_e*(1. - feh)

    # Number density conservation
    n = n_e

    # Define T_mean and S for electrons
    T_mean = (n_ec * T_ec_eV + n_eh * T_eh_eV) / n
    S = (n_ec * np.sqrt(T_ec_eV) + n_eh * np.sqrt(T_eh_eV)) / n

    # Function to solve for electrons
    def equations_electron(vars):
        T_e_eV, kappa_e = vars
        f1 = T_e_eV - T_mean * (1. - 3. / (2. * kappa_e))
        f2 = S - np.sqrt(T_e_eV * kappa_e) * gamma(kappa_e - 1.) / gamma(kappa_e - 0.5)
        return [f1, f2]

    # Initial guesses for electrons
    """
    if i>0:
        initial_guess = [ Te_electrons[i-1],kappa_electrons[i-1]]
    else:
        initial_guess = [T_mean, 5.0]
    """
    
    if R[i]>7:
       initial_guess = [T_mean, 4.0]
    else:
        initial_guess = [T_ec_eV, 4.]#kappa_guess_electrons]
    #initial_guess = [T_ec_eV, kappa_guess_electrons]
        
    

    # Solve equations numerically
    sol_electron, info, ier, mesg = fsolve(equations_electron, initial_guess, full_output=True)

    T_e_eV, kappa_e = sol_electron

    # Check if solution is valid
    if ier == 1 and kappa_e > 1.505 and T_e_eV > 0:
        Te_electrons[i] = T_e_eV
        kappa_electrons[i] = kappa_e
        
    else:
        print('electron error', kappa_e,T_e_eV,'set to nans' )
        print(ier,mesg)
        Te_electrons[i] = np.nan
        kappa_electrons[i] = np.nan

    # Now perform the same calculations for O⁺ ions
    # Calculate hot densities for O⁺ and S²⁺
    nop = nop0[i]
    ns2p = ns2p0[i]
    noph = noph0[i]

    if ns2p > 0.:
        ratio = nop / ns2p
        nophot = (noph * ratio) / (ratio + 1.)      # Hot O⁺ density
        ns2phot = noph  / (ratio + 1.)  # Hot S²⁺ density
    else:
        nophot = noph  # Hot O⁺ density
        ns2phot = 0.0  # Hot S²⁺ density

    # For O⁺
    n_oc_cm3 = nop            # Cold O⁺ density in cm^-3
    n_oh_cm3 = nophot         # Hot O⁺ density in cm^-3
    n_O_cm3 = n_oc_cm3 + n_oh_cm3  # Total O⁺ density in cm^-3

    n_oc = n_oc_cm3 * 1e6     # Convert to m^-3
    n_oh = n_oh_cm3 * 1e6
    n_O = n_oc + n_oh

    T_oc_eV = ti0[i]           # Cold ion temperature in eV
    T_oh_eV = toph0[i]         # Hot ion temperature in eV

    # Number density conservation
    n = n_O

    # Define T_mean and S for O⁺ ions
    T_mean_O = (n_oc * T_oc_eV + n_oh * T_oh_eV) / n
    S_O = (n_oc * np.sqrt(T_oc_eV) + n_oh * np.sqrt(T_oh_eV)) / n

    # Function to solve for O⁺ ions
    def equations_O(vars):
        T_e_eV, kappa_e = vars
        f1 = T_e_eV - T_mean_O * (1. - 3. / (2. * kappa_e))
        f2 = S_O - np.sqrt(T_e_eV * kappa_e) * gamma(kappa_e - 1.) / gamma(kappa_e - 0.5)
        return [f1, f2]

    # Initial guesses for O⁺ ions
    """
    if i>0:
        initial_guess_O = [ Te_O[i-1],kappa_O[i-1]]
    else:
        initial_guess_O = [T_mean_O, 5.0]
    """
    #initial_guess_O = [T_oc_eV, kappa_guess_s2p_and_Op]
    initial_guess_O = [T_mean_O, 5.0]
    

    # Solve equations numerically
    sol_O, info, ier, mesg = fsolve(equations_O, initial_guess_O, full_output=True)

    T_O_eV, kappa_O_ion = sol_O

    # Check if solution is valid
    if ier == 1 and kappa_O_ion > 1.505 and T_O_eV > 0:
        Te_O[i] = T_O_eV
        kappa_O[i] = kappa_O_ion
    else:
        print('O+ error', kappa_e,T_e_eV,'set to nans' )
        Te_O[i] = np.nan
        kappa_O[i] = np.nan

    # Now perform the same calculations for S²⁺ ions
    n_s2c_cm3 = ns2p0[i]           # Cold S²⁺ density in cm^-3
    n_s2h_cm3 = ns2phot            # Hot S²⁺ density in cm^-3
    n_S2_cm3 = n_s2c_cm3 + n_s2h_cm3

    n_s2c = n_s2c_cm3 * 1e6        # Convert to m^-3
    n_s2h = n_s2h_cm3 * 1e6
    n_S2 = n_s2c + n_s2h

    T_s2c_eV = ti0[i]              # Cold ion temperature in eV
    T_s2h_eV = toph0[i]            # Hot ion temperature in eV

    # Number density conservation
    n = n_S2

    # Define T_mean and S for S²⁺ ions
    T_mean_S2 = (n_s2c * T_s2c_eV + n_s2h * T_s2h_eV) / n if n > 0 else T_s2c_eV
    S_S2 = (n_s2c * np.sqrt(T_s2c_eV) + n_s2h * np.sqrt(T_s2h_eV)) / n if n > 0 else np.sqrt(T_s2c_eV)

    # Function to solve for S²⁺ ions
    def equations_S2(vars):
        T_e_eV, kappa_e = vars
        f1 = T_e_eV - T_mean_S2 * (1. - 3. / (2. * kappa_e))
        f2 = S_S2 - np.sqrt(T_e_eV * kappa_e) * gamma(kappa_e - 1.) / gamma(kappa_e - 0.5)
        return [f1, f2]

    # Initial guesses for S²⁺ ions
    """
    if i>0:
        initial_guess_S2 = [ Te_S[i-1],kappa_S[i-1]]
    else:
        initial_guess_S2 = [T_mean_S2, 5.0]
    """
    #kappa_guess_s2p_and_Op=4.*((R/6.) ** pwantops2p)
    #initial_guess_S2 = [T_s2c_eV, kappa_guess_s2p_and_Op]
    initial_guess_S2 = [T_mean_S2, 5.0]

    # Solve equations numerically
    if n > 0:
        sol_S2, info, ier, mesg = fsolve(equations_S2, initial_guess_S2, full_output=True)
        T_S2_eV, kappa_S2 = sol_S2

        # Check if solution is valid
        if ier == 1 and kappa_S2 > 1.505 and T_S2_eV > 0:
            Te_S[i] = T_S2_eV
            kappa_S[i] = kappa_S2
        else:
            print('S2+ error', kappa_e,T_e_eV,'set to nans' )
            Te_S[i] = np.nan
            kappa_S[i] = np.nan
    else:
        # If n_S2 is zero, assign NaN
        print('s2+ error', kappa_e,T_e_eV,'set to nans' )
        Te_S[i] = np.nan
        kappa_S[i] = np.nan

print(kappa_electrons[170])
print(kappa_electrons[180])
print(kappa_electrons[190])
for i in range(len(R)):
    if R[i]>5.8:
        Te_electrons[i] = interp1d(R,Te_electrons)(R[i])
        kappa_electrons[i] = interp1d(R,kappa_electrons)(R[i])
        Te_O[i] = interp1d(R,Te_O)(R[i])
        kappa_O[i] = interp1d(R,kappa_O)(R[i])
        Te_S[i] = interp1d(R,Te_S)(R[i])
        kappa_S[i] = interp1d(R,kappa_S)(R[i])
    else:
        if R[i] >5.7:
            Te_electrons[i] = interp1d(R,Te_electrons)(R[i])
            kappa_electrons[i] = interp1d(R,kappa_electrons)(R[i])
            Te_O[i] = interp1d(R,Te_O)(R[i])
            kappa_O[i] = interp1d(R,kappa_O)(R[i])
            Te_S[i] = interp1d(R,Te_S)(R[i])
            kappa_S[i] = interp1d(R,kappa_S)(R[i])
        else:
            if R[i] >5.5:
                
               # vale=interp1d(R,kappa_electrons)(5.8)
                #plvale
                pplvale = np.log(interp1d(R,kappa_electrons)(5.76)/300.)/np.log(5.76/5.5)
                pplvalete = np.log(interp1d(R,Te_electrons)(5.76)/tec0[150])/np.log(5.76/5.5)
                Te_electrons[i] = tec0[150]*((R[i]/5.5)** pplvalete)
                kappa_electrons[i] = 300.*((R[i]/5.5)** pplvale)#300.#vale*((R[i]/5.8) **plvale )
                Te_O[i] = interp1d(R,Te_O)(R[i])
                kappa_O[i] = interp1d(R,kappa_O)(R[i])
                Te_S[i] = interp1d(R,Te_S)(R[i])
                kappa_S[i] = interp1d(R,kappa_S)(R[i])  
            else:
                if R[i] >5.:
                    Te_electrons[i] = tec0[i]
                    kappa_electrons[i] = 300.#interp1d(R,kappa_electrons)(R[i])#*((R[i]/5.8) ** )
                    Te_O[i] = ti0[i]#interp1d(R,Te_O)(R[i])
                    vali = interp1d(R,kappa_O)(5.5)
                    pvali = np.log(300./vali)/np.log(5./5.5) # set to 300. at 5 and vali at 5.5
                    kappa_O[i] =vali*((R[i]/5.5) **pvali ) #interp1d(R,kappa_O)(R[i])
                    Te_S[i] = ti0[i]#interp1d(R,Te_S)(R[i])
                    kappa_S[i] =vali*((R[i]/5.5) **pvali ) #interp1d(R,kappa_S)(R[i])  
                else:
                    Te_electrons[i] =tec0[i]
                    kappa_electrons[i] = 300.
                    Te_O[i] = ti0[i]
                    kappa_O[i] = 300.
                    Te_S[i] = ti0[i]
                    kappa_S[i] = 300.
                    
             
                
Te_electrons_finite = Te_electrons[np.isfinite(Te_electrons)]
R_finite = R[np.isfinite(Te_electrons)]
                
Te_electrons = interp1d(R_finite,Te_electrons_finite)(R)
kappa_electrons_finite = kappa_electrons[np.isfinite(kappa_electrons)]
R_finite = R[np.isfinite(kappa_electrons)]
kappa_electrons = interp1d(R_finite,kappa_electrons_finite)(R)
Te_O = interp1d(R,Te_O)(R)
kappa_O = interp1d(R,kappa_O)(R)
Te_S = interp1d(R,Te_S)(R)
kappa_S = interp1d(R,kappa_S)(R)


            

netot_nom = nec0 + neh0

nec_nom = (1. - feh_nom )* netot_nom
neh_nom = (feh_nom )* netot_nom

np.savetxt('new_nominal_model_Teh_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv', teh_nom, delimiter=',')
np.savetxt('new_nominal_model_nec_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv', nec_nom, delimiter=',')
np.savetxt('new_nominal_model_neh_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv', neh_nom, delimiter=',')


# Save the computed parameters to CSV files
np.savetxt('new_nominal_model_Te_electrons_012_momentsequal_bimax_to_kappa.csv', Te_electrons, delimiter=',')
np.savetxt('new_nominal_model_kappa_electrons_012_momentsequal_bimax_to_kappa.csv', kappa_electrons, delimiter=',')
np.savetxt('Te_O_012_momentsequal_bimax_to_kappa.csv', Te_O, delimiter=',')
np.savetxt('kappa_O_012_momentsequal_bimax_to_kappa.csv', kappa_O, delimiter=',')
np.savetxt('Te_S_012_momentsequal_bimax_to_kappa.csv', Te_S, delimiter=',')
np.savetxt('kappa_S_012_momentsequal_bimax_to_kappa.csv', kappa_S, delimiter=',')

R = np.linspace(4., 10., 601)
# Set default parameters for bold text and thicker tick marks


# Set default parameters for bold text and thicker tick marks
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.titleweight'] = 'bold'
plt.rcParams['axes.linewidth'] = 2
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['ytick.major.width'] = 2

# First plot
plt.figure(figsize=(10, 6))
plt.semilogy(R, tec0, label='Electron Cold Temperature', linestyle='--')
plt.semilogy(R, Te_electrons, label='Electron Kappa Temperature')
plt.semilogy(R, ti0, label='O⁺ & S²⁺ Cold Temperature', linestyle='--')
plt.semilogy(R, Te_O, label='O⁺ Kappa Temperature')
plt.semilogy(R, Te_S, label='S²⁺ Kappa Temperature')
plt.xlabel(r'$\rho_c$ ($R_J$)', fontsize=14, fontweight='bold')
plt.ylabel('Temperature (eV)', fontsize=14, fontweight='bold')
plt.title('Nominal Standard Isotropic Kappa Temperatures', fontsize=16, fontweight='bold')
plt.grid(True)
plt.tight_layout()

# Get current axes
ax = plt.gca()
# Set tick parameters
ax.tick_params(axis='both', which='major', labelsize=12, width=2)
# Set tick label font weight
for label in ax.get_xticklabels():
    label.set_fontweight('bold')
for label in ax.get_yticklabels():
    label.set_fontweight('bold')

# Set legend with bold font
plt.legend(fontsize=12, prop={'weight': 'bold'})

plt.savefig('new_nominal_model_kappa_temperatures_vs_radial_distance_012_momentsequal_bimax_to_kappa.png', dpi=300)
plt.show()

# Load data for the second plot
r_vs_kappa_steffl2004b = np.loadtxt('r_vs_kappa_steffl2004b.csv', delimiter=',', skiprows=1)
r_vs_kappa_steffl2004b_r = r_vs_kappa_steffl2004b[:, 0]
r_vs_kappa_steffl2004b_kappa = r_vs_kappa_steffl2004b[:, 1]

# Second plot
plt.figure(figsize=(10, 6))
plt.semilogy(R, kappa_electrons, label='Electron κ')
plt.semilogy(R, kappa_O, label='O⁺ κ')
plt.semilogy(R, kappa_S, label='S²⁺ κ')
plt.semilogy(r_vs_kappa_steffl2004b_r, r_vs_kappa_steffl2004b_kappa, label='Steffl (2004b) Electron κ')
plt.xlabel(r'$\rho_c$ ($R_J$)', fontsize=14, fontweight='bold')
plt.ylabel('Kappa Value (κ)', fontsize=14, fontweight='bold')
plt.title('Nominal Standard Isotropic Kappa Values', fontsize=16, fontweight='bold')
plt.grid(True)
plt.tight_layout()

# Get current axes
ax = plt.gca()
# Set tick parameters
ax.tick_params(axis='both', which='major', labelsize=12, width=2)
# Set tick label font weight
for label in ax.get_xticklabels():
    label.set_fontweight('bold')
for label in ax.get_yticklabels():
    label.set_fontweight('bold')

# Set legend with bold font
plt.legend(fontsize=12, prop={'weight': 'bold'})

plt.savefig('new_nominal_model_kappa_values_vs_radial_distance_012_momentsequal_bimax_to_kappa.png', dpi=300)
plt.show()


print('Computation completed and results saved to CSV files.')
