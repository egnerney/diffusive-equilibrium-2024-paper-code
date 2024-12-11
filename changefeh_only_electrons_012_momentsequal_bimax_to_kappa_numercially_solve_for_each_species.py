import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gamma
from scipy.optimize import fsolve

# Load data files (assuming they are in the same directory)
nec0 = np.loadtxt('nec_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv', delimiter=',')
neh0 = np.loadtxt('neh_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv', delimiter=',')
tec0 = np.loadtxt('Tec_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv', delimiter=',')
teh0 = np.loadtxt('Teh_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv', delimiter=',')

# Ensure all data arrays have the same length
num_points = len(tec0)

# Physical constants
m_e = 9.10938356e-31        # Electron mass in kg
e = 1.602176634e-19         # Elementary charge in C (J/eV)
pi = np.pi

# Radial positions from 4 to 10 RJ
R = np.linspace(4, 10, num_points)  # Radial distance array

# Load Steffl data
r_vs_kappa_steffl2004b = np.loadtxt('r_vs_kappa_steffl2004b.csv', delimiter=',', skiprows=1)
r_vs_kappa_steffl2004b_r = r_vs_kappa_steffl2004b[:, 0]
r_vs_kappa_steffl2004b_kappa = r_vs_kappa_steffl2004b[:, 1]

# Define feh0 and fehpl lists
feh0_list = [0.002]#[0.001, 0.002, 0.003]
fehpl_list = [5.,6.,7.,8.,9.]#[4., 5., 6.]

# Initialize figures before the loops
# Initialize figure for temperatures
fig_temp, ax_temp = plt.subplots(figsize=(10, 6))
ax_temp.semilogy(R, tec0, label='Electron Cold Temperature', linestyle='--')

# Initialize figure for kappa values
fig_kappa, ax_kappa = plt.subplots(figsize=(10, 6))
ax_kappa.semilogy(r_vs_kappa_steffl2004b_r, r_vs_kappa_steffl2004b_kappa,label='Steffl (2004b) Electron κ', linestyle='--', linewidth=2)


stringgy = r'T$_{eh}$ = 270$\left(\frac{\rho_c}{6}\right)^{1.56}$ eV'
#stringgy = r'T$_{eh}$ = Constant 270 eV'
# Loop over feh0 and fehpl
for feh0 in feh0_list:
    for fehpl in fehpl_list:
        # Initialize arrays to store computed parameters
        Te_electrons = np.zeros(num_points)
        kappa_electrons = np.zeros(num_points)
        

        
        # Loop over all radial positions
        for i in range(num_points):
            # Electrons
            #Teh_ev = 270. * ((R[i]/6.) ** 1.20626)#500.#270. #35. * ((R[i]/6.) ** 4.2) #teh0[i]#500.
            Teh_ev = 270. * ((R[i]/6.) ** 1.5632)#500.#270. #35. * ((R[i]/6.) ** 4.2) #teh0[i]#500.
            T_eh_eV = Teh_ev
            n_ec_cm3 = nec0[i]             # Core electron density in cm^-3
            n_eh_cm3 = neh0[i]             # Hot electron density in cm^-3
            T_ec_eV = tec0[i]              # Core electron temperature in eV

            # Convert densities to m^-3
            n_ec = n_ec_cm3 * 1e6       # Core electron density in m^-3
            n_eh = n_eh_cm3 * 1e6       # Hot electron density in m^-3
            n_e = n_ec + n_eh           # Total electron density in m^-3

            # Recalculate n_eh and n_ec using feh0 and fehpl
            if R[i] >= 5.7:
                feh = feh0 * ((R[i]/6.) ** fehpl)
            else:
                feh = 0.0

            n_eh = n_e * feh
            n_ec = n_e * (1. - feh)

            # Number density conservation
            n = n_e

            # Define T_mean and S for electrons
            T_mean = (n_ec * T_ec_eV + n_eh * T_eh_eV) / n
            S = (n_ec * np.sqrt(T_ec_eV) + n_eh * np.sqrt(T_eh_eV)) / n  # Corrected variable name

            # Function to solve for electrons
            def equations_electron(vars):
                T_e_eV, kappa_e = vars
                f1 = T_e_eV - T_mean * (1. - 3. / (2. * kappa_e))
                f2 = S - np.sqrt(T_e_eV * kappa_e) * gamma(kappa_e - 1.) / gamma(kappa_e - 0.5)
                return [f1, f2]

            # Initial guesses for electrons
            if R[i]>7:
               initial_guess = [T_mean, 5.0]
            else:
                initial_guess = [T_ec_eV, 4.]

            # Solve equations numerically
            sol_electron, info, ier, mesg = fsolve(equations_electron, initial_guess, full_output=True)

            T_e_eV, kappa_e = sol_electron

            # Check if solution is valid
            if ier == 1 and kappa_e > 1.505 and T_e_eV > 0:
                Te_electrons[i] = T_e_eV
                kappa_electrons[i] = kappa_e
            else:
                print(f'Electron error at R={R[i]:.2f}, kappa_e={kappa_e}, T_e_eV={T_e_eV}, setting to NaN')
                print('ier=', ier, 'message:', mesg)
                Te_electrons[i] = np.nan
                kappa_electrons[i] = np.nan

        # Save the computed parameters to CSV files
        #np.savetxt(stringgy + f'Teh={Teh_ev:.0f},feh0={feh0:.3f},fehpl={fehpl:.1f}_Te_electrons.csv', Te_electrons, delimiter=',')
        #np.savetxt(stringgy + f'Teh={Teh_ev:.0f},feh0={feh0:.3f},fehpl={fehpl:.1f}_kappa_electrons.csv', kappa_electrons, delimiter=',')

        # Plot temperatures
        ax_temp.semilogy(R, Te_electrons, label=f'feh0={feh0}, fehpl={fehpl}', linewidth=1)

        # Plot kappa values
        ax_kappa.semilogy(R, kappa_electrons, label=f'feh0={feh0}, fehpl={fehpl}', linewidth=1)



# Finalize temperature plot
ax_temp.set_xlabel(r'$\rho_c$ ($R_J$)', fontsize=14)
ax_temp.set_ylabel('Temperature (eV)', fontsize=14)
ax_temp.set_title('Electron Kappa Temperatures vs Radial Distance, ' + stringgy, fontsize=16)
ax_temp.legend(fontsize=12)
ax_temp.grid(True)
fig_temp.tight_layout()
fig_temp.savefig('kappa_temperatures_vs_radial_distance.png', dpi=300)
# Do not call plt.show() here

# Finalize kappa plot
ax_kappa.set_xlabel(r'$\rho_c$ ($R_J$)', fontsize=14)
ax_kappa.set_ylabel('Kappa Value (κ)', fontsize=14)
ax_kappa.set_title('Electron Kappa Values vs Radial Distance, ' + stringgy, fontsize=16)
ax_kappa.legend(fontsize=12)
ax_kappa.grid(True)
fig_kappa.tight_layout()
fig_kappa.savefig('kappa_values_vs_radial_distance.png', dpi=300)
# Do not call plt.show() here

# Call plt.show() once at the end
plt.show()

print('Computation completed and results saved to CSV files.')
