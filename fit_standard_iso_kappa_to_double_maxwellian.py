import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gamma

# Physical constants
m_e = 9.10938356e-31        # Electron mass in kg
e = 1.602176634e-19         # Elementary charge in C (J/eV)
pi = np.pi



# Given parameters (feel free to change these)
n_ec_cm3 = 30.5             # Core electron density in cm^-3
n_eh_cm3 = 1.52             # Hot electron density in cm^-3
T_ec_eV = 17.79             # Core electron temperature in eV
T_eh_eV = 299.11            # Hot electron temperature in eV
T_e_eV = 17.79#100.0              # Kappa distribution temperature in eV
kappa_e = 3.#2.0               # Kappa parameter

necfrac = round(n_ec_cm3/(n_eh_cm3 + n_ec_cm3), 3)

nehfrac = round(n_eh_cm3/(n_eh_cm3 + n_ec_cm3),3)

# Convert densities to m^-3
n_ec = n_ec_cm3 * 1e6       # Core electron density in m^-3
n_eh = n_eh_cm3 * 1e6       # Hot electron density in m^-3
n_e = n_ec + n_eh           # Total electron density in m^-3

# Convert temperatures to Joules
T_ec = T_ec_eV * e          # Core temperature in J
T_eh = T_eh_eV * e          # Hot temperature in J
T_e = T_e_eV * e            # Kappa temperature in J

# Velocity grid
v_max = 1.5 * np.sqrt(2 * T_eh / m_e)   # Maximum velocity for the grid
v = np.linspace(0, v_max, 1000)       # Linearly spaced velocity array

# To avoid plotting log(0) for f_e(v=0), set the first velocity to a small value
#v[0] = v[1]*1e-3  # Set v[0] to a small non-zero value

# Constants for distributions
A_ec = n_ec * (m_e / (2 * pi * T_ec))**(1.5)
A_eh = n_eh * (m_e / (2 * pi * T_eh))**(1.5)
A_kappa = n_e * (m_e / (2 * pi * T_e))**(1.5) * gamma(kappa_e) / ( np.sqrt(kappa_e) * gamma(kappa_e - 0.5) )

# Double Maxwellian distribution
f_double = A_ec * np.exp(- m_e * v**2 / (2 * T_ec)) + A_eh * np.exp(- m_e * v**2 / (2 * T_eh))

# Kappa distribution
f_kappa = A_kappa * (1 + m_e * v**2 / (2 * kappa_e * T_e))**( - (kappa_e + 1))

# Plotting
plt.figure(figsize=(10, 6))
plt.semilogy((0.5 * m_e * v ** 2.) / e, f_double, label='Double Maxwellian ($n_{ec}/n_{e}$=' + f'{necfrac},' + ' $n_{eh}/n_{e}$=' + f'{nehfrac},' + ' $T_{ec}$=' + f'{T_ec_eV} eV,' + ' $T_{eh}$=' + f'{T_eh_eV} eV)', linewidth=2)
plt.semilogy((0.5 * m_e * v ** 2.) / e, f_kappa, label=f'Kappa Distribution (κₑ={kappa_e}, Tₑ={T_e_eV} eV)', linestyle='--', linewidth=2)
plt.xlabel('Electron KE $1/2 m v^2$ (eV)', fontsize=14)
plt.ylabel('Distribution Function $f_e$', fontsize=14)
plt.title('Comparison of Double Maxwellian and Kappa Distribution', fontsize=16)
plt.legend(fontsize=12)
plt.grid(True, which="both", ls="--")
plt.tight_layout()
plt.show()
