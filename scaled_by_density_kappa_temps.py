#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  9 03:39:15 2024

@author: edne8319
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  9 02:55:04 2024

@author: edne8319
"""


import numpy as np
import matplotlib.pyplot as plt


# Load the radial profiles from CSV files
nop = np.loadtxt('nop_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv', delimiter=',')
ns2p = np.loadtxt('ns2p_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv', delimiter=',')
noph = np.loadtxt('noph_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv', delimiter=',')


ns2ph = 0.5*noph
noph = 0.5*noph 



nc = nop

nh = noph 





tc = np.loadtxt('Tic_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv', delimiter=',')
th = np.loadtxt('Toph_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv', delimiter=',')
kappa_values = np.loadtxt('kappa_vals_current.csv', delimiter=',')


T_fit =( nc *tc + nh*th) / (nc + nh) #np.array(T_fit)
np.savetxt('weighted_temp_temperatures_op.csv', T_fit, delimiter=',')

# Step 5: Generate and Save the Plots

# Generate the rho array
rho = np.linspace(4, 10, 601)

# Plot tc, th, and T_fit vs rho
plt.figure(figsize=(10, 6))
plt.plot(rho, tc, label='T_c', linestyle='-', marker='o')
plt.plot(rho, th, label='T_h', linestyle='-', marker='x')
plt.plot(rho, T_fit, label='T_fit', linestyle='-', marker='s')

# Adding labels and legend
plt.xlabel('rho')
plt.ylabel('Temperature')
plt.title('Temperature Profiles vs. rho')
plt.legend()
plt.grid(True)

# Save the plot to a file
plt.savefig('weighted_temperature_profiles_vs_rho.png')

# Show the plot
plt.show()

# Plot nc and nh vs rho
plt.figure(figsize=(10, 6))
plt.plot(rho, nc, label='n_c', linestyle='-', marker='o')
plt.plot(rho, nh, label='n_h', linestyle='-', marker='x')

# Adding labels and legend
plt.xlabel('rho')
plt.ylabel('Density')
plt.title('Density Profiles vs. rho')
plt.legend()
plt.grid(True)

# Save the plot to a file


# Save the plot to a file
plt.savefig('density_profiles_vs_rho.png')

# Show the plot
plt.show()

