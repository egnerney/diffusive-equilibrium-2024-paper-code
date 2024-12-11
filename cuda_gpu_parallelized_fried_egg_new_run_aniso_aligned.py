import numpy as np
from numba import cuda, float64, int32
import math

# Define constants
RP_JUPITER = 7.1492e7       # Radius of Jupiter in meters
OMEGA_JUPITER = 1.7583586898863878e-4  # Rotation rate in rad/s
GM_JUPITER = 1.2668689059e17  # Gravitational constant * mass (m^3/s^2)

# Device function for linear interpolation
@cuda.jit(device=True)
def linear_interpolate(x, xp, fp):
    idx = 0
    for i in range(len(xp) - 1):
        if xp[i] <= x <= xp[i + 1]:
            idx = i
            break
    x0 = xp[idx]
    x1 = xp[idx + 1]
    y0 = fp[idx]
    y1 = fp[idx + 1]
    slope = (y1 - y0) / (x1 - x0)
    return y0 + slope * (x - x0)

# Device function for net charge density calculation
@cuda.jit(device=True)
def net_charge_density(Phi, deltaU, Bratio, species_m, species_q, species_T, species_A, species_kappa, species_n, species_type):
    nq = 0.0
    for i in range(len(species_m)):
        m = species_m[i]
        q = species_q[i]
        T = species_T[i]
        A = species_A[i]
        kappa = species_kappa[i]
        n0 = species_n[i]
        dist_type = species_type[i]
        exponent = -(m * deltaU + q * Phi) / T

        # Maxwellian density
        n = n0 * math.exp(exponent)

        nq += q * n
    return nq

# Kernel function to process field lines
@cuda.jit
def process_field_line_kernel(
    x_fl, y_fl, z_fl, rho_fl, lat_fl, phi_fl_wlong, r_fl, Btot_fl, r_00,
    nop0, no2p0, nsp0, ns2p0, ns3p0, nhp0, nnap0, noph0, nec0, neh0,
    ti0, thp0, tec0, kappa_Te, kappa_Top, kappa_e, kappa_op, kappa_Thp,
    species_m, species_q, species_T, species_A, species_kappa, species_n, species_type,
    n_out_op, n_out_hp, n_out_elec, x_out, y_out, z_out, B_out
):
    j = cuda.grid(1)
    if j < x_fl.shape[1]:
        # Planet-specific variables
        RP = RP_JUPITER
        Omega = OMEGA_JUPITER
        GM = GM_JUPITER

        # Interpolate field line data
        npointsfieldline = lat_fl.shape[0]
        for i in range(npointsfieldline):
            lat = lat_fl[i, j]
            x_out[j, i] = x_fl[i, j]
            y_out[j, i] = y_fl[i, j]
            z_out[j, i] = z_fl[i, j]
            B_out[j, i] = Btot_fl[i, j]

            # Compute deltaU
            r = r_fl[i, j] * RP
            U = -GM / r - 0.5 * r ** 2 * math.cos(math.radians(lat)) ** 2 * Omega ** 2
            U0 = U  # Assuming U0 is U at equator
            deltaU = U - U0

            # Compute Bratio
            B_eq = Btot_fl[npointsfieldline // 2, j]
            Bratio = Btot_fl[i, j] / B_eq

            # Solve for Phi using bisection method
            Phi_low = -1e3
            Phi_high = 1e3
            tol = 1e-5
            max_iter = 100
            for iter in range(max_iter):
                Phi_mid = 0.5 * (Phi_low + Phi_high)
                nq_mid = net_charge_density(Phi_mid, deltaU, Bratio, species_m, species_q, species_T, species_A, species_kappa, species_n, species_type)
                if nq_mid > 0:
                    Phi_high = Phi_mid
                else:
                    Phi_low = Phi_mid
                if abs(Phi_high - Phi_low) < tol:
                    break
            Phi = Phi_mid

            # Compute densities
            n_op = species_n[0] * math.exp(-(species_m[0] * deltaU + species_q[0] * Phi) / species_T[0])
            n_hp = species_n[5] * math.exp(-(species_m[5] * deltaU + species_q[5] * Phi) / species_T[5])
            n_elec = species_n[9] * math.exp(-(species_m[9] * deltaU + species_q[9] * Phi) / species_T[9])

            n_out_op[j, i] = n_op
            n_out_hp[j, i] = n_hp
            n_out_elec[j, i] = n_elec

# Main function
def main():
    # Load data files
    x_fl = np.transpose(np.loadtxt('int_aligned_x.csv', delimiter=','))
    y_fl = np.transpose(np.loadtxt('int_aligned_y.csv', delimiter=','))
    z_fl = np.transpose(np.loadtxt('int_aligned_z.csv', delimiter=','))
    rho_fl = np.sqrt(x_fl ** 2 + y_fl ** 2)
    r_fl = np.sqrt(x_fl ** 2 + y_fl ** 2 + z_fl ** 2)
    lat_fl = np.degrees(np.arcsin(z_fl / r_fl))
    phi_fl = np.degrees(np.arctan2(y_fl, x_fl))
    phi_fl_wlong = 360.0 - phi_fl

    Btot_fl = np.transpose(np.loadtxt('Btot_int_aligned.csv', delimiter=','))

    # Species parameters
    species_m = np.array([16.0, 1.0, 1.0 / 1837.0])  # O+, H+, e-
    species_q = np.array([1.0, 1.0, -1.0])
    species_T = np.array([79.3, 79.3, 4.6])
    species_A = np.array([1.0, 1.0, 1.0])
    species_kappa = np.array([100.0, 100.0, 100.0])
    species_n = np.array([592.0, 50.6, 2537.5])
    species_type = np.array([0, 0, 0])  # 0 for Maxwellian

    # Pre-allocate output arrays
    nrfls = x_fl.shape[1]
    npointsfieldline = x_fl.shape[0]
    n_out_op = np.zeros((nrfls, npointsfieldline))
    n_out_hp = np.zeros((nrfls, npointsfieldline))
    n_out_elec = np.zeros((nrfls, npointsfieldline))
    x_out = np.zeros((nrfls, npointsfieldline))
    y_out = np.zeros((nrfls, npointsfieldline))
    z_out = np.zeros((nrfls, npointsfieldline))
    B_out = np.zeros((nrfls, npointsfieldline))

    # Transfer data to GPU
    d_x_fl = cuda.to_device(x_fl)
    d_y_fl = cuda.to_device(y_fl)
    d_z_fl = cuda.to_device(z_fl)
    d_rho_fl = cuda.to_device(rho_fl)
    d_lat_fl = cuda.to_device(lat_fl)
    d_phi_fl_wlong = cuda.to_device(phi_fl_wlong)
    d_r_fl = cuda.to_device(r_fl)
    d_Btot_fl = cuda.to_device(Btot_fl)

    d_species_m = cuda.to_device(species_m)
    d_species_q = cuda.to_device(species_q)
    d_species_T = cuda.to_device(species_T)
    d_species_A = cuda.to_device(species_A)
    d_species_kappa = cuda.to_device(species_kappa)
    d_species_n = cuda.to_device(species_n)
    d_species_type = cuda.to_device(species_type)

    d_n_out_op = cuda.to_device(n_out_op)
    d_n_out_hp = cuda.to_device(n_out_hp)
    d_n_out_elec = cuda.to_device(n_out_elec)
    d_x_out = cuda.to_device(x_out)
    d_y_out = cuda.to_device(y_out)
    d_z_out = cuda.to_device(z_out)
    d_B_out = cuda.to_device(B_out)

    # Launch kernel
    threadsperblock = 64
    blockspergrid = (nrfls + (threadsperblock - 1)) // threadsperblock

    process_field_line_kernel[blockspergrid, threadsperblock](
        d_x_fl, d_y_fl, d_z_fl, d_rho_fl, d_lat_fl, d_phi_fl_wlong, d_r_fl, d_Btot_fl, None,
        None, None, None, None, None, None, None, None, None, None, None,
        None, None, None, None, None, None, None,
        d_species_m, d_species_q, d_species_T, d_species_A, d_species_kappa, d_species_n, d_species_type,
        d_n_out_op, d_n_out_hp, d_n_out_elec, d_x_out, d_y_out, d_z_out, d_B_out
    )

    # Copy results back to host
    n_out_op = d_n_out_op.copy_to_host()
    n_out_hp = d_n_out_hp.copy_to_host()
    n_out_elec = d_n_out_elec.copy_to_host()
    x_out = d_x_out.copy_to_host()
    y_out = d_y_out.copy_to_host()
    z_out = d_z_out.copy_to_host()
    B_out = d_B_out.copy_to_host()

    # Save results
    np.savez_compressed('n_out_op.npz', n_out_op=n_out_op)
    np.savez_compressed('n_out_hp.npz', n_out_hp=n_out_hp)
    np.savez_compressed('n_out_elec.npz', n_out_elec=n_out_elec)

if __name__ == '__main__':
    main()
