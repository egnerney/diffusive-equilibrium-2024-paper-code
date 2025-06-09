# -*- coding: utf-8 -*-
"""
Created on Sat Apr 19 16:56:47 2025

@author: Owner
"""

import numpy as np


r_eq_vals_outer = 0.01 * np.arange(70) + 4.0

# Each of shape (70,1000)



x_fl_inner_ = np.loadtxt("x_io_aligned_fls_jrm33+con20_int_70_field_lines_0-69_4-4.69RJ.csv", delimiter=",")



y_fl_inner_ = np.loadtxt("y_io_aligned_fsl_jrm33+con20_int_70_field_lines_0-69_4-4.69RJ.csv",  delimiter=",")

z_fl_inner_ = np.loadtxt("z_io_aligned_fls_jrm33+con20_int_70_field_lines_0-69_4-4.69RJ.csv",  delimiter=",")


Bx_fl_inner_ = np.loadtxt("Bx_io_aligned_fls_jrm33+con20_int_70_field_lines_0-69_4-4.69RJ.csv",  delimiter=",")

By_fl_inner_ = np.loadtxt("By_io_aligned_fls_jrm33+con20_int_70_field_lines_0-69_4-4.69RJ.csv", delimiter=",")

Bz_fl_inner_ = np.loadtxt("Bz_io_aligned_fls_jrm33+con20_int_70_field_lines_0-69_4-4.69RJ.csv",  delimiter=",")

Btot_fl_inner_ = np.loadtxt("B_io_aligned_fls_jrm33+con20_int_70_field_lines_0-69_4-4.69RJ.csv",  delimiter=",")


x_fl_inner = np.transpose(x_fl_inner_,(1,0))
y_fl_inner = np.transpose(y_fl_inner_,(1,0))
z_fl_inner = np.transpose(z_fl_inner_,(1,0))


Bx_fl_inner = np.transpose(Bx_fl_inner_,(1,0))
By_fl_inner = np.transpose(By_fl_inner_,(1,0))
Bz_fl_inner = np.transpose(Bz_fl_inner_,(1,0))


Btot_fl_inner = np.transpose(Btot_fl_inner_,(1,0))



r_eq_vals_outer = 0.01 * np.arange(29) + 9.72
# Each of shape (29,1000)
x_fl_outer_ = np.loadtxt("x_io_aligned_fls_jrm33+con20_int_29_field_lines_572-600_9.72-10.00RJ_phi0=75.85-0.52=75.33.csv", delimiter=",")

y_fl_outer_ = np.loadtxt("y_io_aligned_fsl_jrm33+con20_int_29_field_lines_572-600_9.72-10.00RJ_phi0=75.85-0.52=75.33.csv",  delimiter=",")

z_fl_outer_ = np.loadtxt("z_io_aligned_fls_jrm33+con20_int_29_field_lines_572-600_9.72-10.00RJ_phi0=75.85-0.52=75.33.csv",  delimiter=",")


Bx_fl_outer_ = np.loadtxt("Bx_io_aligned_fls_jrm33+con20_int_29_field_lines_572-600_9.72-10.00RJ_phi0=75.85-0.52=75.33.csv",  delimiter=",")

By_fl_outer_ = np.loadtxt("By_io_aligned_fls_jrm33+con20_int_29_field_lines_572-600_9.72-10.00RJ_phi0=75.85-0.52=75.33.csv", delimiter=",")

Bz_fl_outer_ = np.loadtxt("Bz_io_aligned_fls_jrm33+con20_int_29_field_lines_572-600_9.72-10.00RJ_phi0=75.85-0.52=75.33.csv",  delimiter=",")

Btot_fl_outer_ = np.loadtxt("B_io_aligned_fls_jrm33+con20_int_29_field_lines_572-600_9.72-10.00RJ_phi0=75.85-0.52=75.33.csv",  delimiter=",")


x_fl_outer = np.transpose(x_fl_outer_,(1,0))
y_fl_outer = np.transpose(y_fl_outer_,(1,0))
z_fl_outer = np.transpose(z_fl_outer_,(1,0))


Bx_fl_outer = np.transpose(Bx_fl_outer_,(1,0))
By_fl_outer = np.transpose(By_fl_outer_,(1,0))
Bz_fl_outer = np.transpose(Bz_fl_outer_,(1,0))


Btot_fl_outer = np.transpose(Btot_fl_outer_,(1,0))



r_eq_vals = 0.01 * np.arange(601) + 4.0
# Each of shape (601,1000)
x_fl = np.transpose(np.loadtxt('int_aligned_x_jrm33+con20integral_601_4-10_properly_aligned_phi0=75.85_perfectaligned_for_4.70-9.71.csv', delimiter=','))
y_fl = np.transpose(np.loadtxt('int_aligned_y_jrm33+con20integral_601_4-10_properly_aligned_phi0=75.85_perfectaligned_for_4.70-9.71.csv', delimiter=','))
z_fl = np.transpose(np.loadtxt('int_aligned_z_jrm33+con20integral_601_4-10_properly_aligned_phi0=75.85_perfectaligned_for_4.70-9.71.csv', delimiter=','))


Bx_fl = np.transpose(np.loadtxt('Bx_int_aligned_jrm33+con20integral_601_4-10_properly_aligned_phi0=75.85_perfectaligned_for_4.70-9.71.csv', delimiter=','))

By_fl = np.transpose(np.loadtxt('By_int_aligned_jrm33+con20integral_601_4-10_properly_aligned_phi0=75.85_perfectaligned_for_4.70-9.71.csv', delimiter=','))

Bz_fl = np.transpose(np.loadtxt('Bz_int_aligned_jrm33+con20integral_601_4-10_properly_aligned_phi0=75.85_perfectaligned_for_4.70-9.71.csv', delimiter=','))

Btot_fl = np.transpose(np.loadtxt('Btot_int_aligned_jrm33+con20integral_601_4-10_properly_aligned_phi0=75.85_perfectaligned_for_4.70-9.71.csv', delimiter=','))



x_new_fl = np.empty_like(x_fl)           # keeps dtype, avoids filling with 1’s
y_new_fl = np.empty_like(y_fl) 
z_new_fl = np.empty_like(z_fl) 
Bx_new_fl = np.empty_like(Bx_fl) 
By_new_fl = np.empty_like(By_fl) 
Bz_new_fl = np.empty_like(Bz_fl) 
Btot_new_fl = np.empty_like(Btot_fl) 

x_new_fl[:, 0:70]   = x_fl_inner         # columns 0 … 69   (70 columns)
x_new_fl[:, 70:572] = x_fl[:, 70:572]    # columns 70 … 571 (502 columns)
x_new_fl[:, 572:]   = x_fl_outer         # columns 572 … 600 (29 columns)


y_new_fl[:, 0:70]   = y_fl_inner         # columns 0 … 69   (70 columns)
y_new_fl[:, 70:572] = y_fl[:, 70:572]    # columns 70 … 571 (502 columns)
y_new_fl[:, 572:]   = y_fl_outer         # columns 572 … 600 (29 columns)


z_new_fl[:, 0:70]   = z_fl_inner         # columns 0 … 69   (70 columns)
z_new_fl[:, 70:572] = z_fl[:, 70:572]    # columns 70 … 571 (502 columns)
z_new_fl[:, 572:]   = z_fl_outer         # columns 572 … 600 (29 columns)



Bx_new_fl[:, 0:70]   = Bx_fl_inner         # columns 0 … 69   (70 columns)
Bx_new_fl[:, 70:572] = Bx_fl[:, 70:572]    # columns 70 … 571 (502 columns)
Bx_new_fl[:, 572:]   = Bx_fl_outer         # columns 572 … 600 (29 columns)


By_new_fl[:, 0:70]   = By_fl_inner         # columns 0 … 69   (70 columns)
By_new_fl[:, 70:572] = By_fl[:, 70:572]    # columns 70 … 571 (502 columns)
By_new_fl[:, 572:]   = By_fl_outer         # columns 572 … 600 (29 columns)


Bz_new_fl[:, 0:70]   = Bz_fl_inner         # columns 0 … 69   (70 columns)
Bz_new_fl[:, 70:572] = Bz_fl[:, 70:572]    # columns 70 … 571 (502 columns)
Bz_new_fl[:, 572:]   = Bz_fl_outer         # columns 572 … 600 (29 columns)


Btot_new_fl[:, 0:70]   = Btot_fl_inner         # columns 0 … 69   (70 columns)
Btot_new_fl[:, 70:572] = Btot_fl[:, 70:572]    # columns 70 … 571 (502 columns)
Btot_new_fl[:, 572:]   = Btot_fl_outer         # columns 572 … 600 (29 columns)

np.savetxt('1000x601_shape_int_x_jrm33+con20integral_601_4-10_properly_aligned_whole_range_by_shifting_ends+-0.52_from_phi0=75.85_which_is_only_goodfor_4.70-9.71_without_shifts_at_ends.csv',x_new_fl, delimiter=',')

np.savetxt('1000x601_shape_int_y_jrm33+con20integral_601_4-10_properly_aligned_whole_range_by_shifting_ends+-0.52_from_phi0=75.85_which_is_only_goodfor_4.70-9.71_without_shifts_at_ends.csv',y_new_fl, delimiter=',')

np.savetxt('1000x601_shape_int_z_jrm33+con20integral_601_4-10_properly_aligned_whole_range_by_shifting_ends+-0.52_from_phi0=75.85_which_is_only_goodfor_4.70-9.71_without_shifts_at_ends.csv',z_new_fl, delimiter=',')

np.savetxt('1000x601_shape_int_Bx_jrm33+con20integral_601_4-10_properly_aligned_whole_range_by_shifting_ends+-0.52_from_phi0=75.85_which_is_only_goodfor_4.70-9.71_without_shifts_at_ends.csv',Bx_new_fl, delimiter=',')

np.savetxt('1000x601_shape_int_By_jrm33+con20integral_601_4-10_properly_aligned_whole_range_by_shifting_ends+-0.52_from_phi0=75.85_which_is_only_goodfor_4.70-9.71_without_shifts_at_ends.csv',By_new_fl, delimiter=',')

np.savetxt('1000x601_shape_int_Bz_jrm33+con20integral_601_4-10_properly_aligned_whole_range_by_shifting_ends+-0.52_from_phi0=75.85_which_is_only_goodfor_4.70-9.71_without_shifts_at_ends.csv',Bz_new_fl, delimiter=',')

np.savetxt('1000x601_shape_int_Btot_jrm33+con20integral_601_4-10_properly_aligned_whole_range_by_shifting_ends+-0.52_from_phi0=75.85_which_is_only_goodfor_4.70-9.71_without_shifts_at_ends.csv',Btot_new_fl, delimiter=',')

