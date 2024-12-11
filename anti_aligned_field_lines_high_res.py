
import JupiterMag as jm

import numpy as np

#import csv

#import matplotlib.pyplot as plt

#import DateTimeTools

#import RecarrayTools

#import PyFileIO

import time
start_time = time.time()


#be sure to configure external field model prior to tracing
jm.Con2020.Config(equation_type='integral')
#this may also become necessary with internal models in future, e.g.
#setting the model degree

#create some starting positions
n = 701

#m = 36
theta = np.pi/2.0
tiltedephi = 165.12#164.62#159.2

#phi = (tiltedephi - 90.0)*np.pi/180.0 # aligned Elong phi eq for jrm33 #(155.8 + 90.0)*np.pi/180.0 #np.linspace(0,359,m)*np.pi/180.0
phi = tiltedephi *np.pi/180.0 # aligned Elong phi eq for jrm33 #(155.8 + 90.0)*np.pi/180.0 #np.linspace(0,359,m)*np.pi/180.0

r = np.linspace(3.5,10.5,n)
x0 = r*np.cos(phi)*np.sin(theta)
y0 = r*np.sin(phi)*np.sin(theta)
z0 = r*np.cos(theta)

#create trace objects, pass starting position(s) x0,y0,z0
#T0 = jm.TraceField(x0,y0,z0,Verbose=True,IntModel='jrm33',ExtModel='none')
T0 = jm.TraceField(x0,y0,z0,Verbose=True,IntModel='jrm33',ExtModel='Con2020')

#plot a trace
ax = T0.PlotRhoZ(label='JRM33 + Con2020 Anti Aligned',color='black')
#ax = T1.PlotRhoZ(fig=ax,label='JRM33 + Con2020',color='red')

ax.set_xlim(-2.0,15.0)
ax.set_ylim(-6.0,6.0)


#print(T0.x[0,:][np.logical_not(np.isnan(T0.x[0,:]))])

#xout=np.empty([])

#i=0
#while i <= n:
#    print(T0.x[i,:][np.logical_not(np.isnan(T0.x[i,:]))].shape)
#    print(T0.y[i,:][np.logical_not(np.isnan(T0.y[i,:]))].shape)
#    print(T0.z[i,:][np.logical_not(np.isnan(T0.z[i,:]))].shape)
#    i += 1
B = np.sqrt(T0.Bx**2. + T0.By**2. + T0.Bz**2. )

#plt.rcParams['figure.dpi'] = 3000

#plt.savefig('testing.png')
#plt.savefig('testing.pdf')

print("--- %s seconds ---" % (time.time() - start_time))


#np.savetxt("x_vip4+con20int_aligned_eq_phi0=69.2_Elong_501pts_5-10.csv", T0.x, delimiter=",")

#np.savetxt("y_vip4+con20int_aligned_eq_phi0=69.2_Elong_501pts_5-10.csv", T0.y, delimiter=",")

#np.savetxt("z_vip4+con20int_aligned_eq_phi0=69.2_Elong_501pts_5-10.csv", T0.z, delimiter=",")


np.savetxt("anti_aligned_x_jrm33+con20_integral_701_3.5-10.5.csv", T0.x, delimiter=",")



np.savetxt("anti_aligned_y_jrm33+con20_integral_701_3.5-10.5.csv", T0.y, delimiter=",")



np.savetxt("anti_aligned_z_jrm33+con20_integral_701_3.5-10.5.csv", T0.z, delimiter=",")


np.savetxt("Btot_anti_aligned_y_jrm33+con20_integral_701_3.5-10.5.csv", B, delimiter=",")