pro ionosphere_only_eddie_run_diffeq_from_ceq_refs_highres_field_line_interp_include_feh_into4_highest_res_aligned_only_integral


planet='Jupiter'    ;Can also be used for Saturn
fcor  =1.0      ;azimuthal velocity relative to corotation
RJ = 1.; my field line R_fl vals already in RJ so RJ = 1.;7.1492e7 meters     ;Just defines units, set RJ to the radius of
; the planet in whatever units you like
if_psh = 0      ;Whether or not to assume phase space holes
; planetward of the choke point

;Set planet-specific variabies
if ( define_planet(planet, RP, Omega, GM) ) then stop,'Planet '+$
  planet+' is not supported'
Omega *= fcor

;Define species and their values at a reference point (usually the equator)

;Number of species
nspec = 2;9 with feh

;Create an array of structures
species0 = replicate({species,name:'O+',m:16.0d,q:1.0d,T:50.0d,A:1.0d,$
  kappa:100.d,n:1800.d,type:'Aniso_Maxwellian'},$
  nspec)

;Name species, cold/core electrons should always be last species
species0.name= ['H+','e-']
;Core electrons should always be the last species

;Define species' mass, charge and temperature, for distributions whose
; temperature varies with latitude, e.g. anisotropic maxwellians
; or kappa distributions, this the temperature is the temperature
; at the reference point

species0.m=[1.d,1d/1837.d]
species0.q=[1.d,-1.d]

;Define species' densities and temperature at the equator
species0.T=[0.31d,0.31d]
species0.n=[2d5,2d5]
species0.A=[5.0d,5.0d]
;Dougherty et al., 2017 (tables 2 and 3) for L=6

;species0[8].kappa = 2.4       ;Meyer-Vernet et al., 1995
;species0[8].type = 'Kappa'
;species0[7].A=species0[7].T/species0[0].T
;species0[7].T=species0[0].T
;species0[7].type='Aniso_Maxwellian'
;species0[8].type='Aniso_Maxwellian'
;Give hot O+ the same parallel temperature as cold O+



longitude0 = 74.62d ; value of reference point longitudeat z=0

longitude0w = 360.d -  longitude0 

nrfls = 601L;long64(291) ; numbers of field lines from 5-10 RJ (51=0.1RJ uniform resolution, 501points = 0.01RJ uniform resolution)
nphifls = 1L; long64(360) ;




openr,1,'int_aligned_x_jrm33+con20analytic_601_4-10.csv'
x_fl=dblarr(1000, nphifls*nrfls ) ;  1000 points along field line including nans...... and nfls field lines
readf,1,x_fl ;RJ
close,1


openr,1,'int_aligned_y_jrm33+con20analytic_601_4-10.csv'
y_fl=dblarr(1000,  nphifls*nrfls ) ;  1000 points along field line including nans...... and nfls field lines
readf,1,y_fl ;RJ
close,1



openr,1,'int_aligned_z_jrm33+con20analytic_601_4-10.csv'
z_fl=dblarr(1000,   nphifls*nrfls ) ;  1000 points along field line including nans...... and nfls field lines
readf,1,z_fl ;RJ
close,1

;CD, '/Users/enerney/Desktop/diff_eq_code_09-22-22_FJC_BETA'

;idx_nonfinite=where(finite(x_fl) eq 0) ; nonfinite elements
;x_fl[idx_nonfinite] = 0.0

;idx_nonfinite=where(finite(y_fl) eq 0) ; nonfinite elements
;y_fl[idx_nonfinite] = 0.0

;idx_nonfinite=where(finite(z_fl) eq 0) ; nonfinite elements
;z_fl[idx_nonfinite] = 0.0

;write_csv,'x_jrm33+con20int_aligned_eq_phi0=74.62_Elong_501pts_5-10_0fills.csv',x_fl
;write_csv,'y_jrm33+con20int_aligned_eq_phi0=74.62_Elong_501pts_5-10_0fills.csv',y_fl
;write_csv,'z_jrm33+con20int_aligned_eq_phi0=74.62_Elong_501pts_5-10_0fills.csv',z_fl


;stop
;print,x_fl(*,0)
;idx_finite=where(finite(x_fl(*,0)) eq 1)

;x_fl = x_fl(idx_finite)
;y_fl = y_fl(idx_finite)
;z_fl = z_fl(idx_finite)
;print,'here'
;print,x_fl(idx_finite,0)
;stop
rho_fl = Sqrt(x_fl^2 + y_fl^2) ;rho systemIII RJ
R_fl =  Sqrt(x_fl^2 + y_fl^2 + z_fl^2) ;spherical r systemIII RJ
lat_fl = (180.d/!pi)*asin(z_fl/R_fl) ; sysIII lat degrees
phi_fl = (180.d/!pi)*atan(y_fl,x_fl) ; elong degrees



help,x_fl
help,y_fl
help,z_fl

help,rho_fl

help,r_fl


help, phi_fl



stop
phi_fl_wlong = 360.d - phi_fl ; wlong degrees


nnn=601

L00=4.0d



r_00 = 0.01d*dindgen(601) + 4.0d;[0.01*findgen(150) + 4.5, 0.1*findgen(41) + 6.];

dlatfl = 0.1d

lat_int_fl_max = 70.0d

lat_int_fl_min = -lat_int_fl_max



npointsfieldline = Round((lat_int_fl_max - lat_int_fl_min)/dlatfl) + 1

;print,npointsfieldline 

;stop
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



lat_int_fl = lat_int_fl_min + dlatfl*findgen(npointsfieldline)

;i_lat = [10*indgen(20)+200,indgen(200)+391,10*indgen(21)+591]

;lat_int_fl=lat_int_fl[i_lat]
;npointsfieldline = n_elements(lat_int_fl)



print,lat_int_fl
print,npointsfieldline
stop
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


n_out = replicate({densities,sp:0.0,s2p:0.0,s3p:0.0,op:0.0,oph:0.0,o2p:0.0,hp:0.0,nap:0.0,eh:0.0,elec:0.0},$
nrfls*nphifls,npointsfieldline) ;
 
 T_out = replicate({Temps,sp:0.0,s2p:0.0,s3p:0.0,op:0.0,oph:0.0,o2p:0.0,hp:0.0,nap:0.0,eh:0.0,elec:0.0},$
  nrfls*nphifls,npointsfieldline) ; 


rho_out = fltarr(nrfls*nphifls,npointsfieldline)

z_out = fltarr(nrfls*nphifls,npointsfieldline)


x_out = fltarr(nrfls*nphifls,npointsfieldline)

y_out = fltarr(nrfls*nphifls,npointsfieldline)


Bz_out = fltarr(nrfls*nphifls,npointsfieldline)


Bx_out = fltarr(nrfls*nphifls,npointsfieldline)

By_out = fltarr(nrfls*nphifls,npointsfieldline)




openr,1,'Bx_int_aligned_jrm33+con20analytic_601_4-10.csv'
Bx_fl=dblarr(1000, nphifls*nrfls ) ;  1000 points along field line including nans...... and nfls field lines
readf,1,Bx_fl ;RJ
close,1

openr,1,'By_int_aligned_jrm33+con20analytic_601_4-10.csv'
By_fl=dblarr(1000, nphifls*nrfls ) ;  1000 points along field line including nans...... and nfls field lines
readf,1,By_fl ;RJ
close,1

openr,1,'Bz_int_aligned_jrm33+con20analytic_601_4-10.csv'
Bz_fl=dblarr(1000, nphifls*nrfls ) ;  1000 points along field line including nans...... and nfls field lines
readf,1,Bz_fl ;RJ
close,1



;nrfls,nphifls
for j = 0, nrfls*nphifls - 1 do begin
  
;for k = 0, nphifls-1 do begin
  
  
;nrfls = 601L; numbers of field lines from 4-10 RJ  601points = 0.01RJ uniform resolution
;nphifls = 1L;

j = 191; rho =  0.01*findgen(601) + 4.0, rho(191) = 5.91


  
print,'j = ', j , ' of ',nrfls*nphifls - 1;,'k = ', k , ' of ',nphifls-1 

  
  
  idx_finite=where(finite(x_fl(*,j)) eq 1) ;finite elements
  
  
  
  
; radius, latitude, west longitude, |B|, with latitude
  ; and longitude in degrees

rho_int_fl = interpol(Rho_fl(idx_finite,j),lat_fl(idx_finite,j),lat_int_fl)

z_int_fl = interpol(z_fl(idx_finite,j),lat_fl(idx_finite,j),lat_int_fl)
x_int_fl = interpol(x_fl(idx_finite,j),lat_fl(idx_finite,j),lat_int_fl)
y_int_fl = interpol(y_fl(idx_finite,j),lat_fl(idx_finite,j),lat_int_fl)

Bz_int_fl = interpol(Bz_fl(idx_finite,j),lat_fl(idx_finite,j),lat_int_fl)
Bx_int_fl = interpol(Bx_fl(idx_finite,j),lat_fl(idx_finite,j),lat_int_fl)
By_int_fl = interpol(By_fl(idx_finite,j),lat_fl(idx_finite,j),lat_int_fl)

B_int_fl = Sqrt(Bx_int_fl^2d + By_int_fl^2d + Bz_int_fl^2d)

r_int_fl = interpol(R_fl(idx_finite,j),lat_fl(idx_finite,j),lat_int_fl)

phi_int_fl_wlong = interpol(phi_fl_wlong(idx_finite,j),lat_fl(idx_finite,j),lat_int_fl)

;Poss = [R_fl(idx_finite,j),lat_fl(idx_finite,j),phi_fl_wlong(idx_finite,j),Replicate(1.,n_elements(R_fl(idx_finite,j)))]
;Poss = [[R_fl(*,j)],[lat_fl(*,j)],[phi_fl_wlong(*,j)],[Replicate(1.,n_elements(R_fl(*,j)))]]
Poss = [[r_int_fl],[lat_int_fl],[phi_int_fl_wlong],[Replicate(1.,n_elements(r_int_fl))]]


idx_ionosph = where(abs((r_int_fl - 1d)) lt 0.01d)

idx_ionosph = 1380;idx_ionosph[0] ; picking first place where true

r_fl_ref = r_int_fl[idx_ionosph]

;if (r_fl_ceq ge 4.) and (r_fl_ceq le 10.) then begin

lat_fl_ref = lat_int_fl[idx_ionosph]

longw_fl_ref = phi_int_fl_wlong[idx_ionosph]

B_fl_ref = B_int_fl[idx_ionosph]

poss = transpose(poss)
Pos0 = [r_fl_ref,lat_fl_ref,longw_fl_ref, B_fl_ref];,1.]; [r_00[j], 0., longitude0w, 1.] ; ref points !!!!!!!!!!!!!!




;;;; wrong fix below

;      ['O+','O++','S+','S++','S+++','H+','Na+','O+(hot)','e-']
;dens and temp at ref points for me at cent-eq

;ti000 = interpol(ti0,r_00,r_fl_ceq)
;species0.T=[ti000,ti000,ti000,ti000,ti000,interpol(thp0,r_00,r_fl_ceq),ti000,interpol(toph0,r_00,r_fl_ceq),interpol(teh0,r_00,r_fl_ceq),interpol(tec0,r_00,r_fl_ceq)] ; correspond to r_00[j]



;species0.n=[interpol(nop0,r_00,r_fl_ceq),interpol(no2p0,r_00,r_fl_ceq),interpol(nsp0,r_00,r_fl_ceq),interpol(ns2p0,r_00,r_fl_ceq),interpol(ns3p0,r_00,r_fl_ceq),interpol(nhp0,r_00,r_fl_ceq),interpol(nnap0,r_00,r_fl_ceq),interpol(noph0,r_00,r_fl_ceq),interpol(neh0,r_00,r_fl_ceq),interpol(nec0,r_00,r_fl_ceq)] ; correspond to r_00[j]




;;;

npoints = n_elements(Poss[0,*])
n_ions = fltarr(nspec,npoints)
phis = fltarr(npoints)

BRATIO = B_int_fl/B_fl_ref; replicate(1.,npoints) ; just so it doesn't break for now bc using isotropic so this doesn't matter A=1 so B doesn't matter
;stop
;Set electron density at reference point for charge neutrality
; This isn't really necessary since diff_eq() does it internally,
; but it couldn't hurt to make sure the variable in main() is also
; charge neutral...
species0[nspec-1].n=$
  -total(species0[0:nspec-2].q*species0[0:nspec-2].n)$ ; was -2 now -3 because not inclding hot e in this because something something hot protons
  /species0[nspec-1].q
;species0[nspec-1].n=$
;  -total(species0[0:nspec-2].q*species0[0:nspec-2].n)$
 ; /species0[nspec-1].q

;Set species0.A and species0.kappa in a similar way if you're using one of
; those distributions.
;Set species0.type to something else if you don't want a Maxwellian. THis can
; be anything, but there must be a proceedure <type>_density.pro which
; returns density of a species

;Print out the species defined above
;for i=0,nspec-1 do print,species0[i]

;Get the positions where you want to calculate the density. These points
; must all be along the same field line. You can replace this with
; any field tracing routing you like. Poss is a 4xN assay, with
; elements of radius, latitude, west longitude, |B|, with latitude
; and longitude in degrees. Pos0 is the reference point, usually
; the equator, and bratio is the local field strength over the field
; strength at the reference point.
; Note: If if_psh=1 then Poss[1,*] (latitude) must be monotonicly
; increasing and there should be no repeated points
;restore,'field_line_L6.sav'





;get the density at all points along the
for i = 0,npoints - 1 do begin
  if (Poss[0,i] gt 1.) then begin
    
  
  n_ions[*,i]=diff_eq_eddie_ionosphere(Poss[0:2,i],Pos0[0:2],Bratio[i],species0,$
    phi,posc,nchoke,phic,planet=planet,fcor=fcor,$
    if_psh=if_psh,RJ=RJ)
  ;See diff_eq for definition of inputs and outours
  phis[i]=phi
  endif else begin
    n_ions[*,i] = replicate(0.,nspec,1)
  endelse
  
endfor




;      ['O+','O++','S+','S++','S+++','H+','Na+','O+(hot)','e-']
;rho_out[j,k,*] = rho_int_fl

z_out[j,*] = z_int_fl
x_out[j,*] = x_int_fl
y_out[j,*] = y_int_fl

Bz_out[j,*] = Bz_int_fl
Bx_out[j,*] = Bx_int_fl
By_out[j,*] = By_int_fl

;write_csv,'z_outlat-70to70_601_3D_integral_only5.91fl_aligned.txt',reform(z_out[j,*])
;write_csv,'x_outlat-70to70_601_3D_integral_only5.91fl_aligned.txt',reform(x_out[j,*])
;write_csv,'y_outlat-70to70_601_3D_integral_only5.91fl_aligned.txt',reform(y_out[j,*])


write_csv,'phi_volts_ionosphere_only_A=5_both.txt',phis




write_csv,'nhp_ionosphere_only_A=5_both.txt',n_ions[0,*]
write_csv,'nec_ionosphere_only_A=5_both.txt',n_ions[1,*]


stop

;p=plot(poss[1,*],n_ions[nspec-1,*],yrange=[1e-1,1e4],xtitle='Latitude [deg]',$
;   ytitle='Density [cm$^-3$]',color='blue',/ylog,xrange=[-90,90])
 ; p2=plot(poss[1,*],n_ions[0,*],color='red',overplot=p)
 ; p3=plot(poss[1,*],n_ions[1,*],overplot=p)
  ;stop;,'.c to continue (and close plot)'

endfor



;write_csv,'r_out_voyager_diffeq_jr09_can_frank_aligned.txt',r_lat_out.r

;write_csv,'lat_out_voyager_diffeq_jr09_can_frank_aligned.txt',r_lat_out.lat

;write_csv,'Lshells_for_voyager_diffeq_jr09_can_frank_aligned.txt',r_00



;help,n_out.sp

;T_out[j,*].elec 
idx_nonfinite=where(finite(T_out.elec) eq 0) ; nonfinite elements
if (idx_nonfinite(0) ge 0) then T_out[idx_nonfinite].elec = 0.0

idx_nonfinite=where(finite(T_out.eh) eq 0) ; nonfinite elements
if (idx_nonfinite(0) ge 0) then T_out[idx_nonfinite].eh = 0.0

idx_nonfinite=where(finite(T_out.sp) eq 0) ; nonfinite elements
if (idx_nonfinite(0) ge 0) then T_out[idx_nonfinite].sp = 0.0




idx_nonfinite=where(finite(T_out.oph) eq 0) ; nonfinite elements
if (idx_nonfinite(0) ge 0) then T_out[idx_nonfinite].oph = 0.0

idx_nonfinite=where(finite(T_out.hp) eq 0) ; nonfinite elements
if (idx_nonfinite(0) ge 0) then T_out[idx_nonfinite].hp = 0.0

idx_nonfinite=where(finite(n_out.sp) eq 0) ; nonfinite elements
if (idx_nonfinite(0) ge 0) then n_out[idx_nonfinite].sp = 0.0

idx_nonfinite=where(finite(n_out.s2p) eq 0) ; nonfinite elements
if (idx_nonfinite(0) ge 0) then n_out[idx_nonfinite].s2p = 0.0

idx_nonfinite=where(finite(n_out.s3p) eq 0) ; nonfinite elements
if (idx_nonfinite(0) ge 0) then n_out[idx_nonfinite].s3p = 0.0

idx_nonfinite=where(finite(n_out.op) eq 0) ; nonfinite elements
if (idx_nonfinite(0) ge 0) then n_out[idx_nonfinite].op = 0.0

idx_nonfinite=where(finite(n_out.o2p) eq 0) ; nonfinite elements
if (idx_nonfinite(0) ge 0) then n_out[idx_nonfinite].o2p = 0.0

idx_nonfinite=where(finite(n_out.oph) eq 0) ; nonfinite elements
if (idx_nonfinite(0) ge 0) then n_out[idx_nonfinite].oph = 0.0

idx_nonfinite=where(finite(n_out.nap) eq 0) ; nonfinite elements
if (idx_nonfinite(0) ge 0) then n_out[idx_nonfinite].nap = 0.0

idx_nonfinite=where(finite(n_out.hp) eq 0) ; nonfinite elements
if (idx_nonfinite(0) ge 0) then n_out[idx_nonfinite].hp = 0.0

idx_nonfinite=where(finite(n_out.elec) eq 0) ; nonfinite elements
if (idx_nonfinite(0) ge 0) then n_out[idx_nonfinite].elec = 0.0

idx_nonfinite=where(finite(n_out.eh) eq 0) ; nonfinite elements
if (idx_nonfinite(0) ge 0) then n_out[idx_nonfinite].eh = 0.0


idx_nonfinite=where(finite(z_out) eq 0) ; nonfinite elements
if (idx_nonfinite(0) ge 0) then z_out[idx_nonfinite] = 0.0


idx_nonfinite=where(finite(x_out) eq 0) ; nonfinite elements
if (idx_nonfinite(0) ge 0) then x_out[idx_nonfinite] = 0.0


idx_nonfinite=where(finite(y_out) eq 0) ; nonfinite elements
if (idx_nonfinite(0) ge 0) then y_out[idx_nonfinite] = 0.0

idx_nonfinite=where(finite(Bz_out) eq 0) ; nonfinite elements
if (idx_nonfinite(0) ge 0) then Bz_out[idx_nonfinite] = 0.0


idx_nonfinite=where(finite(Bx_out) eq 0) ; nonfinite elements
if (idx_nonfinite(0) ge 0) then Bx_out[idx_nonfinite] = 0.0


idx_nonfinite=where(finite(By_out) eq 0) ; nonfinite elements
if (idx_nonfinite(0) ge 0) then By_out[idx_nonfinite] = 0.0




;idx_nonfinite=where(finite(rho_out) eq 0) ; nonfinite elements
;rho_out[idx_nonfinite] = 0.0


;write_csv,'rho_out_lat-50to50_0.1deginterps_501x360_3D_analytic.txt',rho_out

write_csv,'z_outlat-50to50_601_3D_integral.txt',z_out
write_csv,'x_outlat-50to50_601_3D_integral.txt',x_out
write_csv,'y_outlat-50to50_601_3D_integral.txt',y_out

write_csv,'nsp_out_mymodel1_diffeq_jrm33+Con2020_0fillslat-50to50_601_3D_integral_including_feh.txt',n_out.sp

write_csv,'ns2p_out_mymodel1_diffeq_jrm33+Con2020_0fillslat-50to50_601_3D_integral_including_feh.txt',n_out.s2p

write_csv,'ns3p_out_mymodel1_diffeq_jrm33+Con2020_0fillslat-50to50_601_3D_integral_including_feh.txt',n_out.s3p

write_csv,'nop_out_mymodel1_diffeq_jrm33+Con2020_0fillslat-50to50_601_3D_integral_including_feh.txt',n_out.op

write_csv,'noph_out_mymodel1_diffeq_jrm33+Con2020_0fillslat-50to50_601_3D_integral_including_feh.txt',n_out.oph

write_csv,'no2p_out_mymodel1_diffeq_jrm33+Con2020_0fillslat-50to50_601_3D_integral_including_feh.txt',n_out.o2p

write_csv,'nhp_out_mymodel1_diffeq_jrm33+Con2020_0fillslat-50to50_601_3D_integral_including_feh.txt',n_out.hp

write_csv,'nnap_out_mymodel1_diffeq_jrm33+Con2020_0fillslat-50to50_601_3D_integral_including_feh.txt',n_out.nap

write_csv,'nelec_out_mymodel1_diffeq_jrm33+Con2020_0fillslat-50to50_601_3D_integral_including_feh.txt',n_out.elec

write_csv,'neh_out_mymodel1_diffeq_jrm33+Con2020_0fillslat-50to50_601_3D_integral_including_feh.txt',n_out.eh

write_csv,'Telec_out_mymodel1_diffeq_jrm33+Con2020_0fillslat-50to50_601_3D_integral_including_feh.txt',T_out.elec

write_csv,'Tic_out_mymodel1_diffeq_jrm33+Con2020_0fillslat-50to50_601_3D_integral_including_feh.txt',T_out.sp

write_csv,'Toph_out_mymodel1_diffeq_jrm33+Con2020_0fillslat-50to50_601_3D_integral_including_feh.txt',T_out.oph

write_csv,'Thp_out_mymodel1_diffeq_jrm33+Con2020_0fillslat-50to50_601_3D_integral_including_feh.txt',T_out.hp

write_csv,'Teh_out_mymodel1_diffeq_jrm33+Con2020_0fillslat-50to50_601_3D_integral_including_feh.txt',T_out.eh

;p=plot(poss[1,*],n_ions[nspec-1,*],yrange=[1e-1,1e4],xtitle='Latitude [deg]',$
 ; ytitle='Density [cm$^-3$]',color='blue',/ylog,xrange=[-90,90])
;p2=plot(poss[1,*],n_ions[0,*],color='red',overplot=p)
;p3=plot(poss[1,*],n_ions[1,*],overplot=p)
;,'.c to continue (and close plot)'
;p.close
stop
end