function diffeq_fit_bill_kurth_fingers_vary_only_nhpmix0_and_Thp0,p,x=x,y=y,err=err
COMMON block1,nspmix0,ns2pmix0,ns3pmix0,nopmix0,no2pmix0,nophmix0,nnapmix0,nhpmix0,nfac,tfac,nec0, species0, nspec, lathave,latwant, Poss, Pos0, Bratio, j_ctp

planet = 'Jupiter'    ;Can also be used for Saturn ;
fcor  =1.0d      ;azimuthal velocity relative to corotation
RJ = 1.d; my field line R_fl vals already in RJ so RJ = 1.;7.1492e7 meters     ;Just defines units, set RJ to the radius of
; the planet in whatever units you like
if_psh = 0      ;Whether or not to assume phase space holes
; planetward of the choke point

;Set planet-specific variabies
if ( define_planet(planet, RP, Omega, GM) ) then stop,'Planet '+$
  planet+' is not supported'
Omega *= fcor

;Define species and their values at a reference point (usually the equator)

;Number of species
nspec = 10;9 with feh, was 3 with only sp and op and nec

;Create an array of structures
species0 = replicate({species,name:'O+',m:16.0d,q:1.0d,T:50.0d,A:1.0d,$
  kappa:100.d,n:1800.d,type:'Maxwellian'},$
  nspec)

;Name species, cold/core electrons should always be last species
species0.name= ['O+','O++','S+','S++','S+++','H+','Na+','O+(hot)','eh-','e-']
;species0.name= ['O+','S+','e-']
;Core electrons should always be the last species

;Define species' mass, charge and temperature, for distributions whose
; temperature varies with latitude, e.g. anisotropic maxwellians
; or kappa distributions, this the temperature is the temperature
; at the reference point

species0.m=[16.d,16.d,32.d,32.d,32.d,1.d,23.d,16.d,1d/1837.d,1d/1837.d]
species0.q=[1.d,2.d,1.d,2.d,3.d,1.d,1.d,1.d,-1.d,-1.d]

;species0.m=[16.d,32.d,1d/1837.d]
;species0.q=[1.d,1.d,-1.d]

;Define species' densities and temperature at the equator
species0.T=[79.3d,79.3d,79.3d,79.3d,79.3d,79.3d,94.1d,362.d,46.d,4.6d]
species0.n=[592.d,76.3d,163.d,538.d,90.7d,50.6d,97.2d,134.d,2200.d*0.001d,2200.d]
species0.A=[1.0d,1.0d,1.0d,1.0d,1.0d,1.0d,1.0d,1.0d,1.0d,1.0d]



nhpmixneww = p[0]

thpneww = p[1]

necneww = p[2]

deltanhpmix0 = nhpmixneww - nhpmix0


nopmixneww = nopmix0 - deltanhpmix0

;;!!!!!!!!!!!!!!!!!!!!!!!!!!!!
;;;!!!!!!!!!!!!!!!!!!!!!!!!!!!!
;;;!!!!!!!!!!!!!!!!!!!!!!!!!!!!
;;;!!!!!!!!!!!!!!!!!!!!!!!!!!!!
species0.T = [1d,1d,1d,1d,1d,thpneww,1d,1d,1d,1d]
;species0.T = replicate(1d,nspec);[ti000,ti000,ti000,ti000,ti000,interpol(thp0,r_00,r_fl_ceq),ti000,tempfac*interpol(toph0,r_00,r_fl_ceq),interpol(teh0,r_00,r_fl_ceq),tempfac*interpol(tec0,r_00,r_fl_ceq)] ; correspond to r_00[j]
;species0[5].T = thpnew ;interpol(thp0,r_00,5.27d)
;comment back in!!!!!!!!!

;species0.n=[interpol(nfac*nop0,r_00,r_fl_ceq),interpol(nfac*no2p0,r_00,r_fl_ceq),interpol(nfac*nsp0,r_00,r_fl_ceq),interpol(nfac*ns2p0,r_00,r_fl_ceq),interpol(nfac*ns3p0,r_00,r_fl_ceq),interpol(nfac*nhp0,r_00,r_fl_ceq),interpol(nfac*nnap0,r_00,r_fl_ceq),interpol(nfac*noph0,r_00,r_fl_ceq),interpol(neh0,r_00,r_fl_ceq),interpol(nfac*nec0,r_00,r_fl_ceq)] ; correspond to r_00[j]
;species0.name= ['O+','O++','S+','S++','S+++','H+','Na+','O+(hot)','eh-','e-']


;species0.n = [nfac*nopmix0*nec0[j_ctp],nfac*no2pmix0*nec0[j_ctp],nfac*nspmix0*nec0[j_ctp],nfac*ns2pmix0*nec0[j_ctp],nfac*ns3pmix0*nec0[j_ctp],nfac*nhpmix0*nec0[j_ctp],nfac*nnapmix0*nec0[j_ctp],nfac*nophmix0*nec0[j_ctp],0d,nfac*nec0[j_ctp]]
species0.n = [(nopmixneww + nophmix0 + 2d*ns2pmix0 + 3d*ns3pmix0)*necneww,no2pmix0*necneww,nspmix0*necneww,0d,0d,nhpmixneww*necneww,nnapmix0*necneww,0d,0d,necneww]

species0[nspec-1].n=$
  -total(species0[0:nspec-3].q*species0[0:nspec-3].n)$ ; was -2 now -3 because not inclding hot e in this because something something hot protons
  /species0[nspec-1].q

  npoints = n_elements(Poss[0,*])
  n_ions = dblarr(nspec,npoints)
  phis = dblarr(npoints)

  ;get the density at all points along the
  for i = 0,npoints - 1 do begin
    if (Poss[0,i] gt 1d) then begin
      ;!!!!!!!!!!!!!!!!!!
      ; !!!!!!set back to diff_eq_eddie for standard no hot e in charge neutrality !!!!!!!
      n_ions[*,i]=diff_eq_eddie(Poss[0:2,i],Pos0[0:2],Bratio[i],species0,$
        phi,posc,nchoke,phic,planet='Jupiter',fcor=1d,$
        if_psh=0,RJ=1d)
      ;See diff_eq for definition of inputs and outours
      phis[i]=phi
    endif else begin
      n_ions[*,i] = replicate(0d,nspec,1)
    endelse

  endfor

  ;n_out[j,*].op = reform(reform( n_ions[0,*]),1,npointsfieldline)
 ; n_out[j,*].o2p = reform(reform( n_ions[1,*]),1,npointsfieldline)
 ; n_out[j,*].sp =  reform(reform( n_ions[2,*]),1,npointsfieldline);reform(reform(n_ions[2,*]),1,npointsfieldline)
 ; n_out[j,*].s2p = reform(reform( n_ions[3,*]),1,npointsfieldline)
 ; n_out[j,*].s3p = reform(reform( n_ions[4,*]),1,npointsfieldline)
 ; n_out[j,*].hp =  reform(reform(n_ions[5,*]),1,npointsfieldline)
 ; n_out[j,*].nap = reform(reform( n_ions[6,*]),1,npointsfieldline)
 ; n_out[j,*].oph = reform(reform( n_ions[7,*]),1,npointsfieldline)
  ;n_out[j,*].eh = reform( reform( n_ions[8,*]),1,npointsfieldline)
  ;n_out[j,*].elec = reform(reform( n_ions[9,*]),1,npointsfieldline);reform( reform( n_ions[9,*]),1,npointsfieldline)

model0 = reform(n_ions[9,*])
model = interpol(model0, lathave,latwant)


 result = ( y - model)/err
 
 print,total(result^2d)
 
 ;stop
 return, result
end


pro EDDIE_RUN_DIFFEQ_TEST_COMPARISON
COMMON block1,nspmix0,ns2pmix0,ns3pmix0,nopmix0,no2pmix0,nophmix0,nnapmix0,nhpmix0,nfac,tfac,nec0, species0, nspec, lathave,latwant, Poss, Pos0, Bratio, j_ctp


planet = 'Jupiter'    ;Can also be used for Saturn ;
fcor  =1.0d      ;azimuthal velocity relative to corotation
RJ = 1.d; my field line R_fl vals already in RJ so RJ = 1.;7.1492e7 meters     ;Just defines units, set RJ to the radius of
; the planet in whatever units you like
if_psh = 0      ;Whether or not to assume phase space holes
; planetward of the choke point

;Set planet-specific variabies
if ( define_planet(planet, RP, Omega, GM) ) then stop,'Planet '+$
  planet+' is not supported'
Omega *= fcor

;Define species and their values at a reference point (usually the equator)

;Number of species
nspec = 10;9 with feh, was 3 with only sp and op and nec

;Create an array of structures
species0 = replicate({species,name:'O+',m:16.0d,q:1.0d,T:50.0d,A:1.0d,$
  kappa:100.d,n:1800.d,type:'Maxwellian'},$
  nspec)

;Name species, cold/core electrons should always be last species
species0.name= ['O+','O++','S+','S++','S+++','H+','Na+','O+(hot)','eh-','e-']
;species0.name= ['O+','S+','e-']
;Core electrons should always be the last species

;Define species' mass, charge and temperature, for distributions whose
; temperature varies with latitude, e.g. anisotropic maxwellians
; or kappa distributions, this the temperature is the temperature
; at the reference point

species0.m=[16.d,16.d,32.d,32.d,32.d,1.d,23.d,16.d,1d/1837.d,1d/1837.d]
species0.q=[1.d,2.d,1.d,2.d,3.d,1.d,1.d,1.d,-1.d,-1.d]

;species0.m=[16.d,32.d,1d/1837.d]
;species0.q=[1.d,1.d,-1.d]

;Define species' densities and temperature at the equator
species0.T=[79.3d,79.3d,79.3d,79.3d,79.3d,79.3d,94.1d,362.d,46.d,4.6d]
species0.n=[592.d,76.3d,163.d,538.d,90.7d,50.6d,97.2d,134.d,2.5375d,2537.5d]
species0.A=[1.0d,1.0d,1.0d,1.0d,1.0d,1.0d,1.0d,1.0d,1.0d,1.0d]
;species0.A=[1.0d,1.0d,1.0d,1.0d,1.0d,5.0d,1.0d,5.0d,5.0d,1.0d]
;
;species0.T=[35d,50d,5d]
;species0.n=[1750d,250d,2000d]
;species0.A=[2.0d,1.0d,1d]
;Dougherty et al., 2017 (tables 2 and 3) for L=6

;species0[8].kappa = 2.4       ;Meyer-Vernet et al., 1995
;species0[8].type = 'Kappa'
;species0[7].A=species0[7].T/species0[0].T
;species0[7].T=species0[0].T
;species0[0].type='Aniso_Maxwellian'

;species0[*].A = 1d
;species0[*].type='fried_egg'
;species0[5].type='Aniso_Maxwellian'
;species0[7].type='Aniso_Maxwellian'
;species0[8].type='Aniso_Maxwellian'
;Give hot O+ the same parallel temperature as cold O+



longitude0 = 74.62d ; value of reference point longitudeat z=0

longitude0w = 360.d -  longitude0 

nrfls = 601L;long64(291) ; numbers of field lines from 5-10 RJ (51=0.1RJ uniform resolution, 501points = 0.01RJ uniform resolution)
nphifls = 1L; long64(360) ;




openr,1,'int_aligned_x_jrm33+con20integral_601_4-10_properly_aligned_phi0=75.85_perfectaligned_for_4.70-9.71.csv'
x_fl=dblarr(1000, nphifls*nrfls ) ;  1000 points along field line including nans...... and nfls field lines
readf,1,x_fl ;RJ
close,1


openr,1,'int_aligned_y_jrm33+con20integral_601_4-10_properly_aligned_phi0=75.85_perfectaligned_for_4.70-9.71.csv'
y_fl=dblarr(1000,  nphifls*nrfls ) ;  1000 points along field line including nans...... and nfls field lines
readf,1,y_fl ;RJ
close,1



openr,1,'int_aligned_z_jrm33+con20integral_601_4-10_properly_aligned_phi0=75.85_perfectaligned_for_4.70-9.71.csv'
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



;stop
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


lat_int_fl = lat_int_fl_min + dlatfl*dindgen(npointsfieldline)


;i_lat = [10*indgen(20)+200,indgen(200)+391,10*indgen(21)+591]

;lat_int_fl=lat_int_fl[i_lat]
;npointsfieldline = n_elements(lat_int_fl)



;print,lat_int_fl
;print,npointsfieldline
;stop
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


n_out = replicate({densities,sp:0.0d,s2p:0.0d,s3p:0.0,op:0.0d,oph:0.0d,o2p:0.0d,hp:0.0d,nap:0.0d,eh:0.0d,elec:0.0d},$
nrfls*nphifls,npointsfieldline) ;
 
 T_out = replicate({Temps,sp:0.0d,s2p:0.0d,s3p:0.0d,op:0.0d,oph:0.0d,o2p:0.0d,hp:0.0d,nap:0.0d,eh:0.0d,elec:0.0d},$
  nrfls*nphifls,npointsfieldline) ; 


rho_out = dblarr(nrfls*nphifls,npointsfieldline)

z_out = dblarr(nrfls*nphifls,npointsfieldline)


x_out = dblarr(nrfls*nphifls,npointsfieldline)

y_out = dblarr(nrfls*nphifls,npointsfieldline)


Bz_out = dblarr(nrfls*nphifls,npointsfieldline)


Bx_out = dblarr(nrfls*nphifls,npointsfieldline)

By_out = dblarr(nrfls*nphifls,npointsfieldline)

B_out = dblarr(nrfls*nphifls,npointsfieldline)


openr,1,'nop_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv'
nop0=dblarr(601) 
readf,1,nop0 
close,1

openr,1,'no2p_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv'
no2p0=dblarr(601)
readf,1,no2p0 
close,1


openr,1,'nsp_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv'
nsp0=dblarr(601)
readf,1,nsp0
close,1

openr,1,'ns2p_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv'
ns2p0=dblarr(601)
readf,1,ns2p0
close,1

openr,1,'ns3p_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv'
ns3p0=dblarr(601)
readf,1,ns3p0
close,1

openr,1,'nhp_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv'
nhp0=dblarr(601)
readf,1,nhp0
close,1

openr,1,'noph_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv'
noph0=dblarr(601)
readf,1,noph0
close,1


openr,1,'nnap_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv'
nnap0=dblarr(601)
readf,1,nnap0
close,1


openr,1,'neh_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv'
neh0=dblarr(601)
readf,1,neh0
close,1

openr,1,'nec_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv'
nec0=dblarr(601)
readf,1,nec0
close,1

openr,1,'Tic_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv'
ti0=dblarr(601)
readf,1,ti0
close,1


openr,1,'Thp_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv'
thp0=dblarr(601)
readf,1,thp0
close,1

openr,1,'Toph_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv'
toph0=dblarr(601)
readf,1,toph0
close,1

openr,1,'Tec_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv'
tec0=dblarr(601)
readf,1,tec0
close,1

openr,1,'Teh_mymodel1_601pts_4-10_densgaussianinside5_temps1inside5.csv'
teh0=dblarr(601)
readf,1,teh0
close,1

openr,1,'Bx_int_aligned_jrm33+con20integral_601_4-10_properly_aligned_phi0=75.85_perfectaligned_for_4.70-9.71.csv'
Bx_fl=dblarr(1000, nphifls*nrfls ) ;  1000 points along field line including nans...... and nfls field lines
readf,1,Bx_fl ;RJ
close,1

openr,1,'By_int_aligned_jrm33+con20integral_601_4-10_properly_aligned_phi0=75.85_perfectaligned_for_4.70-9.71.csv'
By_fl=dblarr(1000, nphifls*nrfls ) ;  1000 points along field line including nans...... and nfls field lines
readf,1,By_fl ;RJ
close,1

openr,1,'Bz_int_aligned_jrm33+con20integral_601_4-10_properly_aligned_phi0=75.85_perfectaligned_for_4.70-9.71.csv'
Bz_fl=dblarr(1000, nphifls*nrfls ) ;  1000 points along field line including nans...... and nfls field lines
readf,1,Bz_fl ;RJ
close,1

openr,1,'Btot_int_aligned_jrm33+con20integral_601_4-10_properly_aligned_phi0=75.85_perfectaligned_for_4.70-9.71.csv'
Btot_fl=dblarr(1000, nphifls*nrfls ) ;  1000 points along field line including nans...... and nfls field lines
readf,1,Btot_fl ;RJ
close,1




sed_data = READ_CSV('Bill_Kurth_fingers_centlat_vs_ne_fig10e.csv', HEADER=SedHeader, $
  N_TABLE_HEADER=1, TABLE_HEADER=SedTableHeader)


Bill_Kurth_fingers_lat_ceq_deg = sed_data.field1

Bill_Kurth_fingers_ne = sed_data.field2


p1=plot(r_00,(nsp0 + 2d*ns2p0 + nop0 )/nec0,xtitle='$\rho_c$ ($R_J$)',ytitle='$\Sigma Z_i n_i / n_{ec}$',title='Only $S^+$, $S^{++}$, and $O^+$ from Nominal Model',xrange=[4.5,7.0])

rho_ceq_for_totchmix = r_00
totchmix_for_apo = (nsp0 + 2d*ns2p0 + nop0 )/nec0

;Pos0 = [r_fl_ceq,lat_fl_ceq,longw_fl_ceq, B_fl_ceq];,1.]; [r_00[j], 0., longitude0w, 1.] ; ref points !!!!!!!!!!!!!!
poss = [6.0d, 10.0d, 74.62d,1d] 
pos0 = [6.0d, 0.0d, 74.62d,1d] 




;n_elements(Poss[0,*])
n_ions = dblarr(nspec)
phis = 1d

BRATIO = 1d; just so it doesn't break for now bc using isotropic so this doesn't matter A=1 so B doesn't matter
;stop
;Set electron density at reference point for charge neutrality
; This isn't really necessary since diff_eq() does it internally,
; but it couldn't hurt to make sure the variable in main() is also
; charge neutral...
species0[nspec-1].n=$
  -total(species0[0:nspec-3].q*species0[0:nspec-3].n)$ ; was -2 now -3 because not inclding hot e in this because something something hot protons
  /species0[nspec-1].q  


  n_ions=diff_eq_eddie(Poss[0:2],Pos0[0:2],1d,species0,$
    phi,posc,nchoke,phic,planet=planet,fcor=fcor,$
    if_psh=if_psh,RJ=RJ)
  
print,n_ions
stop
end
