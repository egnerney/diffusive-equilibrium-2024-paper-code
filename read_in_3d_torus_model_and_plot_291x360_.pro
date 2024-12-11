pro read_in_3d_torus_model_and_plot_291x360_
tic

 numflsss = 291L*360L
  openr,1,'nelec_out_mymodel1_diffeq_jrm33+Con2020_0fillslat-30to30_0.1deginterps_outto20_291x360_3D_analytic_including_feh_fixed_shifted.txt'
  nel=fltarr(numflsss,241)
  readf,1,nel
  close,1
  
  openr,1,'neh_out_mymodel1_diffeq_jrm33+Con2020_0fillslat-30to30_0.1deginterps_outto20_291x360_3D_analytic_including_feh_fixed_shifted.txt'
  neh=fltarr(numflsss,241)
  readf,1,neh
  close,1
  
  openr,1,'nsp_out_mymodel1_diffeq_jrm33+Con2020_0fillslat-30to30_0.1deginterps_outto20_291x360_3D_analytic_including_feh_fixed_shifted.txt'
  nsp=fltarr(numflsss,241)
  readf,1,nsp
  close,1
  
  openr,1,'ns2p_out_mymodel1_diffeq_jrm33+Con2020_0fillslat-30to30_0.1deginterps_outto20_291x360_3D_analytic_including_feh_fixed_shifted.txt'
  ns2p=fltarr(numflsss,241)
  readf,1,ns2p
  close,1
  
  openr,1,'ns3p_out_mymodel1_diffeq_jrm33+Con2020_0fillslat-30to30_0.1deginterps_outto20_291x360_3D_analytic_including_feh_fixed_shifted.txt'
  ns3p=fltarr(numflsss,241)
  readf,1,ns3p
  close,1
  
  
  openr,1,'nop_out_mymodel1_diffeq_jrm33+Con2020_0fillslat-30to30_0.1deginterps_outto20_291x360_3D_analytic_including_feh_fixed_shifted.txt'
  nop=fltarr(numflsss,241)
  readf,1,nop
  close,1
  
  openr,1,'noph_out_mymodel1_diffeq_jrm33+Con2020_0fillslat-30to30_0.1deginterps_outto20_291x360_3D_analytic_including_feh_fixed_shifted.txt'
  noph=fltarr(numflsss,241)
  readf,1,noph
  close,1
  
  openr,1,'no2p_out_mymodel1_diffeq_jrm33+Con2020_0fillslat-30to30_0.1deginterps_outto20_291x360_3D_analytic_including_feh_fixed_shifted.txt'
  no2p=fltarr(numflsss,241)
  readf,1,no2p
  close,1
  
  openr,1,'nhp_out_mymodel1_diffeq_jrm33+Con2020_0fillslat-30to30_0.1deginterps_outto20_291x360_3D_analytic_including_feh_fixed_shifted.txt'
  nhp=fltarr(numflsss,241)
  readf,1,nhp
  close,1
  
  openr,1,'nnap_out_mymodel1_diffeq_jrm33+Con2020_0fillslat-30to30_0.1deginterps_outto20_291x360_3D_analytic_including_feh_fixed_shifted.txt'
  nnap=fltarr(numflsss,241)
  readf,1,nnap
  close,1
  
  
  
  openr,1,'Telec_out_mymodel1_diffeq_jrm33+Con2020_0fillslat-30to30_0.1deginterps_outto20_291x360_3D_analytic_including_feh_fixed_shifted.txt'
  tec=fltarr(numflsss,241)
  readf,1,tec
  close,1
  
  openr,1,'Teh_out_mymodel1_diffeq_jrm33+Con2020_0fillslat-30to30_0.1deginterps_outto20_291x360_3D_analytic_including_feh_fixed_shifted.txt'
  teh=fltarr(numflsss,241)
  readf,1,teh
  close,1
  
  openr,1,'Tic_out_mymodel1_diffeq_jrm33+Con2020_0fillslat-30to30_0.1deginterps_outto20_291x360_3D_analytic_including_feh_fixed_shifted.txt'
  tico=fltarr(numflsss,241)
  readf,1,tico
  close,1
  
  openr,1,'Thp_out_mymodel1_diffeq_jrm33+Con2020_0fillslat-30to30_0.1deginterps_outto20_291x360_3D_analytic_including_feh_fixed_shifted.txt'
  thp=fltarr(numflsss,241)
  readf,1,thp
  close,1
  
  openr,1,'Toph_out_mymodel1_diffeq_jrm33+Con2020_0fillslat-30to30_0.1deginterps_outto20_291x360_3D_analytic_including_feh_fixed_shifted.txt'
  toph=fltarr(numflsss,241)
  readf,1,toph
  close,1
  
  


  
  
  openr,1,'x_outlat-30to30_0.1deginterps_outto20_291x360_3D_analytic_fixed.txt'
  x=fltarr(numflsss,241)
  readf,1,x
  close,1
  
  openr,1,'y_outlat-30to30_0.1deginterps_outto20_291x360_3D_analytic_fixed.txt'
  y=fltarr(numflsss,241)
  readf,1,y
  close,1
  
  openr,1,'z_outlat-30to30_0.1deginterps_outto20_291x360_3D_analytic_fixed.txt'
  z=fltarr(numflsss,241)
  readf,1,z
  close,1
  
  rho = Sqrt(x^2 + y^2) ;rho systemIII RJ
  r =  Sqrt(x^2 + y^2 + z^2) ;spherical r systemIII RJ
  lat = (180./!pi)*asin(z/r) ; sysIII lat degrees
  phi = (180./!pi)*atan(y,x) ; elong degrees


;for i =0,35 do p1=plot(lat[i*10,*],nel[i*10,*])
;j=0. ; rho = 5.0 , j=100 rho = 6.0, j=500 rho 10.0

;want to remap to 3D array instead of 2
 ;for i=0,10 do print,r(360.*j + i*10.,500)
 nrfls = 291L
 nphifls = 360L
 
 r_rebin =  fltarr(nrfls,nphifls,241)
phi_rebin = fltarr(nrfls,nphifls,241)
rho_rebin =  fltarr(nrfls,nphifls,241)
lat_rebin =  fltarr(nrfls,nphifls,241)


x_rebin =  fltarr(nrfls,nphifls,241)
y_rebin =  fltarr(nrfls,nphifls,241)
z_rebin =   fltarr(nrfls,nphifls,241)

ne_rebin = fltarr(nrfls,nphifls,241)
neh_rebin = fltarr(nrfls,nphifls,241)
nsp_rebin = fltarr(nrfls,nphifls,241)
ns2p_rebin = fltarr(nrfls,nphifls,241)
ns3p_rebin = fltarr(nrfls,nphifls,241)
nop_rebin = fltarr(nrfls,nphifls,241)
no2p_rebin = fltarr(nrfls,nphifls,241)
noph_rebin = fltarr(nrfls,nphifls,241)
nhp_rebin = fltarr(nrfls,nphifls,241)
nnap_rebin = fltarr(nrfls,nphifls,241)

Tec_rebin = fltarr(nrfls,nphifls,241)
Teh_rebin = fltarr(nrfls,nphifls,241)
Tic_rebin = fltarr(nrfls,nphifls,241)
Thp_rebin = fltarr(nrfls,nphifls,241)
Toph_rebin = fltarr(nrfls,nphifls,241)

;fl_idx = 80000

;p1=plot(lat(fl_idx,*),nel(fl_idx,*),xtitle='$Lat_{III}$ (Degrees)',ytitle='$n_{ec}$ ($cm^{-3}$)')
;p2=plot(lat(fl_idx,*),nel(fl_idx,*) - (nsp(fl_idx,*) + 2d*ns2p(fl_idx,*) + 3d*ns3p(fl_idx,*) + nop(fl_idx,*) + 2d*no2p(fl_idx,*) + noph(fl_idx,*) + nhp(fl_idx,*) + nnap(fl_idx,*)),xtitle='$Lat_{III}$ (Degrees)',ytitle='$n_{ec}$ - $\Sigma n_{i}q_{i}$ ($cm^{-3}$)')
;p3=plot(lat(fl_idx,*),(nel(fl_idx,*) - (nsp(fl_idx,*) + 2d*ns2p(fl_idx,*) + 3d*ns3p(fl_idx,*) + nop(fl_idx,*) + 2d*no2p(fl_idx,*) + noph(fl_idx,*) + nhp(fl_idx,*) + nnap(fl_idx,*)))/nel(fl_idx,*),xtitle='$Lat_{III}$ (Degrees)',ytitle='delta_n/n')


;stop

 for i=0L,290L do begin
 
  for j=0L, 359L do begin
    r_rebin(i,j,*) = r(360L*i + j,*)
    
    rho_rebin(i,j,*) = rho(360L*i + j,*)
    phi_rebin(i,j,*) = phi(360L*i + j,*)
    lat_rebin(i,j,*) = lat(360L*i + j,*)
    x_rebin(i,j,*) = x(360L*i + j,*)
    y_rebin(i,j,*) = y(360L*i + j,*)
    z_rebin(i,j,*) = z(360L*i + j,*)
    ne_rebin(i,j,*) = nel(360L*i + j,*)
    neh_rebin(i,j,*) = neh(360L*i + j,*)
    nsp_rebin(i,j,*) = nsp(360L*i + j,*)
    ns2p_rebin(i,j,*) = ns2p(360L*i + j,*)
    ns3p_rebin(i,j,*) = ns3p(360L*i + j,*)
    nop_rebin(i,j,*) = nop(360L*i + j,*)
    no2p_rebin(i,j,*) = no2p(360L*i + j,*)
    noph_rebin(i,j,*) = noph(360L*i + j,*)
    nhp_rebin(i,j,*) = nhp(360L*i + j,*)
    nnap_rebin(i,j,*) = nnap(360L*i + j,*)
    
    Tec_rebin(i,j,*) = Tec(360L*i + j,*)
    Teh_rebin(i,j,*) = Teh(360L*i + j,*)
     Tic_rebin(i,j,*) = Tico(360L*i + j,*)
      Thp_rebin(i,j,*) = Thp(360L*i + j,*)
       Toph_rebin(i,j,*) = Toph(360L*i + j,*)
     
    
    
  endfor
 endfor
; for i=0,10 do print,r(360.*j + i*10.,500)
;for i=0L,10L do print,r(360L*j + i*10L,120)????
;fl_idx = 80000

rho_fl_idx = 251 ; 6ish
phi_fl_idx = 75 ; 75 degrees approx aligned

rho_fl_idx_vec = indgen(40) + 251

;for i=0,39 do begin
  
;rho_fl_idx = rho_fl_idx_vec(i)

;p1=plot(lat_rebin(rho_fl_idx,phi_fl_idx,*),ne_rebin(rho_fl_idx,phi_fl_idx,*),xtitle='$Lat_{III}$ (Degrees)',ytitle='$n_{ec}$ ($cm^{-3}$)')
;p2=plot(lat_rebin(rho_fl_idx,phi_fl_idx,*),ne_rebin(rho_fl_idx,phi_fl_idx,*) - (nsp_rebin(rho_fl_idx,phi_fl_idx,*) + 2d*ns2p_rebin(rho_fl_idx,phi_fl_idx,*) + 3d*ns3p_rebin(rho_fl_idx,phi_fl_idx,*) + nop_rebin(rho_fl_idx,phi_fl_idx,*) + 2d*no2p_rebin(rho_fl_idx,phi_fl_idx,*) + noph_rebin(rho_fl_idx,phi_fl_idx,*) + nhp_rebin(rho_fl_idx,phi_fl_idx,*) + nnap_rebin(rho_fl_idx,phi_fl_idx,*)),xtitle='$Lat_{III}$ (Degrees)',ytitle='$n_{ec}$ - $\Sigma n_{i}q_{i}$ ($cm^{-3}$)')
;  p3=plot(lat_rebin(rho_fl_idx,phi_fl_idx,*),(ne_rebin(rho_fl_idx,phi_fl_idx,*) - (nsp_rebin(rho_fl_idx,phi_fl_idx,*) + 2d*ns2p_rebin(rho_fl_idx,phi_fl_idx,*) + 3d*ns3p_rebin(rho_fl_idx,phi_fl_idx,*) + nop_rebin(rho_fl_idx,phi_fl_idx,*) + 2d*no2p_rebin(rho_fl_idx,phi_fl_idx,*) + noph_rebin(rho_fl_idx,phi_fl_idx,*) + nhp_rebin(rho_fl_idx,phi_fl_idx,*) + nnap_rebin(rho_fl_idx,phi_fl_idx,*)))/ne_rebin(rho_fl_idx,phi_fl_idx,*),xtitle='$Lat_{III}$ (Degrees)',ytitle='delta_n/n')


;stop
;endfor

x = x_rebin
y = y_rebin
z = z_rebin
nel = ne_rebin
neh = neh_rebin
nsp = nsp_rebin
ns2p = ns2p_rebin
ns3p = ns3p_rebin
nop = nop_rebin
no2p = no2p_rebin
noph = noph_rebin
nhp = nhp_rebin
nnap = nnap_rebin

Tec = Tec_rebin
Teh = Teh_rebin
Ti = Tic_rebin
Thp = Thp_rebin
Toph = Toph_rebin


toc
tic
SAVE,x,y,z,nel,neh,nsp,ns2p,ns3p,nop,no2p,noph,nhp,nnap,Tec,teh,Ti,Toph,Thp, FILENAME = 'Output_rebin3D_smaller_291x360x241_including_feh_fixed_shifted.sav',/verbose
toc
stop
end