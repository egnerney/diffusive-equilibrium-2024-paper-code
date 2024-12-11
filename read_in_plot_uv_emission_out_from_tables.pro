pro read_in_plot_UV_emission_out_from_tables

Restore, 'Mymodelv1_3D_JRM33+con20_cartesian_irrcart_mapped_no_interp_2degree_phi_into_rhoIII_4.5_281x281x129.sav',/verbose
  dx1 = 0.1
  nx1 = 40

  dx2 = 0.05;0.025
  nx2 = 78;390 is for 0.01 ; 156 is 0.025 ; 78 for 0.05 ;
  ;429 for x and y for 0.025, 281 for 0.05
  nx3 = 44 ; middle bits at 0.1


  xgrid1d = [ dx1*findgen(nx1) - 10. , dx2*findgen(nx2) - 6. , dx1*findgen(nx3) - 2.1 , dx2*findgen(nx2) + 2.1 , dx1*findgen(nx1+1) + 6. ]



  ;z_min = -2.5;
  ;z_max = 2.5;
  ;z_step_size = 0.05 ;  RJ
  ;n_z = Round((z_max - z_min)/z_step_size) + 1
  dz1 = 0.1
  nz1 = 12

  dz2 = 0.025
  nz2 = 104;104 for 0.025 ; 80 for 0.05; 260 for 0.01 for -1.3 to 1.3 so 2.6RJ total
  ;285 elements for zgrid if dz=0.01 for z= -1.3 to 1.3, 129 if 0.025
  zgrid1d = [dz1*findgen(nz1) - 2.5, dz2*findgen(nz2) - 1.3,dz1*findgen(nz1 + 1) + 1.3 ];z_step_size*findgen(n_z) + z_min

  z_half_step_sizegrid1d =  [ replicate(dz1/2.,nz1) , replicate(dz2/2.,nz2) , replicate(dz1/2.,nz1 +1) ]
  x_half_step_sizegrid1d =  [ replicate(dx1/2.,nx1) , replicate(dx2/2.,nx2) , replicate(dx1/2.,nx3) , replicate(dx2/2.,nx2) , replicate(dx1/2.,nx1+1) ]
  ygrid1d = xgrid1d
  y_half_step_sizegrid1d  = x_half_step_sizegrid1d

  i_z0 = 64
  i_x0 = 139
  i_y0 = i_x0
  ;;;;
 
  xgrid3d = xgrid

  ygrid3d = ygrid

  zgrid3d=zgrid 
 
  xgrid = xgrid1d
  ygrid = ygrid1d
  zgrid=zgrid1d
  
  nel = [dindgen(9) + 1d, 10d*dindgen(9) + 10d, 100d*dindgen(4) + 100d, 250d*dindgen(39) + 500d]

  Tec = [0.1d*dindgen(9) + 0.1d, 0.5d*dindgen(17)+1d,dindgen(10)+10d,5d*dindgen(8) + 20d, 20d*dindgen(7) + 60d, 80d*dindgen(6) + 200d]


  nelhalf = [replicate(0.5d,9), replicate(5d,9), replicate(50d,4), replicate(125d,39)]

  Techalf = [replicate(0.05d,9), replicate(0.25d,17),replicate(0.5d,10),replicate(2.5d,8), replicate(10d,7), replicate(40d,6)]


  ypts_vary_x0_and_z0=dblarr(n_elements(xgrid),n_elements(zgrid), 8)
  ypts_total_vary_x0_and_z0=dblarr(n_elements(xgrid),n_elements(zgrid))

  openr,1,'ypts_total_vary_x0_and_z0_emission_table_281x129.txt'
  ypts_total_vary_x0_and_z0=dblarr(n_elements(xgrid),n_elements(zgrid)) ; ne x tec x discrete wavelength centers of each emission line
  readf,1,ypts_total_vary_x0_and_z0
  close,1
  
  openr,1,'ypts_vary_x0_and_z0_emission_table_281x129x884.txt'
  ypts_vary_x0_and_z0=dblarr(n_elements(xgrid),n_elements(zgrid),884) ; ne x tec x discrete wavelength centers of each emission line
  readf,1,ypts_vary_x0_and_z0
  close,1


idx = where(abs(ypts_total_vary_x0_and_z0) gt 10^6 )
ypts_total_vary_x0_and_z0(idx) = 0.0
print,idx

xgrid2d = fltarr(n_elements(xgrid),n_elements(zgrid))
zgrid2d = fltarr(n_elements(xgrid),n_elements(zgrid))


for i=0,n_elements(xgrid)-1 do begin
  for j=0,n_elements(zgrid)-1 do begin
    xgrid2d(i,j) = xgrid(i) 
    zgrid2d(i,j) = zgrid(j)
  endfor
endfor
print,xgrid2d(idx)
print,zgrid2d(idx)


p1 = jade_spectrogram(ypts_total_vary_x0_and_z0,xgrid,x_half_step_sizegrid1d,zgrid,z_half_step_sizegrid1d)



stop

end