pro read_in_systemIII_IPT_model_for_Todd
  Restore, 'Mymodelv1_3D_JRM33+con20_cartesian_irrcart_mapped_no_interp_2degree_phi_into_rhoIII_4.5_281x281x129.sav',/verbose
  
  i_to_keep = [indgen(161), indgen(118)+163] ; picking out the 279 elements in x and y to keep



  nel_3Dcart = nel_3Dcart(i_to_keep,*,*)
  nel_3Dcart = nel_3Dcart(*,i_to_keep,*)

  nsp_3Dcart = nsp_3Dcart(i_to_keep,*,*)
  nsp_3Dcart = nsp_3Dcart(*,i_to_keep,*)
  ns2p_3Dcart = ns2p_3Dcart(i_to_keep,*,*)
  ns2p_3Dcart = ns2p_3Dcart(*,i_to_keep,*)
  ns3p_3Dcart = ns3p_3Dcart(i_to_keep,*,*)
  ns3p_3Dcart = ns3p_3Dcart(*,i_to_keep,*)
  nop_3Dcart = nop_3Dcart(i_to_keep,*,*)
  nop_3Dcart = nop_3Dcart(*,i_to_keep,*)
  no2p_3Dcart = no2p_3Dcart(i_to_keep,*,*)
  no2p_3Dcart = no2p_3Dcart(*,i_to_keep,*)
  
  noph_3Dcart = noph_3Dcart(i_to_keep,*,*)
  noph_3Dcart = noph_3Dcart(*,i_to_keep,*)
  nhp_3Dcart = nhp_3Dcart(i_to_keep,*,*)
  nhp_3Dcart = nhp_3Dcart(*,i_to_keep,*)
  nnap_3Dcart = nnap_3Dcart(i_to_keep,*,*)
  nnap_3Dcart = nnap_3Dcart(*,i_to_keep,*)
  
  Ti_3Dcart = Ti_3Dcart(i_to_keep,*,*)
  Ti_3Dcart = Ti_3Dcart(*,i_to_keep,*)
  Tec_3Dcart = Tec_3Dcart(i_to_keep,*,*)
  Tec_3Dcart = Tec_3Dcart(*,i_to_keep,*)
  Thp_3Dcart = Thp_3Dcart(i_to_keep,*,*)
  Thp_3Dcart = Thp_3Dcart(*,i_to_keep,*)
  Toph_3Dcart = Toph_3Dcart(i_to_keep,*,*)
  Toph_3Dcart = Toph_3Dcart(*,i_to_keep,*)
  
 
  
  
  dx1 = 0.1
  nx1 = 40

  dx2 = 0.05;0.025
  nx2 = 78;390 is for 0.01 ; 156 is 0.025 ; 78 for 0.05 ;
  ;429 for x and y for 0.025, 281 for 0.05
  nx3 = 42 ; middle bits at 0.1, was 44 before with mistake


  xgrid1d = [ dx1*findgen(nx1) - 10. , dx2*findgen(nx2) - 6. , dx1*findgen(nx3) - 2.1 , dx2*findgen(nx2) + 2.1 , dx1*findgen(nx1+1) + 6. ]



  dz1 = 0.1
  nz1 = 12

  dz2 = 0.025
  nz2 = 104;104 for 0.025 for -1.3 to 1.3 so 2.6RJ total
  
  zgrid1d = [dz1*findgen(nz1) - 2.5, dz2*findgen(nz2) - 1.3,dz1*findgen(nz1 + 1) + 1.3 ];z_step_size*findgen(n_z) + z_min

  z_half_step_sizegrid1d =  [ replicate(dz1/2.,nz1) , replicate(dz2/2.,nz2) , replicate(dz1/2.,nz1 +1) ]
  x_half_step_sizegrid1d =  [ replicate(dx1/2.,nx1) , replicate(dx2/2.,nx2) , replicate(dx1/2.,nx3) , replicate(dx2/2.,nx2) , replicate(dx1/2.,nx1+1) ]
  ygrid1d = xgrid1d
  y_half_step_sizegrid1d  = x_half_step_sizegrid1d

  i_z0 = 64 ; z =0 index
  i_x0 = 139 ; x=0 index
  i_y0 = i_x0 ; y=0 index
  ;;;;
  xgrid3d = xgrid ; 3D version with x value everywhere in carteisan 3D grid 281x281x129
 
  ygrid3d = ygrid ; 3D version with y value everywhere in carteisan 3D grid 281x281x129
 
  zgrid3d = zgrid  ; 3D version with z value everywhere in carteisan 3D grid 281x281x129
  
  xgrid = double(xgrid1d)
  ygrid = double(ygrid1d)
  zgrid=double(zgrid1d)

;p1=jade_spectrogram(reform(NEL_3DCART(*,i_y0,*)),xgrid,x_half_step_sizegrid1d,zgrid,z_half_step_sizegrid1d,xtitle='x ($R_J$)',ytitle='z ($R_J$)',ctitle='$n_e$ $(cm^{-3})$')
;stop
 ; p2=jade_spectrogram(reform(tec_3DCART(*,i_y0,*)),xgrid,x_half_step_sizegrid1d,zgrid,z_half_step_sizegrid1d,xtitle='x ($R_J$)',ytitle='z ($R_J$)',ctitle='$T_{ec}$ $(eV)$')

;
;n_elements(xgrid1d) = 281 ; irregular cartesian grid from -10 to 10 RJ for systemIII x
;n_elements(ygrid1d) = 281 ; irregular cartesian grid from -10 to 10 RJ for systemIII y
;n_elements(zgrid1d) = 129 ; ; irregular cartesian grid from -2.5 to 2.5 RJ in systemIII z 
;n_elements(NEL_3DCART) = 281x281x129 (3D array) = 10,185,969 elements, everything inside of rho_III = 4.5 is set to 0
;0 fills for density and temperatuyres used everywhere no values so be careful of dividing by 0 outside of torus
;NEL_3DCART(*,i_y0,*) is electron density slice in x-z systemIII plane at y=0
 ;same for, NSP_3DCART;NS2P_3DCART;NS3P_3DCART;NOP_3DCART;NO2P_3DCART;NOPH_3DCART;NHP_3DCART;NNAP_3DCART;TEC_3DCART;TI_3DCART;TOPH_3DCART;THP_3DCART








 ;SAVE,xgrid,ygrid,zgrid,nel_3dcart,nsp_3dcart,ns2p_3dcart,ns3p_3dcart,nop_3dcart,no2p_3dcart,noph_3dcart,nhp_3dcart,nnap_3dcart,Tec_3dcart,Ti_3dcart,Toph_3dcart,Thp_3dcart, FILENAME = 'Mymodelv1_3D_JRM33+con20_cartesian_irrcart_mapped_no_interp_2degree_phi_into_rhoIII_4.5_279x279x129.sav',/verbose


 nxgrid = long64(n_elements(xgrid))
 nygrid = long64(n_elements(ygrid))
 nzgrid = long64(n_elements(zgrid))
 
 xgrid_for_3dto2d = dblarr( nxgrid*nygrid)
 ygrid_for_3dto2d = xgrid_for_3dto2d

 nel_3dto2d = dblarr( nxgrid*nygrid, nzgrid)
 nsp_3dto2d = nel_3dto2d
 ns2p_3dto2d = nel_3dto2d
 ns3p_3dto2d = nel_3dto2d
 nop_3dto2d = nel_3dto2d
 no2p_3dto2d = nel_3dto2d
 
 noph_3dto2d = nel_3dto2d
 nhp_3dto2d = nel_3dto2d
 nnap_3dto2d = nel_3dto2d
 
  Ti_3dto2d = nel_3dto2d
  Tec_3dto2d = nel_3dto2d
  Toph_3dto2d = nel_3dto2d
  Thp_3dto2d = nel_3dto2d
 

 ;2nd ->1d mapping for griddata interp
 k=long64(-1)
 for i=0, nxgrid -1 do begin
   for j=0, nygrid -1 do begin
     k = k + long64(1)
     xgrid_for_3dto2d[k] = xgrid[i]
     ygrid_for_3dto2d[k] = ygrid[j]

      nel_3dto2d[k,*] = nel_3dcart[i,j,*]
      nsp_3dto2d[k,*] = nsp_3dcart[i,j,*]
      ns2p_3dto2d[k,*] = ns2p_3dcart[i,j,*]
      ns3p_3dto2d[k,*] = ns3p_3dcart[i,j,*]
      nop_3dto2d[k,*] = nop_3dcart[i,j,*]
      no2p_3dto2d[k,*] = no2p_3dcart[i,j,*]
      noph_3dto2d[k,*] = noph_3dcart[i,j,*]
      nhp_3dto2d[k,*] = nhp_3dcart[i,j,*]
      nnap_3dto2d[k,*] = nnap_3dcart[i,j,*]
      
      Tec_3dto2d[k,*] =Tec_3dcart[i,j,*]
      Ti_3dto2d[k,*] =Ti_3dcart[i,j,*]
      Toph_3dto2d[k,*] =Toph_3dcart[i,j,*]
      Thp_3dto2d[k,*] =Thp_3dcart[i,j,*]
      

     ;print,yptsi_in_2dto1d[k,*]
     ;stop
   endfor
 endfor
 
 
 write_csv,'nel_3dto2d_mapping_77841x129.csv',nel_3dto2d
 write_csv,'nsp_3dto2d_mapping_77841x129.csv',nsp_3dto2d
 write_csv,'ns2p_3dto2d_mapping_77841x129.csv',ns2p_3dto2d
 write_csv,'ns3p_3dto2d_mapping_77841x129.csv',ns3p_3dto2d
 write_csv,'nop_3dto2d_mapping_77841x129.csv',nop_3dto2d
 write_csv,'no2p_3dto2d_mapping_77841x129.csv',no2p_3dto2d
 write_csv,'noph_3dto2d_mapping_77841x129.csv',noph_3dto2d
 write_csv,'nhp_3dto2d_mapping_77841x129.csv',nhp_3dto2d
 write_csv,'nnap_3dto2d_mapping_77841x129.csv',nnap_3dto2d
 write_csv,'Tec_3dto2d_mapping_77841x129.csv',Tec_3dto2d
 write_csv,'Ti_3dto2d_mapping_77841x129.csv',Ti_3dto2d
 write_csv,'Thp_3dto2d_mapping_77841x129.csv',Thp_3dto2d
 write_csv,'Toph_3dto2d_mapping_77841x129.csv',Toph_3dto2d
 
  write_csv,'xgrid_for_3dto2d_mapping_77841.csv',xgrid_for_3dto2d
 write_csv,'ygrid_for_3dto2d_mapping_77841.csv',ygrid_for_3dto2d
 
 
 
 






stop
end