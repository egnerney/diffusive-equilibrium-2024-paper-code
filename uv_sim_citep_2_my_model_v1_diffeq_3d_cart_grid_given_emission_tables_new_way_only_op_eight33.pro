
function simulate_IPT_spectrum_Rayleighs_citep2_ERF_form,x,spec_binsize, xwavi,yptsi, $
  fwhm = fwhm
  ;Tel in eV, nel in #/cm^3, Nsp etc in #/cm^3, max/min/fwhm in angstroms
  ; sample spectra

  ;if keyword_set(fwhm) then fwhm = fwhm else fwhm = 4.d

  ; relation between sigma and fwhm:
  ; => fwhm = 2*sigma*sqrt(2*ln(2))
  ;sigma = fwhm/(2d*sqrt(2d*alog(2d)))
  ;sigma2 = sigma*sigma
  ; a = 1d/(sigma*sqrt(2d*!dpi)) ; normalization
  ;c=1d/(2d*sigma2);
  rootc= 2d*sqrt(alog(2d))/fwhm;=1d/(sqrt(2d)*sigma)


  maxwav = max(x)
  minwav = min(x)

  nj=n_elements(x)


  ni=n_elements(yptsi)

  ypts=dblarr(nj)
  yptsij=dblarr(ni,nj)




  ;for j=0, nj-1 do begin
  for i=0, ni-1 do begin ; going through each discrete emission wavelength center




    ypts += yptsi(i)*0.5d*(Erf((x - xwavi(i) + spec_binsize/2d)*rootc) - Erf((x - xwavi(i) - spec_binsize/2d)*rootc))




  endfor

  ypts /= spec_binsize
  ;stop
  ;ypts(j)=total(yptsij(*,j))

  ;endfor


  return, ypts


end


;;;;;;;;;;;spacer

FUNCTION find_closest, grid, value
  differences = ABS(grid - value)
  min_of_differences = MIN(differences, Min_Subscript)
  RETURN, grid[Min_Subscript]
END


;;;;;;;;;;;spacer
function calc_LOS_n_s_and_Tel_s_given_my_model_3D_cart_grid_diffthreee, xgrid, ygrid, zgrid, norm_vec, slit_pos_vec, n_in, T_in, s=s, tel=tel
  ;, xgrid,ygrid,zgrid,norm_vec,slit_pos_vec,n_in,t_in,L_in,s=s,tel=tel
  epsilon = 1d-3



  ; Compute s values for each direction based on the provided norm_vec components.
  ; Note: If the component of norm_vec is close to zero, the resulting s value for that direction
  ; will be large and out of bounds for the rest of the calculations. We'll filter these out later.

  ds = 0.025d

  num_s =  floor(Sqrt(20.^2. + 20.^2. + 5.^2.) / ds) + 1
  s = ds*dindgen(num_s)

  ;calculate x,y, & zvals along line of sight (LOS)

  x_los_higher_res = slit_pos_vec[0] + s*norm_vec[0]; in RJ
  y_los_higher_res = slit_pos_vec[1] + s*norm_vec[1]; in RJ
  z_los_higher_res = slit_pos_vec[2] + s*norm_vec[2]; in RJ

  idx_want  = WHERE(x_los_higher_res GE MIN(xgrid) AND x_los_higher_res LE MAX(xgrid) AND y_los_higher_res GE MIN(ygrid) AND y_los_higher_res LE MAX(ygrid) AND z_los_higher_res GE MIN(zgrid) AND z_los_higher_res LE MAX(zgrid))
  x_los_higher_res = x_los_higher_res[idx_want]
  y_los_higher_res = y_los_higher_res[idx_want]
  z_los_higher_res = z_los_higher_res[idx_want]
  s_los = s[idx_want]

  num_s = n_elements(x_los_higher_res)


  x_los = dblARR(num_s)
  y_los = dblARR(num_s)
  z_los = dblARR(num_s)



  FOR i=0, num_s - 1 DO BEGIN
    x_los[i] = find_closest(xgrid, x_los_higher_res[i])
    y_los[i] = find_closest(ygrid, y_los_higher_res[i])
    z_los[i] = find_closest(zgrid, z_los_higher_res[i])
  ENDFOR

  ; 3. Filter out duplicates
  unique_points = HASH()
  x_los_unique = []
  y_los_unique = []
  z_los_unique = []
  s_los_unique = []

  FOR i=0, num_s-1 DO BEGIN
    key = STRING(x_los[i]) + '_' + STRING(y_los[i]) + '_' + STRING(z_los[i])
    IF ~unique_points.HASKEY(key) THEN BEGIN
      unique_points[key] = 1
      x_los_unique = [x_los_unique, x_los[i]]
      y_los_unique = [y_los_unique, y_los[i]]
      z_los_unique = [z_los_unique, z_los[i]]
      s_los_unique = [s_los_unique, s_los[i]]
    ENDIF
  ENDFOR

  s = s_los_unique


  ;finite_indices = WHERE(FINITE(s_los_unique) eq 1)

  ;s = s_los_unique;(finite_indices)
  ;print,N_elements(s), n_elements(x_los_unique)
  ;stop

  ;PRINT, 'NUnique x:',n_elements(x_los_unique),'min(x):',min(x_los_unique),'max(x):',max(x_los_unique),'Unique x:', x_los_unique
  ;PRINT, 'NUnique y:',n_elements(y_los_unique),'min(y):',min(y_los_unique),'max(y):',max(y_los_unique),'Unique y:', y_los_unique
  ;PRINT, 'NUnique z:',n_elements(z_los_unique),'min(z):',min(z_los_unique),'max(x):',max(z_los_unique),'Unique z:', z_los_unique
  ;stop



  x_los = x_los_unique
  y_los = y_los_unique
  z_los = z_los_unique
  s = s_los_unique

  n_a=n_elements(x_los)
  rho_los = sqrt(x_los^2 + y_los^2) ; in RJ
  phi_los = (180d/!dpi)*atan(y_los,x_los)
  ;stop
  ;print,x_los
  ;stop
  ;print, y_los
  ;stop
  ;
  ;rho_los = sqrt(x_los^2 + y_los^2) ; in RJ
  ;r_los = sqrt(x_los^2 + y_los^2 + z_los^2) ; in RJ
  ;lat_los = (180d/!dpi)*asin(z_los/r_los) ; latitude in degrees systemIII

  ;phi_los = (180d/!dpi)*atan(y_los,x_los); ELONG RIGHTHANDED in degrees systemIII
  ;print,rho_los

  ;print,size(rho_los,/n_dimen)
  ;help,rho_los


  x_index=make_array(n_a,/integer)
  y_index=make_array(n_a,/integer)
  z_index=make_array(n_a,/integer)

  ;T=replicate({Temps_LOS, s: 0.0d, sp: 0.0d, s2p: 0.0d, s3p: 0.0d, s4p:0.0d, $
  ; o: 0.0d, op: 0.0d, o2p: 0.0d, el: 0.0d,elh:0.0d},n_a)

  n=replicate({densities_LOS, s: 0.0d, sp: 0.0d, s2p: 0.0d, s3p: 0.0d, s4p:0.0d, $
    o: 0.0d, op: 0.0d, o2p: 0.0d, el: 0.0d,elh:0.0d},n_a)




  tel = dblarr(n_a)

  epsilon = 1d-5 
  for a = 0, n_a - 1  do begin

    ii = where((xgrid le (x_los(a) + epsilon)) and (xgrid ge (x_los(a) - epsilon)))



    jj = where((ygrid le (y_los(a) + epsilon)) and (ygrid ge (y_los(a) - epsilon)))

    kk = where((zgrid le (z_los(a) + epsilon)) and (zgrid ge (z_los(a) - epsilon)))

    n[a].sp = n_in[ii,jj,kk].sp;griddata(  rho_in, z_in,n_in.sp,xout = rho_los, yout = z_los )
    n[a].s2p = n_in[ii,jj,kk].s2p; griddata( rho_in, z_in,n_in.s2p, xout = rho_los, yout = z_los )
    n[a].s3p = n_in[ii,jj,kk].s3p  ;griddata(rho_in, z_in,n_in.s3p,  xout = rho_los, yout = z_los )
    n[a].op = n_in[ii,jj,kk].op ;griddata( rho_in, z_in,n_in.op,  xout = rho_los, yout = z_los )
    n[a].o2p = n_in[ii,jj,kk].o2p ; griddata( rho_in, z_in,n_in.o2p, xout = rho_los, yout = z_los )
    n[a].el = n_in[ii,jj,kk].el  ;griddata(rho_in, z_in, n_in.el,  xout = rho_los, yout = z_los )
    tel[a] =  T_in[ii,jj,kk].el;griddata(rho_in, z_in, T_in.el,  xout = rho_los, yout = z_los )
  endfor
  ;print, n.el
  ;help, n.el
  ;stop
  ;idx2 = where;where((l_shell_los lt 5.) or (l_shell_los gt 10.) ) ; approximate dipole implied inner and outer edges of fits where everything should be 0 but interp might have other ideas!
  ;n[idx2].el= 0.0
  ; n[idx2].sp= 0.0
  ;  n[idx2].s2p= 0.0
  ; n[idx2].s3p= 0.0
  ;n[idx2].op= 0.0
  ;n[idx2].o2p= 0.0
  ;Tel[idx2] = 0.0
  ;stop
  return,n

end





;;;;;;;;;;spacer








Function new_CITEP_2_given_n_s_and_tel_s_using_emission_tables_sum_eightthirtythree_op_line,s,n,Tel_s,fwhm,spec_binsize,yptsi_,nel_,tec_


n_a=n_elements(tel_s)
B_intgrnd_values=dblarr(n_a)
zero_idxs=where((Tel_s eq 0d) or (n.el eq 0d) )
neq_zero_idxs=where((Tel_s ne 0d) and (n.el ne 0d))
  nel_s=reform(n.el)
   nsp_s=reform(n.sp)
    ns2p_s=reform(n.s2p)
     ns3p_s=reform(n.s3p)
      nop_s=reform(n.op)
       no2p_s=reform(n.o2p)
       
      
 
;stop
if n_elements(zero_idxs) ne n_elements(nel_s) then begin
  
  
    remove, zero_idxs, Tel_s, nel_s, nsp_s, ns2p_s, ns3p_s, nop_s, no2p_s
  
  

yptsi_los = dblarr(n_elements(nel_s))
    ;yptsi_,nel_,tec_,xwavi,yname
 ;for j=0,num_emiss-1 do begin
  
  triangulate,nel_,tec_,tr
  yptsi_LOS = griddata(  nel_, tec_,yptsi_,xout = nel_s, yout = tel_s,/linear,triangles = tr  )
 
  
  yptsi_LOS = nop_s*yptsi_LOS
   
 
 
;stop

  if (zero_idxs(0) ge 0) then for i=0,n_elements(zero_idxs)-1 do B_intgrnd_values(zero_idxs(i))=0d
  if (neq_zero_idxs(0) ge 0) then for i=0,n_elements(neq_zero_idxs)-1 do B_intgrnd_values(neq_zero_idxs(i))=yptsi_LOS(i)



 ; yptsi=dblarr(num_emiss)
  ;for j=0,num_emiss-1 do begin


    ypts=7.1492d3*INT_TABULATED(s, B_intgrnd_values,/double) ; no factor of 2 this time

  ;endfor

 ; ypts=simulate_IPT_spectrum_Rayleighs_citep2_ERF_form(x_wav,spec_binsize, xwavi,yptsi, $
 ;   fwhm = fwhm)

  
endif else begin
  ypts = 0d
endelse
 

  return, ypts

End


Pro UV_sim_citep_2_my_model_v1_diffeq_3D_cart_grid_given_emission_tables_new_way_only_op_eight33

   openr,1,'yptsi_550-2100_Angstroms_vary_ne_te_lookuptable_61x57x884_CHIANTI801.dat'
   yptsi_in=dblarr(61,57,884) ; ne x tec x discrete wavelength centers of each emission line  
   readf,1,yptsi_in
   close,1
   

   
   openr,1,'xwavi_550-2100_Angstroms_CHIANTI801.dat'
   xwavi_in=dblarr(884)
   readf,1,xwavi_in
   close,1
   

   
   openr,1,'yname_550-2100_Angstroms_CHIANTI801.dat'
   yname_in=strarr(884)
   readf,1,yname_in
   close,1
   
  sp_idxs = where(yname_in eq '"S  II"')
  s2p_idxs = where(yname_in eq '"S  III"')
  s3p_idxs = where(yname_in eq '"S  IV"')
  s4p_idxs = where(yname_in eq '"S  V"')
  op_idxs = where(yname_in eq '"O  II"')
  o2p_idxs = where(yname_in eq '"O  III"')
  
  
  idxxx = where((xwavi_in(op_idxs)  le 843) and (xwavi_in(op_idxs)  ge 829) )
 xwavi_in =xwavi_in(op_idxs)
 
 yptsi_in_temp = reform(yptsi_in(28,17,op_idxs))
  print,xwavi_in(idxxx)
  print, yptsi_in_temp(idxxx)
  
 
  yptsi_in_op_eightthirtythree = dblarr(61,57)
 stop
  for iiii = 0, 61 -1 do begin
  for jjjj = 0, 57 -1  do begin ; to make single number to simulate total brightness because linearity of integrals
    
    yptsi_in_temp = reform(yptsi_in(iiii,jjjj,op_idxs))
    
    
    yptsi_in_op_eightthirtythree(iiii,jjjj)  = total(yptsi_in_temp(idxxx),/double)
    
    

  endfor
  endfor
  
  yptsi_in = yptsi_in_s2p_sixeightyline
  
  
  
  



 

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
 ; xgrid3d = xgrid 
  
  ;ygrid3d = ygrid 
  
  ;zgrid3d = zgrid 
  
  xgrid = xgrid1d
  ygrid = ygrid1d
  zgrid = zgrid1d
  
  
  n_in=replicate({densities, s: 0.0, sp: 0.0, s2p: 0.0, s3p: 0.0, s4p:0.0, $
    o: 0.0, op: 0.0,oph: 0.0, o2p: 0.0, el: 0.0,elh:0.0},n_elements(xgrid),n_elements(ygrid),n_elements(zgrid)) ;(cm^-3) 501 L shells 5-10RJ 0.01RJ res 1001 lat -50-50 degrees 0.1 degree res

  T_in=replicate({Temperatures, s: 0.0, sp: 0.0, s2p: 0.0, s3p: 0.0, s4p:0.0, $
    o: 0.0, op: 0.0,oph: 0.0, o2p: 0.0, el: 0.0,elh:0.0},n_elements(xgrid),n_elements(ygrid),n_elements(zgrid)) ;(eV), constant over field lines for maxwellians as this case is
    
  Restore, 'Mymodelv1_3D_JRM33+con20_cartesian_irrcart_mapped_no_interp_2degree_phi_into_rhoIII_4.5_279x279x129_with_hote.sav',/verbose
 ; Restore, 'Mymodelv1_3D_JRM33+con20_cartesian_irrcart_mapped_no_interp_2degree_phi_into_rhoIII_4.5_279x279x129.sav',/verbose


;T_in.sp = Ti_in
T_in.el = Tec_3Dcart
n_in.el = nel_3Dcart

n_in.sp = nsp_3Dcart
n_in.s2p = ns2p_3Dcart
n_in.s3p = ns3p_3Dcart
n_in.op = nop_3Dcart
n_in.o2p = no2p_3Dcart


 ;stop
xmin=550d;3000d;
xmax=2100d;7000d;
xstepsize=1d
nx=Round((xmax-xmin)/xstepsize) + 1d
x=dindgen(nx)*xstepsize + xmin 
;print,x
  

  
;euvbin=0.6049d ;euv binsize cassini UVIS
  binsize=1d;0.6049d ; cassini res
  fwhm=4.47;24.5d ;  UVIS res
  
  ;norm_vec=[0d,1d/sqrt(2d),-1d/sqrt(2d)] ; looking at 45 degree angle do
  norm_vec=[0d,1d,0d] ;looking in yhat direction
  x0=6d
  y0=-10d
  z0=0d
  slit_pos_vec=[x0,y0,z0]
  
 
  
  
  ;imin=uint(reform(min(where(xgrid eq 4.5 )))) ;100 for current res.
  ;imax=uint(reform(min(where(xgrid eq 10.0 )))) 
 ; jmin=uint(reform(min(where(zgrid eq 0.0 ))))
  ;jmax=uint(reform(min(where(zgrid eq 3.0 ))))
  
  
 
   ; ypts_total_vary_x0_and_z0=dblarr(nframes,n_elements(xgrid),n_elements(zgrid))

    nel = [dindgen(9) + 1d, 10d*dindgen(9) + 10d, 100d*dindgen(4) + 100d, 250d*dindgen(39) + 500d]

    Tec = [0.1d*dindgen(9) + 0.1d, 0.5d*dindgen(17)+1d,dindgen(10)+10d,5d*dindgen(8) + 20d, 20d*dindgen(7) + 60d, 80d*dindgen(6) + 200d]

    n_ne = n_elements(nel)

    n_te = n_elements(tec)
    
    nel2dto1d =dblarr(n_ne*n_te)
    tec2dto1d =dblarr(n_ne*n_te)
    yptsi_in_2dto1d = dblarr(n_ne*n_te)

;2nd ->1d mapping for griddata interp
k=-1
for i=0, 61 -1 do begin
  for j=0, 57 -1 do begin
    k = k + 1
    nel2dto1d[k] = nel[i]
    tec2dto1d[k] = tec[j]
    
    yptsi_in_2dto1d[k] = yptsi_in[i,j]



  endfor
endfor


nframes = 45;8d ; proof of concept for now...
ypts_vary_x0_and_z0=dblarr(nframes,n_elements(xgrid),n_elements(zgrid))



dphi = (360d/nframes)*(!dpi/180d)

phi = dphi*dindgen(nframes)



tic


;for k=0, nframes -1 do begin
k=0
  print,k
  print,phi(k)*(180d/!dpi)
  norm_vec=[-sin(phi[k]) , Cos(phi[k]) , 0d]
For i = 0, n_elements(xgrid)-1 do begin ; 
  For j = 0, n_elements(zgrid)-1 do begin ; 
    
    slit_pos_vec = [xgrid(i)*Cos(phi[k]) + 10d*Sin(phi[k]) , xgrid(i)*sin(phi[k]) - 10d*cos(phi[k]) , zgrid(j)];[x0t,y0t,z0];[x0,y0,z0]
    n = calc_LOS_n_s_and_Tel_s_given_my_model_3D_cart_grid_diffthreee( xgrid, ygrid, zgrid, norm_vec, slit_pos_vec, n_in, T_in, s=s, tel=tel)
 
    ypts=new_CITEP_2_given_n_s_and_tel_s_using_emission_tables_sum_sixeighty_s2p_line(s,n,Tel,fwhm,binsize,yptsi_in_2dto1d,nel2dto1d,tec2dto1d)
   
    ypts_vary_x0_and_z0(k,i,j)=ypts
   


  Endfor
Endfor

;write_csv,'ypts_vary_x0_and_z0_emission_table_45x279x129_new_way_UV_only680s2plines_fixed.txt',ypts_vary_x0_and_z0

;endfor



;p1=jade_spectrogram(reform(ypts_total_vary_x0_and_z0(k,*,*)),xgrid,x_half_step_sizegrid1d,zgrid,z_half_step_sizegrid1d)
;stop
p1=jade_spectrogram(reform(ypts_vary_x0_and_z0(k,*,*)),xgrid,x_half_step_sizegrid1d,zgrid,z_half_step_sizegrid1d)

;stop

;write_csv,'ypts_vary_x0_and_z0_minusyhat_from_y=-10_mymodelv1_irr_cart_grid_277x277x101_EUV.txt',ypts_vary_x0_and_z0
;write_csv,'ypts_total_vary_x0_and_z0_minusyhat_from_y=-10_mymodelv1_irr_cart_grid_277x277x101.txt',ypts_total_vary_x0_and_z0
  ;;;;
  toc
  
  ;n=calc_LOS_n_s_and_Tel_s_azimuthally_sym_scale_heights( xgrid,ygrid,zgrid,norm_vec,slit_pos_vec,an,bn,atel,btel,s=s,tel=tel)
  ;s,n,Tel_s,x_wav,fwhm,spec_binsize
  ;toc
; xrange1=[xmin,xmax]
  
  ;p1=plot( x_wav_uvis,y_uvis,xtitle='Wavelength (Å)',ytitle='Rayleighs/Å',title='', NAME='Total',xrange=xrange1)
  ;p2=plot( x_wav_uvis,3d*ypts,xtitle='Wavelength (Å)',ytitle='Rayleighs/Å',title='',color='red',/overplot,transparency=50,xrange=xrange1)
  
  
    ;write_csv,'ypts_total_vary_x0_and_z0_emission_table_279x129_new_way_UV.txt',ypts_total_vary_x0_and_z0
  toc
  stop
  




end




