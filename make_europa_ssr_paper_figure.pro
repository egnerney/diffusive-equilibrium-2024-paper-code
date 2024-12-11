
function simulate_IPT_spectrum_Rayleighs_citep2_ERF_form_,x,spec_binsize__, xwavi,yptsi, fwhm = fwhm
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


;print,spec_binsize__
;stop

  ;for j=0, nj-1 do begin
  for i=0, ni-1 do begin ; going through each discrete emission wavelength center




    ypts += yptsi(i)*0.5d*(Erf((x - xwavi(i) + spec_binsize__/2d)*rootc) - Erf((x - xwavi(i) - spec_binsize__/2d)*rootc))




  endfor

  ypts /= spec_binsize__
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
function calc_LOS_n_s_and_Tel_s_given_my_model_3D_cart_grid_diffthreee_for_faraway, xgrid, ygrid, zgrid, norm_vec, slit_pos_vec, n_in, T_in, s=s, tel=tel
  ;, xgrid,ygrid,zgrid,norm_vec,slit_pos_vec,n_in,t_in,L_in,s=s,tel=tel
  epsilon = 1d-3



  ; Compute s values for each direction based on the provided norm_vec components.
  ; Note: If the component of norm_vec is close to zero, the resulting s value for that direction
  ; will be large and out of bounds for the rest of the calculations. We'll filter these out later.

  ds = 0.025d

  num_s =  floor(Sqrt(40.^2. + 40.^2. + 5.^2.) / ds) + 1
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








Function new_CITEP_2_given_n_s_and_tel_s_using_emission_tables_uvs,s,n,Tel_s,fwhm,spec_binsize_,yptsi_,nel_,tec_,ynamei,x_wave, xwavi

nlines = 4153
n_LOS=n_elements(tel_s)
B_intgrnd_values=dblarr(n_LOS,nlines )
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
  
  

yptsi_los = dblarr(n_elements(neq_zero_idxs),nlines);dblarr(n_LOS,nlines)
    ;yptsi_,nel_,tec_,xwavi,yname
 ;for j=0,num_emiss-1 do begin
   triangulate,nel_,tec_,tr
  for lll = 0 , nlines - 1 do begin
  CASE ynamei(lll) OF
    '"S II"': yptsi_LOS[*,lll] = nsp_s*griddata(  nel_, tec_,yptsi_[*,lll],xout = nel_s, yout = tel_s,/linear,triangles = tr  )
    '"S III"': yptsi_LOS[*,lll] = ns2p_s*griddata(  nel_, tec_,yptsi_[*,lll],xout = nel_s, yout = tel_s,/linear,triangles = tr  )
    '"S IV"': yptsi_LOS[*,lll] = ns3p_s*griddata(  nel_, tec_,yptsi_[*,lll],xout = nel_s, yout = tel_s,/linear,triangles = tr  )
    '"S V"': yptsi_LOS[*,lll] = replicate(0d,n_elements(neq_zero_idxs))
    '"O II"': yptsi_LOS[*,lll] = nop_s*griddata(  nel_, tec_,yptsi_[*,lll],xout = nel_s, yout = tel_s,/linear,triangles = tr  )
    '"O III"': yptsi_LOS[*,lll] = no2p_s*griddata(  nel_, tec_,yptsi_[*,lll],xout = nel_s, yout = tel_s,/linear,triangles = tr  )
    ELSE: stop
  ENDCASE
 endfor
 
 
;stop

  if (zero_idxs(0) ge 0) then for i=0,n_elements(zero_idxs)-1 do B_intgrnd_values(zero_idxs(i),*)= replicate(0d,nlines)
  if (neq_zero_idxs(0) ge 0) then for i=0,n_elements(neq_zero_idxs)-1 do B_intgrnd_values(neq_zero_idxs(i),*)=yptsi_LOS(i,*)



 ; yptsi=dblarr(num_emiss)
  ;for j=0,num_emiss-1 do begin
yptsi = dblarr(nlines)
for lll = 0 , nlines - 1 do begin
    yptsi(lll) = 7.1492d3*INT_TABULATED(s, reform(B_intgrnd_values(*,lll)),/double) ; no factor of 2 this time
endfor


 

;p1=plot(s(neq_zero_idxs),tel_s,layout=[2,4,1])
;p1=plot(s(neq_zero_idxs),nel_s,layout=[2,4,2],/current)
;p1=plot(s(neq_zero_idxs),nsp_s,layout=[2,4,3],/current)
;p1=plot(s(neq_zero_idxs),ns2p_s,layout=[2,4,4],/current)
;p1=plot(s(neq_zero_idxs),ns3p_s,layout=[2,4,5],/current)
;p1=plot(s(neq_zero_idxs),nop_s,layout=[2,4,6],/current)
;p1=plot(s(neq_zero_idxs),no2p_s,layout=[2,4,7],/current)
;p1=plot(xwavi,yptsi)
;stop

;print,spec_binsize_
;stop
ypts = simulate_IPT_spectrum_Rayleighs_citep2_ERF_form_(x_wave,spec_binsize_, xwavi,yptsi, fwhm = fwhm)
 
endif else begin
  ypts = replicate(0d,n_elements(x_wave))
endelse
 

  return, ypts

End


Pro make_europa_ssr_paper_figure


nlines = 4153 ; for CHIANTI 10.1 there are more lines than for 8
   openr,1,'yptsi_550-2100_Angstroms_vary_ne_te_lookuptable_61x57x4153_CHIANTI10_1_noprotons_NOIONREC_NO_RREC.dat'
   yptsi_in=dblarr(61,57,nlines ) ; ne x tec x discrete wavelength centers of each emission line  
   readf,1,yptsi_in
   close,1
   

   
   openr,1,'xwavi_550-2100_Angstroms_CHIANTI10_1_noprotons_NOIONREC_NO_RREC.dat'
   xwavi_in=dblarr(nlines)
   readf,1,xwavi_in
   close,1
   

   
   openr,1,'yname_550-2100_Angstroms_CHIANTI10_1_noprotons_NOIONREC_NO_RREC.dat'
   yname_in=strarr(nlines)
   readf,1,yname_in
   close,1
   
   
   ;openr,1,'y.txt'
   ;yyy=fltarr(1569)
   ;readf,1,yyy
  ; close,1
   
   
   ;openr,1,'x.txt'
   ;xxx=fltarr(1569)
   ;readf,1,xxx
   ;close,1
   
   sp_idxs = where(yname_in eq '"S II"')
   s2p_idxs = where(yname_in eq '"S III"')
   s3p_idxs = where(yname_in eq '"S IV"')
   s4p_idxs = where(yname_in eq '"S V"')
   op_idxs = where(yname_in eq '"O II"')
   o2p_idxs = where(yname_in eq '"O III"')

;help,sp_idxs

 ;  stop
  

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
;stop

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
xstepsize=1d ; 0.6049d  ;4.47d
nx=Round((xmax-xmin)/xstepsize) + 1d
x_wave=dindgen(nx)*xstepsize + xmin 
;print,x
  

  
;euvbin=0.6049d ;euv binsize cassini UVIS
  spec_binsize = xstepsize;0.6049d ; cassini res ; 1d
  fwhm = 6d; JUICE/CLIPPER UVS FWHM best case scenario (Angstrom)... ;4.47 ang was  UVIS res

  
  

    nel = [dindgen(9) + 1d, 10d*dindgen(9) + 10d, 100d*dindgen(4) + 100d, 250d*dindgen(39) + 500d]

    Tec = [0.1d*dindgen(9) + 0.1d, 0.5d*dindgen(17)+1d,dindgen(10)+10d,5d*dindgen(8) + 20d, 20d*dindgen(7) + 60d, 80d*dindgen(6) + 200d]

    n_ne = Long(n_elements(nel))

    n_te = Long(n_elements(tec))
    
    nel2dto1d =dblarr(n_ne*n_te)
    tec2dto1d =dblarr(n_ne*n_te)
    yptsi_in_2dto1d = dblarr(n_ne*n_te,nlines)

;2nd ->1d mapping for griddata interp
k=-1


for i=0, 61 -1 do begin
  for j=0, 57 -1 do begin
    k = k + 1
    nel2dto1d[k] = nel[i]
    tec2dto1d[k] = tec[j]
    
    yptsi_in_2dto1d[k,*] = yptsi_in[i,j,*]



  endfor
endfor



restore,'clipper_spice_airglow_boresight_pointing_num_intervals=45_num_et_per_interval=1801_each_duration=30min_dt=1s.sav',/verbose
;BS_pointing_vector_Clipper_UVS_AP,et_SC,utc_SC,x_SC,y_SC,z_SC
n_intervals = 45;8d ; proof of concept for now...

dt = 1d ; seconds
n_et_each_interval = floor(30d*60d/dt) + 1;1801 ;30 minutes for each interval at 1 second resolution is (30*60/0.01) + 1 = 1801, duration/dt + 1

ypts_out = dblarr(n_et_each_interval,n_intervals,n_elements(x_wave))
;n_intervals = n_elements(intervals) ; 45
;ap_boresight_sys3_out = dblarr(3,n_et_each_interval,n_intervals); 3d pointing vector x 1801 et values x n_intervals
;x_out = dblarr(n_et_each_interval,n_intervals)
;y_out = dblarr(n_et_each_interval,n_intervals)
;z_out = dblarr(n_et_each_interval,n_intervals)
;et_out = dblarr(n_et_each_interval,n_intervals)
;utc_out = strarr(n_et_each_interval,n_intervals)

indicies_for_jjj_proxy = [0,900,1800] ; Beginning of et_interval, middle, and end.

strings_for_et_intervals = [', Beginning of  30 min Observation',', Middle of 30 min Observation',', End of 30 min Observation' ]

utc_outz =strarr(45,3)
    norm_vec_outz = dblarr(3,45,3)
    slit_pos_vec_outz = dblarr(3,45,3)
s_outz = replicate(!Values.F_NAN,45,3,2272) 
x_outz = s_outz
y_outz = s_outz
z_outz = s_outz
tec_outz = s_outz 
nec_outz = tec_outz
nsp_outz = tec_outz
ns2p_outz = tec_outz
ns3p_outz = tec_outz
nop_outz = tec_outz
no2p_outz = tec_outz

restore,'Europa_UVS_EffArea_20220331.sav'
effective_area_clipper = interpol(effarea_ap,w*10d,x_wave)

nspec_bins = n_elements(x_wave)

strds_per_pixel = 0.3d*0.1d*((!dpi/180d)^2d)
factor = ((10d^6d)/(4d*!dpi))*strds_per_pixel*effective_area_clipper



phi = (!dpi/180d)*68.3d ; alignedish

wave_idx_want = where(x_wave le 2060d)

;for iii = 0, n_intervals - 1 do begin
;  for jjj_proxy = 0, 2 do begin;n_et_each_interval - 1 do begin
    
    
  norm_vec = [-sin(phi) , Cos(phi) , 0d]; [0d,1d,0d];
    slit_pos_vec =  [5.91d*Cos(phi) + 5.91d*Sin(phi) , 5.91d*sin(phi) - 10d*cos(phi) , 0d];[x_SC[jjj,iii],y_SC[jjj,iii],z_SC[jjj,iii]];[x0t,y0t,z0];[x0,y0,z0] ; [6d,-10d,0d];
    n = calc_LOS_n_s_and_Tel_s_given_my_model_3D_cart_grid_diffthreee_for_faraway( xgrid, ygrid, zgrid, norm_vec, slit_pos_vec, n_in, T_in, s=s, tel=tel)
    
    
    
    ypts = new_CITEP_2_given_n_s_and_tel_s_using_emission_tables_uvs(s,n,Tel,fwhm,spec_binsize,yptsi_in_2dto1d,nel2dto1d,tec2dto1d,yname_in,x_wave, xwavi_in)
   ; ypts_out(jjj,iii,*) = ypts
   ypts1 = ypts
  ; p2 = plot(xxx,yyy,xtitle='Wavelength (Å)',YTITLE='Rayleighs/Å',title='UVIS in Black, Model in Red')
   ; p1= plot(x_wave,((3d/1.5d)^2d)*ypts,xtitle='Wavelength (Å)',YTITLE='Rayleighs/Å');,/overplot,color='firebrick');,title='test for iii =  ' + string(iii) + ' and jjj  = ' + string(jjj) )
    
   p1= plot(x_wave(wave_idx_want)/10d,10d*((3d/2.2d)^2d)*ypts(wave_idx_want),xtitle='Wavelength (nm)',YTITLE='Rayleighs/nm');,/overplot,color='firebrick');,title='test for iii =  ' + string(iii) + ' and jjj  = ' + string(jjj) )

   ;p1= plot(x_wave(wave_idx_want),((3d/2d)^2d)*ypts(wave_idx_want),xtitle='Wavelength (A)',YTITLE='Rayleighs/A');,/overplot,color='firebrick');,title='test for iii =  ' + string(iii) + ' and jjj  = ' + string(jjj) )

   
  ; p1= plot(x_wave(wave_idx_want),()*((3d/2.2d)^2d)*ypts(wave_idx_want),xtitle='Wavelength (A)',YTITLE='Rayleighs/A');,/overplot,color='firebrick');,title='test for iii =  ' + string(iii) + ' and jjj  = ' + string(jjj) )

   


   counts1 = dblarr(nspec_bins)
   for mmm=0, nspec_bins - 1 do begin
     counts1(mmm) = 1800d*factor[mmm]*((3d/2.2d)^2d)*ypts[mmm]
   endfor
   
   
   
   
   ;;;;;;;;
   
   
   norm_vec = [-sin(phi) , Cos(phi) , 0d]; [0d,1d,0d];
   slit_pos_vec =  [9.4*Cos(phi) + 10d*Sin(phi) , 9.4*sin(phi) - 10d*cos(phi) , 0d];[x_SC[jjj,iii],y_SC[jjj,iii],z_SC[jjj,iii]];[x0t,y0t,z0];[x0,y0,z0] ; [6d,-10d,0d];
   n = calc_LOS_n_s_and_Tel_s_given_my_model_3D_cart_grid_diffthreee_for_faraway( xgrid, ygrid, zgrid, norm_vec, slit_pos_vec, n_in, T_in, s=s, tel=tel)



   ypts = new_CITEP_2_given_n_s_and_tel_s_using_emission_tables_uvs(s,n,Tel,fwhm,spec_binsize,yptsi_in_2dto1d,nel2dto1d,tec2dto1d,yname_in,x_wave, xwavi_in)
   ; ypts_out(jjj,iii,*) = ypts
ypts2 = ypts
   ; p2 = plot(xxx,yyy,xtitle='Wavelength (Å)',YTITLE='Rayleighs/Å',title='UVIS in Black, Model in Red')
   ; p1= plot(x_wave,((3d/1.5d)^2d)*ypts,xtitle='Wavelength (Å)',YTITLE='Rayleighs/Å');,/overplot,color='firebrick');,title='test for iii =  ' + string(iii) + ' and jjj  = ' + string(jjj) )

   p1= plot(x_wave(wave_idx_want)/10d,10d*((3d/2.2d)^2d)*ypts(wave_idx_want),xtitle='Wavelength (nm)',YTITLE='Rayleighs/nm');,/overplot,color='firebrick');,title='test for iii =  ' + string(iii) + ' and jjj  = ' + string(jjj) )

   ;p1= plot(x_wave(wave_idx_want),((3d/2d)^2d)*ypts(wave_idx_want),xtitle='Wavelength (A)',YTITLE='Rayleighs/A');,/overplot,color='firebrick');,title='test for iii =  ' + string(iii) + ' and jjj  = ' + string(jjj) )


   ; p1= plot(x_wave(wave_idx_want),()*((3d/2.2d)^2d)*ypts(wave_idx_want),xtitle='Wavelength (A)',YTITLE='Rayleighs/A');,/overplot,color='firebrick');,title='test for iii =  ' + string(iii) + ' and jjj  = ' + string(jjj) )




   counts2 = dblarr(nspec_bins)
   for mmm=0, nspec_bins - 1 do begin
     counts2(mmm) = 1800d*factor[mmm]*((3d/2.2d)^2d)*ypts[mmm]
   endfor
   
   p001=plot(x_wave(wave_idx_want)/10d,counts1(wave_idx_want),xtitle='Wavelength (nm)',YTITLE='Counts/Pixel',title='Simulated 30 min Europa-UVS Spectral Image of Io Plasma Torus',/ylog,xrange=[50,210],yrange=[1e-4,1e4],name='At Io ($\rho_c = 5.9 R_J$)' )

 p002=plot(x_wave(wave_idx_want)/10d,counts2(wave_idx_want),/overplot,/ylog,xrange=[50,210],color='firebrick',name='At Europa ($\rho_c = 9.4 R_J$)',yrange=[1e-4,1e4] )
 
 leg = LEGEND(TARGET=[p001,p002], POSITION=[2.1e2,1e4], /DATA, /AUTO_TEXT_COLOR)
 
 
 
   
   leg.save,'Io_vs_europa_Counts_Io_torus_emission_30min_exposure_best_case_scenario.png',resolution=500
    ;stop
 ; endfor
 ;write_csv,'CLIPPER_UVS_Rayleighs_per_ang_num_et_each_interval=1801_num_intervals=45_num_xwavebins=1501.csv',ypts_out
;endfor

;save,slit_pos_vec_outz,norm_vec_outz,s_outz,x_outz,y_outz,z_outz,utc_outz,tec_outz,nec_outz,nsp_outz,ns2p_outz,ns3p_outz,nop_outz,no2p_outz,filename='positions_and_densities_to_plot_for_3_et_intervals_for_each_of_45_Clipper_observations_v1.sav',/verbose

write_csv,'wave_euvs.csv',x_wave(wave_idx_want)/10d

write_csv,'io_counts_euvs.csv',counts1(wave_idx_want)


write_csv,'europa_counts_euvs.csv',counts2(wave_idx_want)



write_csv,'io_raypernm_euvs.csv',10d*((3d/2.2d)^2d)*ypts1(wave_idx_want)


write_csv,'europa_raypernm_euvs.csv',10d*((3d/2.2d)^2d)*ypts2(wave_idx_want)







  stop
  




end




