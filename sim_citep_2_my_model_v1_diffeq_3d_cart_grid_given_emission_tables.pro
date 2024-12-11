
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
function calc_LOS_n_s_and_Tel_s_given_my_model_3D_cart_grid, xgrid, ygrid, zgrid, norm_vec, slit_pos_vec, n_in, T_in, s=s, tel=tel
;, xgrid,ygrid,zgrid,norm_vec,slit_pos_vec,n_in,t_in,L_in,s=s,tel=tel

  ;tic
  IF (norm_vec[0] ne 0d) AND (norm_vec[1] ne 0d) AND (norm_vec[2] ne 0d) THEN BEGIN
    s_x = (xgrid - slit_pos_vec[0] )/norm_vec[0]
    s_y = (ygrid - slit_pos_vec[1] )/norm_vec[1]
    s_z = (zgrid - slit_pos_vec[2] )/norm_vec[2]
    ss=[s_x,s_y,s_z]
    sf=max(ss)
    sm=min(ss)
  ENDIF ELSE Begin


    IF (norm_vec[0] eq 0d) AND (norm_vec[1] ne 0d) and (norm_vec[2] ne 0d) THEN BEGIN
      s_y = (ygrid - slit_pos_vec[1] )/norm_vec[1]
      s_z = (zgrid - slit_pos_vec[2] )/norm_vec[2]
      ss=[s_y,s_z]
      sf=max(ss)
      sm=min(ss)
    ENDIF

    IF (norm_vec[0] ne 0d) AND (norm_vec[1] eq 0d) and (norm_vec[2] ne 0d) THEN BEGIN
      s_x = (xgrid - slit_pos_vec[0] )/norm_vec[0]
      s_z = (zgrid - slit_pos_vec[2] )/norm_vec[2]
      ss=[s_x,s_z]
      sf=max(ss)
      sm=min(ss)
    ENDIF

    IF (norm_vec[0] ne 0d) AND (norm_vec[1] ne 0d) and (norm_vec[2] eq 0d) THEN BEGIN
      s_x = (xgrid - slit_pos_vec[0] )/norm_vec[0]
      s_y = (ygrid - slit_pos_vec[1] )/norm_vec[1]
      ss=[s_x,s_y]
      sf=max(ss)
      sm=min(ss)
    ENDIF


    IF (norm_vec[0] eq 0d) AND (norm_vec[1] eq 0d) and (norm_vec[2] ne 0d) THEN BEGIN
      s_z = (zgrid - slit_pos_vec[2] )/norm_vec[2]
      ss=[s_z]
      sf=max(ss)
      sm=min(ss)
    ENDIF

    IF (norm_vec[0] ne 0d) AND (norm_vec[1] eq 0d) and (norm_vec[2] eq 0d) THEN BEGIN
      s_x = (xgrid - slit_pos_vec[0] )/norm_vec[0]
      ss=[s_x]
      sf=max(ss)
      sm=min(ss)
    ENDIF

    IF (norm_vec[0] eq 0d) AND (norm_vec[1] ne 0d) and (norm_vec[2] eq 0d) THEN BEGIN
      s_y = (ygrid - slit_pos_vec[1] )/norm_vec[1]
      ss=[s_y]
      sf=max(ss)
      sm=min(ss)
    ENDIF

  ENDELSE



  s=ss[UNIQ(ss, SORT(ss))]


  ;Print,sf,sm
  ;print,s

  ;stop



  n_a=n_elements(s)

  ;calculate x,y, & zvals along line of sight (LOS)

  x_los=slit_pos_vec[0] + s*norm_vec[0]; in RJ
  y_los=slit_pos_vec[1] + s*norm_vec[1]; in RJ
  z_los=slit_pos_vec[2] + s*norm_vec[2]; in RJ

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





  ;B_intgrnd_values=dblarr(n_a,num_emiss)

  ;print,n_elements(emiss_R0)
  ;stop
  ;toc;
  ;stop
  ;;;;;;;;;;;;
  ;tic
  

      
      ;ii_= where((L_in ge (rho - eps)) and (L_in le (rho + eps)) )
      
      ;ii_=min(where( abs(l_in - rho) eq  min(abs(l_in - rho))) )
      ;L_shells= 5.+ 0.01*findgen(501)
      
      ;l_shell_los = r_los/cos((!dpi/180d)*lat_los)^2   ; only works for aligned equators as this is for magnetic latitude
       
       tel = fltarr(n_a)
       
       for a = 0, n_a - 1  do begin
       
      ii = where(xgrid eq x_los(a))
    
      
      jj = where(ygrid eq y_los(a))
      
      kk = where(zgrid eq z_los(a))
      
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








Function new_CITEP_2_given_n_s_and_tel_s_using_emission_tables,s,n,Tel_s,x_wav,fwhm,spec_binsize,yptsi_,nel_,tec_,xwavi,yname

wav_min=min(x_wav)
wav_max=max(x_wav)
n_a=n_elements(tel_s)
num_emiss=n_elements(xwavi)
B_intgrnd_values=dblarr(n_a,num_emiss)
zero_idxs=where((Tel_s eq 0d) or (n.el eq 0d) )
neq_zero_idxs=where((Tel_s ne 0d) and (n.el ne 0d))
  nel_s=reform(n.el)
   nsp_s=reform(n.sp)
    ns2p_s=reform(n.s2p)
     ns3p_s=reform(n.s3p)
      nop_s=reform(n.op)
       no2p_s=reform(n.o2p)
       
       
       sp_idxs = where(yname eq '"S  II"')
       s2p_idxs = where(yname eq '"S  III"')
       s3p_idxs = where(yname eq '"S  IV"')
       s4p_idxs = where(yname eq '"S  V"')
       op_idxs = where(yname eq '"O  II"')
       o2p_idxs = where(yname eq '"O  III"')
 
;stop
if n_elements(zero_idxs) ne n_elements(nel_s) then begin
  
  
    remove, zero_idxs, Tel_s, nel_s, nsp_s, ns2p_s, ns3p_s, nop_s, no2p_s
  
  

yptsi_los = dblarr(n_elements(nel_s),num_emiss)
    ;yptsi_,nel_,tec_,xwavi,yname
 for j=0,num_emiss-1 do begin
  yptsi_LOS(*,j) = griddata(  nel_, tec_,reform(yptsi_(*,j)),xout = nel_s, yout = tel_s )
 
  CASE yname(j) OF
    '"S  II"': yptsi_LOS(*,j) = nsp_s*yptsi_LOS(*,j)
    '"S  III"': yptsi_LOS(*,j) = ns2p_s*yptsi_LOS(*,j)
    '"S  IV"': yptsi_LOS(*,j) = ns3p_s*yptsi_LOS(*,j)
    '"S  V"': yptsi_LOS(*,j) = 0.*yptsi_LOS(*,j)
    '"O  II"': yptsi_LOS(*,j) = nop_s*yptsi_LOS(*,j)
    '"O  III"': yptsi_LOS(*,j) = no2p_s*yptsi_LOS(*,j)
  ENDCASE
 
 Endfor

;stop

  if (zero_idxs(0) ge 0) then for i=0,n_elements(zero_idxs)-1 do B_intgrnd_values(zero_idxs(i),*)=replicate(0d,num_emiss)
  if (neq_zero_idxs(0) ge 0) then for i=0,n_elements(neq_zero_idxs)-1 do B_intgrnd_values(neq_zero_idxs(i),*)=reform(yptsi_LOS(i,*))



  yptsi=dblarr(num_emiss)
  for j=0,num_emiss-1 do begin


    yptsi(j)=7.1492d3*INT_TABULATED(s, reform(B_intgrnd_values(*,j)),/double) ; no factor of 2 this time

  endfor

  ypts=simulate_IPT_spectrum_Rayleighs_citep2_ERF_form(x_wav,spec_binsize, xwavi,yptsi, $
    fwhm = fwhm)

  
endif else begin
  ypts=replicate(0d,n_elements(x_wav))
endelse
 

  return, ypts

End


Pro sim_citep_2_my_model_v1_diffeq_3D_cart_grid_given_emission_tables

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


  Restore, 'Mymodelv1_3D_JRM33+con20_cartesian_irrcart_mapped_no_interp_2degree_phi_into_rhoIII_4.5_281x281x129.sav',/verbose
;;

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

z_half_step_sizegrid1 =  [ replicate(dz1/2.,nz1) , replicate(dz2/2.,nz2) , replicate(dz1/2.,nz1 +1) ]
x_half_step_sizegrid1 =  [ replicate(dx1/2.,nx1) , replicate(dx2/2.,nx2) , replicate(dx1/2.,nx3) , replicate(dx2/2.,nx2) , replicate(dx1/2.,nx1+1) ]
ygrid1d = xgrid1d
y_half_step_sizegrid1  = x_half_step_sizegrid1

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
  
  
  n_in=replicate({densities, s: 0.0, sp: 0.0, s2p: 0.0, s3p: 0.0, s4p:0.0, $
    o: 0.0, op: 0.0,oph: 0.0, o2p: 0.0, el: 0.0,elh:0.0},n_elements(xgrid),n_elements(ygrid),n_elements(zgrid)) ;(cm^-3) 501 L shells 5-10RJ 0.01RJ res 1001 lat -50-50 degrees 0.1 degree res

  T_in=replicate({Temperatures, s: 0.0, sp: 0.0, s2p: 0.0, s3p: 0.0, s4p:0.0, $
    o: 0.0, op: 0.0,oph: 0.0, o2p: 0.0, el: 0.0,elh:0.0},n_elements(xgrid),n_elements(ygrid),n_elements(zgrid)) ;(eV), constant over field lines for maxwellians as this case is
    
  

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
  
 ypts_vary_x0_and_z0=dblarr(n_elements(xgrid),n_elements(zgrid), n_elements(x))
    ypts_total_vary_x0_and_z0=dblarr(n_elements(xgrid),n_elements(zgrid))

    nel = [dindgen(9) + 1d, 10d*dindgen(9) + 10d, 100d*dindgen(4) + 100d, 250d*dindgen(39) + 500d]

    Tec = [0.1d*dindgen(9) + 0.1d, 0.5d*dindgen(17)+1d,dindgen(10)+10d,5d*dindgen(8) + 20d, 20d*dindgen(7) + 60d, 80d*dindgen(6) + 200d]

    n_ne = n_elements(nel)

    n_te = n_elements(tec)
    
    nel2dto1d =dblarr(n_ne*n_te)
    tec2dto1d =dblarr(n_ne*n_te)
    yptsi_in_2dto1d = dblarr(n_ne*n_te,884)

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


tic
For i = 0, n_elements(xgrid)-1 do begin ; for x= 4.5 to 10
  For j = 0, n_elements(zgrid)-1 do begin ; for z=0 to z= 3

    x0=xgrid(i)
    z0=zgrid(j)
    ;y0=-10
     print,x0,y0,z0
    slit_pos_vec=[x0,y0,z0]
    
    
 
    
;tic
    
    
    n = calc_LOS_n_s_and_Tel_s_given_my_model_3D_cart_grid( xgrid, ygrid, zgrid, norm_vec, slit_pos_vec, n_in, T_in, s=s, tel=tel)
   ; n = calc_LOS_n_s_and_Tel_s_given_my_model_3D_cart_grid(xgrid,ygrid,zgrid,norm_vec,slit_pos_vec,n_in,T_in,s=s,tel=tel)
   ; toc
    ;stop
  ;  tic
    ypts=new_CITEP_2_given_n_s_and_tel_s_using_emission_tables(s,n,Tel,x,fwhm,binsize,yptsi_in_2dto1d,nel2dto1d,tec2dto1d,xwavi_in,yname_in)
    ;print,imin,jmin
    
    ;p0=plot(x_wav_uvis,y_uvis,xtitle='Wavelength (Å)',ytitle='Rayleighs/Å',title='UVIS 6 RJ in Black, CITEP2 Model1 in Red',xrange=[550d,2100d])
   ;p1=plot(x,ypts,xtitle='Wavelength (Å)',ytitle='Rayleighs/Å',title='x0=6,y0=-10,z0=0, pointing in yhat direction');,xrange=[1200,1800]);,xrange=[550d,2100d]);,/overplot,color='red')
   ; toc
   ; worked!@!!!!!!!!!!! in  seconds at 0.1RJ no hote in 24.81 seconds from terminal
    
   ; stop
    ypts_vary_x0_and_z0(i,j ,*)=ypts
   ypts_total_vary_x0_and_z0(i,j)=total(ypts);,/double)
   ; stop


  Endfor
Endfor

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
  
  write_csv,'ypts_vary_x0_and_z0_emission_table_281x129x884.txt',ypts_vary_x0_and_z0
    write_csv,'ypts_total_vary_x0_and_z0_emission_table_281x129.txt',ypts_total_vary_x0_and_z0
  toc
  stop
  




end




