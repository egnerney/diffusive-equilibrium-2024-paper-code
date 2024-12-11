function calculate_IPT_emiss_te_same_size_ne, Tel, nel, Nsp, Ns2p, Ns3p,Ns4p, Nop, No2p, $
  min = min, max = max, xwavi=xwavi
  ;Tel in eV, nel in #/cm^3, Nsp etc in #/cm^2, max/min/fwhm in angstroms
  ; sample spectra
  if keyword_set(max) then maxwav = max else maxwav = 1800d
  if keyword_set(min) then minwav = min else minwav = 550d

  num_tel=n_elements(tel)
  conversion=11604.5221d ; conversion between eV and Kelvin from NIST
  p_ratio=replicate(0.1d,num_tel)
  log10tel=alog10(Tel*conversion)
  log10nel=alog10(nel)

  ;tic

  ; emission line calculation: Sulfur ions
  s2em = emiss_calc_egn(16, 2, temp = log10tel, dens = log10nel, /no_de, radt = 1.d, /quiet,proton_ratio=p_ratio)
  s3em = emiss_calc_egn(16, 3, temp = log10tel, dens = log10nel, /no_de, radt = 1.d, /quiet,proton_ratio=p_ratio)
  s4em = emiss_calc_egn(16, 4, temp = log10tel, dens = log10nel, /no_de, radt = 1.d, /quiet,proton_ratio=p_ratio)
  s5em = emiss_calc_egn(16, 5, temp = log10tel, dens = log10nel, /no_de, radt = 1.d, /quiet,proton_ratio=p_ratio)


  ; emission line calculation: oxygen ions
  o2em = emiss_calc_egn(8, 2, temp = log10tel, dens = log10nel, /no_de, radt = 1.d, /quiet,proton_ratio=p_ratio)
  o3em = emiss_calc_egn(8, 3, temp = log10tel, dens = log10nel, /no_de, radt = 1.d, /quiet,proton_ratio=p_ratio)
  ;toc



  ;tic


  xwavi = [s2em.lambda, s3em.lambda, s4em.lambda,s5em.lambda, o2em.lambda, o3em.lambda]
  wsort = sort(xwavi) ; sort by wavelength
  xwavi = xwavi[wsort]
  yname = [s2em.ion_name, s3em.ion_name, s4em.ion_name, s5em.ion_name, o2em.ion_name, o3em.ion_name]
  yname = yname[wsort]

  ; cut down arrays to necessary range
  avgwav = (minwav + maxwav)/2.d
  wrange = where(abs(xwavi - avgwav) le avgwav - minwav)
  xwavi = xwavi[wrange]

  yname = yname[wrange]

  ;tic



  yptsi=dblarr(num_tel,n_elements(xwavi))
  for i=0, num_tel - 1 do begin



    yptsi_temp = [reform(s2em.em[i,i,*])*Nsp[i],reform(s3em.em[i,i,*])*Ns2p[i], reform(s4em.em[i,i,*])*Ns3p[i], reform(s5em.em[i,i,*])*Ns4p[i], reform(o2em.em[i,i,*])*Nop[i], reform(o3em.em[i,i,*])*No2p[i]]
    yptsi_temp= yptsi_temp[wsort]
    yptsi[i,*] =yptsi_temp[wrange]

  endfor

  ;toc
  ;stop

  return, yptsi

  ; note: xwavi is also an output, placed as a keyword

end

function calculate_IPT_emiss_citep2_output_yname, Tel, nel, Nsp, Ns2p, Ns3p,Ns4p, Nop, No2p, $
  min = min, max = max, xwavi=xwavi, yname = yname
  ;Tel in eV, nel in #/cm^3, Nsp etc in #/cm^2, max/min/fwhm in angstroms
  ; sample spectra
  if keyword_set(max) then maxwav = max else maxwav = 1800d
  if keyword_set(min) then minwav = min else minwav = 550d

  ;p_ratio=replicate(0.1d,n_elements(tel))
  conversion=11604.5221d ; conversion between eV and Kelvin from NIST
  log10tel=alog10(Tel*conversion)
  log10nel=alog10(nel)

  ;tic
  ; emission line calculation: Sulfur ions
  s2em = emiss_calc_egn(16, 2, temp = log10tel, dens = log10nel, /no_de, radt = 1.d, /quiet,proton_ratio=0.1d)
  s3em = emiss_calc_egn(16, 3, temp = log10tel, dens = log10nel, /no_de, radt = 1.d, /quiet,proton_ratio=0.1d)
  s4em = emiss_calc_egn(16, 4, temp = log10tel, dens = log10nel, /no_de, radt = 1.d, /quiet,proton_ratio=0.1d)
  s5em = emiss_calc_egn(16, 5, temp = log10tel, dens = log10nel, /no_de, radt = 1.d, /quiet,proton_ratio=0.1d)


  ; emission line calculation: oxygen ions
  o2em = emiss_calc_egn(8, 2, temp = log10tel, dens = log10nel, /no_de, radt = 1.d, /quiet,proton_ratio=0.1d)
  o3em = emiss_calc_egn(8, 3, temp = log10tel, dens = log10nel, /no_de, radt = 1.d, /quiet,proton_ratio=0.1d)
  ;toc
  ;help,s3em.em
  ;stop


  ;puts in rayleighs, nsp etc. in local density units #/cm^3
  s2emiss = reform(s2em.em)*Nsp
  s3emiss = reform(s3em.em)*Ns2p
  s4emiss = reform(s4em.em)*Ns3p
  s5emiss=  reform(s5em.em)*Ns4p
  o2emiss = reform(o2em.em)*Nop
  o3emiss = reform(o3em.em)*No2p

  ; discrete
  xwavi = [s2em.lambda, s3em.lambda, s4em.lambda,s5em.lambda, o2em.lambda, o3em.lambda]
  wsort = sort(xwavi) ; sort by wavelength
  xwavi = xwavi[wsort]

  yptsi = [s2emiss, s3emiss, s4emiss, s5emiss, o2emiss, o3emiss]
  yptsi = yptsi[wsort]

  yname = [s2em.ion_name, s3em.ion_name, s4em.ion_name, s5em.ion_name, o2em.ion_name, o3em.ion_name]
  yname = yname[wsort]

  ; cut down arrays to necessary range
  avgwav = (minwav + maxwav)/2.d
  wrange = where(abs(xwavi - avgwav) le avgwav - minwav)
  xwavi = xwavi[wrange]
  yptsi = yptsi[wrange]
  yname = yname[wrange]

  return, yptsi

  ; note: xwavi is also an output, placed as a keyword

end

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




function calc_LOS_n_s_and_Tel_s_azimuthally_sym_scale_heights, xgrid,ygrid,zgrid,norm_vec,slit_pos_vec,an,bn,atel,btel,s=s,tel=tel

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

  x_los=slit_pos_vec[0] + s*norm_vec[0]
  y_los=slit_pos_vec[1] + s*norm_vec[1]
  z_los=slit_pos_vec[2] + s*norm_vec[2]



  

  x_index=make_array(n_a,/integer)
  y_index=make_array(n_a,/integer)
  z_index=make_array(n_a,/integer)
  
  T=replicate({Temps_LOS, s: 0.0d, sp: 0.0d, s2p: 0.0d, s3p: 0.0d, s4p:0.0d, $
    o: 0.0d, op: 0.0d, o2p: 0.0d, el: 0.0d,elh:0.0d},n_a)
    
    n=replicate({densities_LOS, s: 0.0d, sp: 0.0d, s2p: 0.0d, s3p: 0.0d, s4p:0.0d, $
    o: 0.0d, op: 0.0d, o2p: 0.0d, el: 0.0d,elh:0.0d},n_a)


  


  ;B_intgrnd_values=dblarr(n_a,num_emiss)

  ;print,n_elements(emiss_R0)
  ;stop
  ;toc;
  ;stop
  ;;;;;;;;;;;;
  ;tic
  for a=0,n_a-1 do begin
    
   rho=sqrt(x_los(a)^2d + y_los(a)^2d)
    
    If (rho ge 6d) and (rho le 10d) then begin 
      T[a].sp=-2992.71d + 1401.02d*rho - 269.568d*(rho^2d) + 27.3165d*(rho^3d) - 1.12374d*(rho^4d)
      T[a].s2p=-4650.83d + 2122.14d*rho  - 351.473d*(rho^2d) + 25.9979d*(rho^3d)  - 0.732711d*(rho^4d)
      T[a].s3p=-2826.d + 1232.17d*rho  - 188.612d*(rho^2d) + 12.4472d*(rho^3d)  - 0.299455d*(rho^4d)
      T[a].op=-4042.21d + 1885.72d*rho  - 321.068d*(rho^2d) + 24.6705d*(rho^3d)  - 0.728725d*(rho^4d)
      T[a].o2p=-1766.47d + 695.67d*rho  - 88.4991d*(rho^2d) + 4.25863d*(rho^3d)  - 0.0512119d*(rho^4d)

   
      T[a].el=atel*((sqrt(x_los(a)^2d + y_los(a)^2d)/6d)^btel)
      H=get_scale_heights_citep2_inRJ_single(T[a])
      
      n[a].sp=an[0]*((sqrt(x_los(a)^2d + y_los(a)^2d)/6d)^bn[0])*(Exp(-(z_los(a)/H.sp)^2d))
      n[a].s2p=an[1]*((sqrt(x_los(a)^2d + y_los(a)^2d)/6d)^bn[1])*(Exp(-(z_los(a)/H.s2p)^2d))
      n[a].s3p=an[2]*((sqrt(x_los(a)^2d + y_los(a)^2d)/6d)^bn[2])*(Exp(-(z_los(a)/H.s3p)^2d))
      n[a].op=an[3]*((sqrt(x_los(a)^2d + y_los(a)^2d)/6d)^bn[3])*(Exp(-(z_los(a)/H.op)^2d))
      n[a].o2p=an[4]*((sqrt(x_los(a)^2d + y_los(a)^2d)/6d)^bn[4])*(Exp(-(z_los(a)/H.o2p)^2d))

   
    endif 
    
 

 endfor
 
    
    
    n.el=(n.sp + 2d*n.s2p + 3d*n.s3p + n.op + 2d*n.o2p)/0.9d
    tel=reform(T.el)
    return,n

end

function calc_LOS_n_s_and_Tel_s_azimuthally_sym_scale_heights_exact_same_before, xgrid,ygrid,zgrid,norm_vec,slit_pos_vec,an,bn,atel,btel,s=s,tel=tel

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

  x_los=slit_pos_vec[0] + s*norm_vec[0]
  y_los=slit_pos_vec[1] + s*norm_vec[1]
  z_los=slit_pos_vec[2] + s*norm_vec[2]





  x_index=make_array(n_a,/integer)
  y_index=make_array(n_a,/integer)
  z_index=make_array(n_a,/integer)

  T=replicate({Temps_LOS, s: 0.0d, sp: 0.0d, s2p: 0.0d, s3p: 0.0d, s4p:0.0d, $
    o: 0.0d, op: 0.0d, o2p: 0.0d, el: 0.0d,elh:0.0d},n_a)

  n=replicate({densities_LOS, s: 0.0d, sp: 0.0d, s2p: 0.0d, s3p: 0.0d, s4p:0.0d, $
    o: 0.0d, op: 0.0d, o2p: 0.0d, el: 0.0d,elh:0.0d},n_a)





  ;B_intgrnd_values=dblarr(n_a,num_emiss)

  ;print,n_elements(emiss_R0)
  ;stop
  ;toc;
  ;stop
  ;;;;;;;;;;;;
  ;tic
  
  anel=2200d
  anel2=anel*((7.8d/6d)^(-5.4d))
  for a=0,n_a-1 do begin

    rho=sqrt(x_los(a)^2d + y_los(a)^2d)

    If (rho ge 6d) and (rho le 10d) then begin
      T[a].sp=-2992.71d + 1401.02d*rho - 269.568d*(rho^2d) + 27.3165d*(rho^3d) - 1.12374d*(rho^4d) ;eV
      T[a].s2p=-4650.83d + 2122.14d*rho  - 351.473d*(rho^2d) + 25.9979d*(rho^3d)  - 0.732711d*(rho^4d)
      T[a].s3p=-2826.d + 1232.17d*rho  - 188.612d*(rho^2d) + 12.4472d*(rho^3d)  - 0.299455d*(rho^4d)
      T[a].op=-4042.21d + 1885.72d*rho  - 321.068d*(rho^2d) + 24.6705d*(rho^3d)  - 0.728725d*(rho^4d)
      T[a].o2p=-1766.47d + 695.67d*rho  - 88.4991d*(rho^2d) + 4.25863d*(rho^3d)  - 0.0512119d*(rho^4d)


      T[a].el=atel*((sqrt(x_los(a)^2d + y_los(a)^2d)/6d)^btel)
      H=get_scale_heights_citep2_inRJ_single(T[a])

       If (rho le 7.8d) then n[a].el = anel*((rho/6d)^(-5.4d))
       
       If (rho gt 7.8d) then n[a].el = anel2*((rho/7.8d)^(-12d))
       
      n[a].sp=n[a].el*an[0]*((rho/6d)^bn[0])*(Exp(-(z_los(a)/H.sp)^2d))
      n[a].s2p=n[a].el*an[1]*((rho/6d)^bn[1])*(Exp(-(z_los(a)/H.s2p)^2d))
      n[a].s3p=n[a].el*an[2]*((rho/6d)^bn[2])*(Exp(-(z_los(a)/H.s3p)^2d))
      n[a].op=n[a].el*an[3]*((rho/6d)^bn[3])*(Exp(-(z_los(a)/H.op)^2d))
      n[a].o2p=n[a].el*an[4]*((rho/6d)^bn[4])*(Exp(-(z_los(a)/H.o2p)^2d))


    endif



  endfor



  n.el=(n.sp + 2d*n.s2p + 3d*n.s3p + n.op + 2d*n.o2p)/0.9d
  tel=reform(T.el)
  return,n

end




Function new_CITEP_2_given_n_s_and_tel_s,s,n,Tel_s,x_wav,fwhm,spec_binsize

wav_min=min(x_wav)
wav_max=max(x_wav)
  n_a=n_elements(tel_s)

  feh0=0.0025d
  tel0= 5.0d
  nel0=2000.d
  nsp0=0.06d*2000.d
  ns2p0=0.21d*2000.d
  ns3p0=0.03d*2000.d

  nop0=0.26d*2000.d
  no2p0=0.03d*2000.d


  ;fwhm=4.
  emiss_R0=calculate_IPT_emiss_citep2( Tel0, nel0, Nsp0, Ns2p0, Ns3p0,0d, Nop0, No2p0, $
    min = wav_min, max = wav_max, xwavi=xwavi) ; just to get proper length of array and xwavi values

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
 
;stop
if n_elements(zero_idxs) ne n_elements(nel_s) then begin
  
  
    remove, zero_idxs, Tel_s, nel_s, nsp_s, ns2p_s, ns3p_s, nop_s, no2p_s
  
  



  yptsi_LOS=calculate_IPT_emiss_te_same_size_ne(Tel_s, nel_s, nsp_s, ns2p_s, ns3p_s,0d*nel_s, nop_s, no2p_s, $
    min = wav_min, max = wav_max,xwavi=xwavi)


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

FUNCTION MYFITTER_cassini_euv_and_fuv_citep_2,p,x=x,y=y,err=err
  ; Parameter values are passed in "p"
  
  x_min=-10d;
  x_max=10d;
  x_step_size=0.2d ;  RJ
  n_x=Round((x_max - x_min)/x_step_size) + 1
  xgrid=x_step_size*Dindgen(n_x) + x_min

  y_min=-10d;
  y_max=10d;
  y_step_size=x_step_size ;  RJ
  n_y=Round((y_max - y_min)/y_step_size) + 1
  ygrid=y_step_size*Dindgen(n_y) + y_min

  z_min=-4d;
  z_max=4d;
  z_step_size=x_step_size ;  RJ
  n_z=Round((z_max - z_min)/z_step_size) + 1
  zgrid=z_step_size*Dindgen(n_z) + z_min
  
  fwhm=4.47d
  ;x_wav=spec_binsize*dindgen(n_wav) + wav_min
  ;;;

  ;;;
  norm_vec=[0d,1d,0d]
  x0=6d
  y0=-10d
  z0=0d
  slit_pos_vec=[x0,y0,z0]
  

  euv_bin_size=0.6049d ; for euv Cassini UVIS (550-1150A)
  fuv_bin_size=0.7794d;

  n=calc_LOS_n_s_and_Tel_s_azimuthally_sym_scale_heights( xgrid,ygrid,zgrid,norm_vec,slit_pos_vec,[p[0],p[1],p[2],p[3],p[4]],[p[5],p[6],p[7],p[8],p[9]],p[10],p[11],s=s,tel=tel)
  ;toc; xgrid,ygrid,zgrid,norm_vec,slit_pos_vec,an,bn,atel,btel,s=s,tel=tel

  ;tic
  model=new_CITEP_2_given_n_s_and_tel_s(s,n,Tel,x,fwhm,euv_bin_size);s,n,Tel_s,x_wav,fwhm,spec_binsize
  ;toc
  ;help,ypts,euv_wav,x_wave
  ;tic
  idx_fuv_min=where((x ge 1150d - euv_bin_size/2d) and (x le 1150d + euv_bin_size/2d))

  ;print,x_wav_uvis(idx_fuv_min)
  ;stop
  model[idx_fuv_min + 1:n_elements(x)-1]=model[idx_fuv_min + 1:n_elements(x)-1]*(euv_bin_size/fuv_bin_size)

  ;toc


  ;model=ypts

  return, (y-model)/err
END

FUNCTION MYFITTER_cassini_fuv_citep_2,p,x=x,y=y,err=err
  ; Parameter values are passed in "p"

  x_min=-10d;
  x_max=10d;
  x_step_size=0.2d ;  RJ
  n_x=Round((x_max - x_min)/x_step_size) + 1
  xgrid=x_step_size*Dindgen(n_x) + x_min

  y_min=-10d;
  y_max=10d;
  y_step_size=x_step_size ;  RJ
  n_y=Round((y_max - y_min)/y_step_size) + 1
  ygrid=y_step_size*Dindgen(n_y) + y_min

  z_min=-4d;
  z_max=4d;
  z_step_size=x_step_size ;  RJ
  n_z=Round((z_max - z_min)/z_step_size) + 1
  zgrid=z_step_size*Dindgen(n_z) + z_min

 

  ;;;
  norm_vec=[0d,1d,0d]
  x0=6d
  y0=-10d
  z0=0d
  slit_pos_vec=[x0,y0,z0]


  euv_bin_size=0.6049d ; for euv Cassini UVIS (550-1150A)
  fuv_bin_size=0.7794d;

  n=calc_LOS_n_s_and_Tel_s_azimuthally_sym_scale_heights( xgrid,ygrid,zgrid,norm_vec,slit_pos_vec,[p[0],p[1],p[2],p[3],p[4]],[p[5],p[6],p[7],p[8],p[9]],p[10],p[11],s=s,tel=tel)
  ;toc; xgrid,ygrid,zgrid,norm_vec,slit_pos_vec,an,bn,atel,btel,s=s,tel=tel

  ;tic
  model=new_CITEP_2_given_n_s_and_tel_s(s,n,Tel,x,fwhm,fuv_bin_size);s,n,Tel_s,x_wav,fwhm,spec_binsize
  ;toc
  ;help,ypts,euv_wav,x_wave
  ;tic
  ;idx_fuv_min=where((x ge 1150d - euv_bin_size/2d) and (x le 1150d + euv_bin_size/2d))

  ;print,x_wav_uvis(idx_fuv_min)
  ;stop
  ;model[idx_fuv_min + 1:n_elements(x)-1]=model[idx_fuv_min + 1:n_elements(x)-1]*(euv_bin_size/fuv_bin_size)

  ;toc


  ;model=ypts

  return, (y-model)/err
END


Pro sim_citep_2_CHIANTI8_tables

tic
 

 ;stop
xmin=550d;3000d;
xmax=2100d;7000d;


nel = [dindgen(9) + 1d, 10d*dindgen(9) + 10d, 100d*dindgen(4) + 100d, 250d*dindgen(39) + 500d]

Tec = [0.1d*dindgen(9) + 0.1d, 0.5d*dindgen(17)+1d,dindgen(10)+10d,5d*dindgen(8) + 20d, 20d*dindgen(7) + 60d, 80d*dindgen(6) + 200d]

n_ne = n_elements(nel)

n_te = n_elements(tec)

nxwavi = 884 ; for 550 - 2100 CHIANTI 8.0.1

yptsi_vary_nel_tec = dblarr(n_ne, n_te, nxwavi)
for i=0,n_ne -1 do begin
  print,'i = ',i
for j=0, n_te -1 do begin  
  
  neli = nel(i)
  tecj = tec(j)
  
    yptsi=calculate_IPT_emiss_citep2_output_yname( Tecj, neli, 1d, 1d, 1d,1d, 1d, 1d, $
    min = xmin, max = xmax, xwavi=xwavi, yname = yname)
    
    
    
   yptsi_vary_nel_tec(i,j ,*)=yptsi
  ; print,xwavi
  ; print,yptsi
  ; print,n_elements(xwavi)
 stop

Endfor
endfor

write_csv,'yptsi_550-2100_Angstroms_vary_ne_te_lookuptable_61x57x884_CHIANTI801.dat',yptsi_vary_nel_tec
write_csv,'xwavi_550-2100_Angstroms_CHIANTI801.dat',xwavi
write_csv,'yname_550-2100_Angstroms_CHIANTI801.dat',yname

sp_idxs = where(yname eq 'S  II') 
s2p_idxs = where(yname eq 'S  III') 
s3p_idxs = where(yname eq 'S  IV') 
s4p_idxs = where(yname eq 'S  V') 
op_idxs = where(yname eq 'O  II') 
o2p_idxs = where(yname eq 'O  III') 
 
  
  
  
  toc
  stop


end




