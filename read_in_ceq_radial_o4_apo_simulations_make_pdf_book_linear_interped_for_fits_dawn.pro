
function given_number_LOS_dens_temp_for_each_layer_simulate_rayleighs_dawn_or_dusk_all_at_once_8APO_lines_fixed_proper_nintegration_now,Tec_a,nec_a,nsp_a,ns2p_a,nop_a,rho_max,rho_min, drho, yptsi_,nel_,tec_

  num_emiss = 8
  num_elements = n_elements(nec_a)
  epsilon_a = dblarr(num_elements,num_emiss)
  ;yptsi_,nel_,tec_,xwavi,yname

  ; stop
  for j=0,num_emiss -1 do begin
    triangulate,nel_,tec_,tr
    epsilon_a(*,j) = griddata(  nel_, tec_,reform(yptsi_(*,j)),xout = nec_a, yout = tec_a,/linear,triangles = tr )

  Endfor


  epsilon_a(*,0) = 1d-6*ns2p_a*epsilon_a(*,0)
  epsilon_a(*,1) = 1d-6*nop_a*epsilon_a(*,1)
  epsilon_a(*,2) = 1d-6*nop_a*epsilon_a(*,2)
  epsilon_a(*,3) = 1d-6*nsp_a*epsilon_a(*,3)
  epsilon_a(*,4) = 1d-6*nsp_a*epsilon_a(*,4)
  epsilon_a(*,5) = 1d-6*ns2p_a*epsilon_a(*,5)
  epsilon_a(*,6) = 1d-6*nsp_a*epsilon_a(*,6)
  epsilon_a(*,7) = 1d-6*nsp_a*epsilon_a(*,7)

  ; rho_c min of line of sight must be decreasing for this algorithm for increasing index fyi
  ;epsilon_a = dblarr(num_elements,8)

  rayleighs = dblarr(num_elements,8)

  ;sigma_n = dblarr(num_elements)


  conversion = 7.1492d9 ; 1RJ  in cm
  ; 4.6 - 7 finite tor
  n_elem_s = floor(2d*rho_max/drho) + 1
  s_LOS = drho*dindgen(n_elem_s)

  n_elem_x0_LOS = floor((rho_max - rho_min)/drho) + 1
  y0 = - rho_max
  x_LOS = drho*dindgen(n_elem_x0_LOS) + rho_min ; x = x_0
  y_LOS = y0 + s_los

  rho_LOS = dblarr(n_elem_x0_LOS,n_elem_s);sqrt( y_LOS^2d + x_LOS^2d )
  idx_wantz = dblarr(n_elem_x0_LOS,n_elem_s)

  x_los = reverse(x_los)

  for i = 0 , n_elem_x0_LOS - 1 do rho_los(i,*) = sqrt( x_LOS(i)^2d + y_LOS^2d )

  for idx = 0,   num_elements - 1 do begin
    idx_want = where((reform(rho_LOS(idx,*)) le rho_max ) and (reform(rho_LOS(idx,*)) ge rho_min ))
    if ((idx_want(0) ge 0) and (n_elements( idx_want) ge 2 )) then begin


      ;where(x_los eq rho_LOS(idx,*))
      for linenum = 0, 7 do begin
        emiss_s = interpol(reform(epsilon_a(*,linenum)),x_LOS,reform(rho_LOS(idx,idx_want))) ; reform(epsilon_a(idx_want,linenum))
        rayleighs(idx,linenum) = conversion*INT_TABULATED( s_LOS(idx_want),emiss_s, /DOUBLE)
      endfor
    endif else begin
      rayleighs(idx,*) = replicate(0d,8)
    endelse


  endfor

  return, rayleighs
  ;sigma_n is given as optional output via keyword argument above

end



pro read_in_ceq_radial_o4_APO_simulations_make_pdf_book_linear_interped_for_fits_dawn
  restore, filename ='grids.sav', /verbose


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  nvalues=5

  xmin_want = 4.6d;5.
  xmax_want = 7.0d

  dx_want = 0.04d


  ;ninterpsteps = round((xinterpmax - xinterpmin)/xinterpstep) + 1
  factor = round(dx_want/0.02d) ; dx_want/og_dx_had
  n_want = round(n_elements(x_grid)/factor)
  idx_want = factor*indgen(n_want)
  ;x_want = ;xinterpstep*dindgen(ninterpsteps) + xinterpmin

  ;ninterpsteps = n_elements(x_want)


  ; interped_dawn_grid=dblarr(ninterpsteps,8)
  ; interped_dusk_grid=dblarr(ninterpsteps,8)
  ; interped_err_dawn_grid=dblarr(ninterpsteps,8)
  ; interped_err_dusk_grid=dblarr(ninterpsteps,8)


  dawn_grid_new = dblarr(n_elements(x_grid),8)
  dusk_grid_new = dblarr(n_elements(x_grid),8)

  err_dawn_grid_new = dblarr(n_elements(x_grid),8)
  err_dusk_grid_new = dblarr(n_elements(x_grid),8)

  for i=0,7 do begin
    for j = 2,n_elements(x_grid)-3 do begin

      ;5 points running average and associated errors/uncertainties assuming Gaussian

      dawn_grid_new[j,i] = (dawn_grid[j - 2,i] + dawn_grid[j - 1,i] + dawn_grid[j,i] + dawn_grid[j + 1,i] + dawn_grid[j + 2,i])/double(nvalues) ; ts_smooth(dawn_grid[*,i],nvalues)

      dusk_grid_new[j,i] =  (dusk_grid[j - 2,i] + dusk_grid[j - 1,i] + dusk_grid[j,i] + dusk_grid[j + 1,i] + dusk_grid[j + 2,i])/double(nvalues) ;ts_smooth(dusk_grid[*,i],nvalues)


      if (dawn_grid_new[j,i] lt 0d) then dawn_grid_new[j,i] = 0.1d

      if (dusk_grid_new[j,i] lt 0d) then dusk_grid_new[j,i] = 0.1d

      ;interped_dawn_grid[j,i] = ;interpol(reform(dawn_grid[*,i]),x_grid, x_interps)

      ;interped_dusk_grid[j,i] = ;interpol(dusk_grid[*,i],x_grid, x_interps)


      err_dawn_grid_new[j,i] = Sqrt( err_dawn_grid[j - 2,i]^2d + err_dawn_grid[j - 1,i]^2d + err_dawn_grid[j,i]^2d + err_dawn_grid[j + 1,i]^2d + err_dawn_grid[j + 2,i]^2d )/double(nvalues);ts_smooth(err_dawn_grid[*,i],nvalues)

      err_dusk_grid_new[j,i] =  Sqrt( err_dusk_grid[j - 2,i]^2d + err_dusk_grid[j - 1,i]^2d + err_dusk_grid[j,i]^2d + err_dusk_grid[j + 1,i]^2d + err_dusk_grid[j + 2,i]^2d )/double(nvalues);ts_smooth(err_dusk_grid[*,i],nvalues)



      ; if (dawn_grid_new[j,i] le 0d) then err_dawn_grid_new[j,i] = 10000d ; big number

      ;if (dusk_grid_new[j,i] le 0d) then err_dusk_grid_new[j,i] = 10000d


      ; interped_err_dawn_grid[j,i] = interpol(err_dawn_grid[*,i],x_grid, x_interps)

      ; interped_err_dusk_grid[j,i] = interpol(err_dusk_grid[*,i],x_grid, x_interps)
    endfor
  endfor

  dawn_grid = dawn_grid_new[idx_want,*] ;interped_dawn_grid
  dusk_grid = dusk_grid_new[idx_want,*];interped_dusk_grid

  err_dawn_grid = err_dawn_grid_new[idx_want,*];interped_err_dawn_grid/sqrt(double(nvalues))  ; assuming same then moving 4 avg then 1/sqrt(4)=1/2 smaller on avg
  err_dusk_grid = err_dusk_grid_new[idx_want,*];interped_err_dusk_grid/sqrt(double(nvalues)) ; assuming same then moving 4 avg then 1/sqrt(4)=1/2 smaller on avg



  x_grid = x_grid[idx_want];x_interps


  idx_want = where((x_grid ge xmin_want) and (x_grid le xmax_want))

  dawn_grid = dawn_grid[idx_want,*] ;interped_dawn_grid
  dusk_grid = dusk_grid[idx_want,*];interped_dusk_grid

  err_dawn_grid = err_dawn_grid[idx_want,*];interped_err_dawn_grid/sqrt(double(nvalues))  ; assuming same then moving 4 avg then 1/sqrt(4)=1/2 smaller on avg
  err_dusk_grid = err_dusk_grid[idx_want,*];interped_err_dusk_grid/sqrt(double(nvalues)) ; assuming same then moving 4 avg then 1/sqrt(4)=1/2 smaller on avg



  x_grid = x_grid[idx_want];x_interps




  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



  ;legday.save,'5point_Moving_avg_0.04res_dawn_and_dusk_ratios_proper_uncs_4.6_to_7.0.png',resolution=300



  ;stop
  ;;;;;;;;;;;;;;;;;;;;;

  ; indx=where((x_grid ge 4.6d) and (x_grid le 7.0d) );indx=where((x_grid ge 5.) and (x_grid le 6.5) )
  ; x_grid=x_grid(indx) ; values between 5-7 RJ

  ;dawn_grid = dawn_grid(indx,*)

  ;dusk_grid = dusk_grid(indx,*)

  ;err_dawn_grid = err_dawn_grid(indx,*)

  ;err_dusk_grid = err_dusk_grid(indx,*)



  ;err_dawn_grid(*,3) =  5d*err_dawn_grid(*,3)
  ;err_dusk_grid(*,3) =  5d*err_dusk_grid(*,3)
  ;err_dawn_grid(*,4) =  20d*err_dawn_grid(*,4)
  ;err_dusk_grid(*,4) =  20d*err_dusk_grid(*,4)






  ;stop






nframes= 90 
p_out = OBJARR(nframes)

;openr,1,'linear_interpedhigher_res_simulated_B6731A_in_each_bin_o4_first_try_90framesx376_ceq_radialbins.csv'
;B=fltarr(nframes,376)
;readf,1,B
;close,1


restore,filename='linear_interpedhigher_res_simulated_apo_8lines_DAWN_90_frames_given_1.2xionsfit_90framesx376_ceq_radialbins_8lines.sav',/verbose
;B_SIM           DOUBLE    = Array[90, 376, 8]
;stop


nframes = double(nframes)


dPhi = 360d/nframes

CMLdeg = dphi*findgen(nframes) + 270d

CMLdeg =  CMLdeg Mod 360d

xRight = 7.5d
xLeft = -7.5d
spacing = 0.04d

n_ceq_radial_bins = Floor((xRight - xLeft)/spacing) + 1;

rho_ceq = spacing*findgen(n_ceq_radial_bins) +  xLeft

p11=errorplot(x_grid,dawn_grid(*,2),ERR_dusk_grid(*,2),NAME='OII (O+) 3729 Angstroms',ERRORBAR_COLOR="blue",color="blue",title='1.2xIon Densities from Dusk Fit',xtitle='Dawn $\rho_c$ ($R_J$)',ytitle='Rayleighs',SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01,yrange=[0,400],xrange=[4.5,7])

p22=errorplot(x_grid,dawn_grid(*,6),ERR_dusk_grid(*,6),NAME='SII (S+) 6716 Angstroms',ERRORBAR_COLOR="red",color="red",/overplot,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)
p33=errorplot(x_grid,dawn_grid(*,7),ERR_dusk_grid(*,7),NAME='SII (S+) 6731 Angstroms',ERRORBAR_COLOR="brown",color="brown",/overplot,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)

p44=errorplot(x_grid,dawn_grid(*,0),ERR_dusk_grid(*,0),NAME='SIII (S++) 3722 Angstroms',ERRORBAR_COLOR="pink",color="pink",/overplot,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)
p55=errorplot(x_grid,dawn_grid(*,1),ERR_dusk_grid(*,1),NAME='OII (O+) 3726 Angstroms',ERRORBAR_COLOR="purple",color="purple",/overplot,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)

p66=errorplot(x_grid,dawn_grid(*,5),ERR_dusk_grid(*,5),NAME='SIII (S++) 6312 Angstroms',ERRORBAR_COLOR="orange",color="orange",/overplot,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)

p77=errorplot(x_grid,dawn_grid(*,3),ERR_dusk_grid(*,3),NAME='SII (S+) 4069 Angstroms',ERRORBAR_COLOR="light Green",color="light Green",/overplot,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)
p88=errorplot(x_grid,dawn_grid(*,4),ERR_dusk_grid(*,4),NAME='SII (S+) 4076 Angstroms',ERRORBAR_COLOR="dark green",color="dark green",/overplot,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)




p1=plot(rho_ceq,mean(reform(B_sim(*,*,2)),dimension = 1),NAME='$O^+$ 3729 Å',color="blue",/overplot)

p2=plot(rho_ceq,mean(reform(B_sim(*,*,6)),dimension = 1),NAME='$S^+$ 6716 Å',color="red",/overplot);,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)
p3=plot(rho_ceq,mean(reform(B_sim(*,*,7)),dimension = 1),NAME='$S^+$ 6731 Å',color="brown",/overplot);,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)

p4=plot(rho_ceq,mean(reform(B_sim(*,*,0)),dimension = 1),NAME='$S^{++}$ 3722 Å',color="pink",/overplot);,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)
p5=plot(rho_ceq,mean(reform(B_sim(*,*,1)),dimension = 1),NAME='$O^+$ 3726 Å',color="purple",/overplot);,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)

p6=plot(rho_ceq,mean(reform(B_sim(*,*,5)),dimension = 1),NAME='$S^{++}$ 6312 Å',color="orange",/overplot);,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)

p7=plot(rho_ceq,mean(reform(B_sim(*,*,3)),dimension = 1),NAME='$S^+$ 4069 Å',color="light Green",/overplot);,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)
p8=plot(rho_ceq,mean(reform(B_sim(*,*,4)),dimension = 1),NAME='$S^+$ 4076 Å',color="dark green",/overplot);,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)


leg2 = LEGEND(TARGET=[p3,p2,p6,p8,p7,p1,p5,p4], POSITION=[7.1,400.], $
  /DATA, /AUTO_TEXT_COLOR)

  leg2.save,'0to400rayleigh_DAWN_no_offset_1.2xions_duskcoaddedfits_mean_interped_each_frame_all_apo_lines.png',resolution=300




stop
openr,1,'yptsi_8_APO_LINES_optical_vary_ne_te_lookuptable_61x57x8_CHIANTI_10.1.dat';'yptsi_8_APO_LINES_optical_vary_ne_te_lookuptable_61x57x8_CHIANTI801.dat'
yptsi_in=dblarr(61,57,8) ; ne x tec x discrete wavelength centers of each emission line
readf,1,yptsi_in
close,1

nel = [dindgen(9) + 1d, 10d*dindgen(9) + 10d, 100d*dindgen(4) + 100d, 250d*dindgen(39) + 500d]

Tec = [0.1d*dindgen(9) + 0.1d, 0.5d*dindgen(17)+1d,dindgen(10)+10d,5d*dindgen(8) + 20d, 20d*dindgen(7) + 60d, 80d*dindgen(6) + 200d]


nelhalf = [replicate(0.5d,9), replicate(5d,9), replicate(50d,4), replicate(125d,39)]

Techalf = [replicate(0.05d,9), replicate(0.25d,17),replicate(0.5d,10),replicate(2.5d,8), replicate(10d,7), replicate(40d,6)]



n_ne = n_elements(nel)

n_te = n_elements(tec)

nel2dto1d =dblarr(n_ne*n_te)
tec2dto1d =dblarr(n_ne*n_te)
yptsi_in_2dto1d = dblarr(n_ne*n_te,8)

;2nd ->1d mapping for griddata interp
k=long64(-1)
for i=0, 61 -1 do begin
  for j=0, 57 -1 do begin
    k = k + long64(1)
    nel2dto1d[k] = nel[i]
    tec2dto1d[k] = tec[j]

    yptsi_in_2dto1d[k,*] = yptsi_in[i,j,*]

    ;print,yptsi_in_2dto1d[k,*]
    ;stop
  endfor
endfor




restore,'DAWN_varyall_layers_after_func_form_only3sp_lines_simult_nec_nsp_fit_11params_int_tab_whole_side_simult_steffl_Tec_then_fit_ns2p_and_nop_7_params.sav',/verbose

TE_FIXEDD = TE_FIXED
NEC_FIXED_FOUNDD = NEC_FIXED_FOUND
NSP_FIXED_FOUNDD = NSP_FIXED_FOUND 

;stop

restore,FILENAME = 'DAWN_funcform_only3sp_lines_simult_nec_nsp_fit_11params_int_tab_whole_side_simult_steffl_Tec_then_fit_ns2p_and_nop_7_params.sav',/verbose

model = given_number_LOS_dens_temp_for_each_layer_simulate_rayleighs_dawn_or_dusk_all_at_once_8APO_lines_fixed_proper_nintegration_now(te_fixedD, nec_fixed_foundD,1.2d*nsp_fixed_foundD,1.2d*ns2p_FUNCFORM_TEST,1.2d*nop_FUNCFORM_TEST, max(x_grid),min(x_grid), 0.04D,yptsi_in_2dto1d,nel2dto1d,tec2dto1d)

;model = given_number_LOS_dens_temp_for_each_layer_simulate_rayleighs_dawn_or_dusk_all_at_once_8APO_lines_fixed_proper_nintegration_now(te_fixed, nec_fixed_found,nsp_fixed_found,ns2p_fixed_found,nop_fixed_found, max(x_grid),min(x_grid), x_step_size,yptsi_in_2dto1d,nel2dto1d,tec2dto1d)


x_grid = reverse(x_grid)

p1=plot(x_grid,model(*,2)/interpol(mean(reform(B_sim(*,*,2)),dimension = 1),rho_ceq,-x_grid),NAME='$O^+$ 3729 Å',color="blue",xtitle='Dawn $\rho_c$ ($R_J$)',ytitle='Ratio',title='Model Aligned/(Model run through process), same densities in eq.',yrange=[1,2])

p2=plot(x_grid,model(*,6)/interpol(mean(reform(B_sim(*,*,6)),dimension = 1),rho_ceq,-x_grid),NAME='$S^+$ 6716 Å',color="red",/overplot);,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)
p3=plot(x_grid,model(*,7)/interpol(mean(reform(B_sim(*,*,7)),dimension = 1),rho_ceq,-x_grid),NAME='$S^+$ 6731 Å',color="brown",/overplot);,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)

p4=plot(x_grid,model(*,0)/interpol(mean(reform(B_sim(*,*,0)),dimension = 1),rho_ceq,-x_grid),NAME='$S^{++}$ 3722 Å',color="pink",/overplot);,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)
p5=plot(x_grid,model(*,1)/interpol(mean(reform(B_sim(*,*,1)),dimension = 1),rho_ceq,-x_grid),NAME='$O^+$ 3726 Å',color="purple",/overplot);,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)

p6=plot(x_grid,model(*,5)/interpol(mean(reform(B_sim(*,*,5)),dimension = 1),rho_ceq,-x_grid),NAME='$S^{++}$ 6312 Å',color="orange",/overplot);,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)

p7=plot(x_grid,model(*,3)/interpol(mean(reform(B_sim(*,*,3)),dimension = 1),rho_ceq,-x_grid),NAME='$S^+$ 4069 Å',color="light Green",/overplot);,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)
p8=plot(x_grid,model(*,4)/interpol(mean(reform(B_sim(*,*,4)),dimension = 1),rho_ceq,-x_grid),NAME='$S^+$ 4076 Å',color="dark green",/overplot);,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)

print,(model(*,5)/interpol(mean(reform(B_sim(*,*,5)),dimension = 1),rho_ceq,-x_grid))/(model(*,0)/interpol(mean(reform(B_sim(*,*,0)),dimension = 1),rho_ceq,-x_grid))
leg2 = LEGEND(TARGET=[p3,p2,p6,p8,p7,p1,p5,p4], POSITION=[5.5,2], $
  /DATA, /AUTO_TEXT_COLOR)



;leg2.save,'ratio_1.2x_DAWN_model_aligned_over_process.png',resolution=300


p1=plot(x_grid,interpol(mean(reform(B_sim(*,*,2)),dimension = 1),rho_ceq,-x_grid),NAME='$O^+$ 3729 Å',color="blue",xtitle='Dawn $\rho_c$ ($R_J$)',ytitle='Ratio',title='Model Aligned/(Model run through process), same densities in eq.',yrange=[0,800])
p11 = plot(x_grid,model(*,2),color="blue",/overplot)
  p2=plot(x_grid,interpol(mean(reform(B_sim(*,*,6)),dimension = 1),rho_ceq,-x_grid),NAME='$S^+$ 6716 Å',color="red",/overplot);,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)
p22 = plot(x_grid,model(*,6),color="red",/overplot)
p3=plot(x_grid,interpol(mean(reform(B_sim(*,*,7)),dimension = 1),rho_ceq,-x_grid),NAME='$S^+$ 6731 Å',color="brown",/overplot);,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)
p33 = plot(x_grid,model(*,7),color="brown",/overplot)
p4=plot(x_grid,interpol(mean(reform(B_sim(*,*,0)),dimension = 1),rho_ceq,-x_grid),NAME='$S^{++}$ 3722 Å',color="pink",/overplot);,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)
p44 = plot(x_grid,model(*,0),color="pink",/overplot)
p5=plot(x_grid,interpol(mean(reform(B_sim(*,*,1)),dimension = 1),rho_ceq,-x_grid),NAME='$O^+$ 3726 Å',color="purple",/overplot);,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)
p55 = plot(x_grid,model(*,1),color="purple",/overplot)
p6=plot(x_grid,interpol(mean(reform(B_sim(*,*,5)),dimension = 1),rho_ceq,-x_grid),NAME='$S^{++}$ 6312 Å',color="orange",/overplot);,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)
p66 = plot(x_grid,model(*,5),color="orange",/overplot)
p7=plot(x_grid,interpol(mean(reform(B_sim(*,*,3)),dimension = 1),rho_ceq,-x_grid),NAME='$S^+$ 4069 Å',color="light Green",/overplot);,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)
p77 = plot(x_grid,model(*,3),color="light green",/overplot)
p8=plot(x_grid,interpol(mean(reform(B_sim(*,*,4)),dimension = 1),rho_ceq,-x_grid),NAME='$S^+$ 4076 Å',color="dark green",/overplot);,SYM_INCREMENT=2,ERRORBAR_CAPSIZE=0.01)
p88 = plot(x_grid,model(*,4),color="dark green",/overplot)

leg2 = LEGEND(TARGET=[p3,p2,p6,p8,p7,p1,p5,p4], POSITION=[7.1,800], $
  /DATA, /AUTO_TEXT_COLOR)

  factor_corrections = dblarr(60,3)

  factor_corrections(*,0) = (model(*,6)/interpol(mean(reform(B_sim(*,*,6)),dimension = 1),rho_ceq,x_grid) + model(*,7)/interpol(mean(reform(B_sim(*,*,7)),dimension = 1),rho_ceq,x_grid))/2d
  factor_corrections(*,1) = model(*,5)/interpol(mean(reform(B_sim(*,*,5)),dimension = 1),rho_ceq,x_grid)
  factor_corrections(*,2) = (model(*,1)/interpol(mean(reform(B_sim(*,*,1)),dimension = 1),rho_ceq,x_grid) + model(*,2)/interpol(mean(reform(B_sim(*,*,2)),dimension = 1),rho_ceq,x_grid))/2d
  write_csv,'factor_corrections_dawn.csv',factor_corrections




stop

for k= 0, nframes -1 do begin
  p_out(k) = plot(rho_ceq,reform(B(k,*)),xtitle='$\rho_{ceq}$',ytitle='Rayleighs',title='Linear Interp., LH Dusk $\lambda_{III}$ = ' + string(sigfig((360d - CMLdeg(k)) - 90d,4)) + '$\deg$',xrange=[4.5,7.7])
endfor

for q=0, nframes -2 do p_out[q].Save, 'LH_labeled_dusk_only_4.5-7.7_linear_interped_APO_6731_simulation_convolved_for_seeing_and_avged_over_each_bin_using_o4_constant_ceq.pdf', /APPEND

p_out[nframes -1].Save, 'LH_labeled_dusk_only_4.5-7.7_linear_interped_APO_6731_simulation_convolved_for_seeing_and_avged_over_each_bin_using_o4_constant_ceq.pdf', /APPEND, /CLOSE

stop

end