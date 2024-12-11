pro read_in_plot_APO_chianti10_1_emission_ratios



  nel = [dindgen(9) + 1d, 10d*dindgen(9) + 10d, 100d*dindgen(4) + 100d, 250d*dindgen(39) + 500d]

  Tec = [0.1d*dindgen(9) + 0.1d, 0.5d*dindgen(17)+1d,dindgen(10)+10d,5d*dindgen(8) + 20d, 20d*dindgen(7) + 60d, 80d*dindgen(6) + 200d]

print,tec(17),nel(28)

nel2 = [1d,100d,500d, 1000d,1500d,2000d,2500d,3000d,3500d,4000d,6000d,8000d,10000d]

;Tec = [0.1d*dindgen(9) + 0.1d, 0.5d*dindgen(17)+1d,dindgen(10)+10d,5d*dindgen(8) + 20d, 20d*dindgen(7) + 60d, 80d*dindgen(6) + 200d]
Tec2 = [0.1,0.5d,1d,2d,3d,4d,5d,6d,8d,10d]



feh2 = [0.0001d,0.001d,0.0025d,0.005d,0.01d,0.05d,0.1d]
Teh2 = [5d*dindgen(4) + 35d, 20d*dindgen(3) + 60d, 100d*dindgen(3) + 200d]

n_ne = n_elements(nel2)

n_te = n_elements(tec2)
n_feh = n_elements(feh2)
n_teh = n_elements(teh2)


feh3 = [0.0025d,0.01d];[0.0001d,0.001d,0.0025d,0.005d,0.01d,0.05d,0.1d]
Teh3 = [45d,300d];[5d*dindgen(4) + 35d, 20d*dindgen(3) + 60d, 100d*dindgen(3) + 200d]


 ; stop

  openr,1,'yptsi_8_APO_LINES_optical_vary_ne_te_lookuptable_61x57x8_CHIANTI_10.1.dat' ;  ;'yptsi_8_APO_LINES_optical_vary_ne_te_lookuptable_61x57x8_CHIANTI801.dat'
  yptsi_in=dblarr(61,57,8) ; ne x tec x discrete wavelength centers of each emission line
  readf,1,yptsi_in
  close,1
  
  
  openr,1,'yptsi_8_APO_LINES_optical_vary_nec_tec_feh_teh_lookuptable_13x10x7x10x8_CHIANTI_10.1_low_res.dat' ;  ;'yptsi_8_APO_LINES_optical_vary_ne_te_lookuptable_61x57x8_CHIANTI801.dat'
  yptsi_in_hote=dblarr(13,10,7,10,8) ; ne x tec x discrete wavelength centers of each emission line
  readf,1,yptsi_in_hote
  close,1
  
  openr,1,'yptsi_8_APO_LINES_optical_vary_nec_tec_feh_teh_lookuptable_Tec_slice_nel=2000percc_57x2x2x8_CHIANTI_10.1.dat' ;  ;'yptsi_8_APO_LINES_optical_vary_ne_te_lookuptable_61x57x8_CHIANTI801.dat'
  yptsi_in_hote_vstec_slice=dblarr(57,2,2,8) ; Ted x feh3 x tec3 discrete wavelength centers of each emission line
  readf,1,yptsi_in_hote_vstec_slice
  close,1
  
  openr,1,'yptsi_8_APO_LINES_optical_vary_nec_tec_feh_teh_lookuptable_nel_slice_Tec=5eV_61x2x2x8_CHIANTI_10.1.dat' ;  ;'yptsi_8_APO_LINES_optical_vary_ne_te_lookuptable_61x57x8_CHIANTI801.dat'
  yptsi_in_hote_vsnel_slice=dblarr(61,2,2,8) ; nel x feh3 x tec3 discrete wavelength centers of each emission line
  readf,1,yptsi_in_hote_vsnel_slice
  close,1
  
  
    
;tec(17) = 5eV


;nel(28) = 2000 cm^{-3}

left = 0.2
right = 0.06
MARGINz0 =[left,0.1,right,0.1]
MARGINz1 =[left,0.1,right,0.0]
MARGINz2 =[left,0.25,right,0.0]

marginz = []   ;[left, bottom, right, top]

;p1  = plot(tec,yptsi_in(28,*,6)/yptsi_in(28,*,7),xtitle='$T_{ec}$ (eV)',ytitle='Emissivity Ratio',title='$n_{ec}$ = 2000 $cm^{-3}$',name='$S^{+}$ 6716Å to 6731Å ',yrange=[0.4,1.6],/xlog,layout=[2,1,1],margin=marginz0)

;p2  = plot(tec,yptsi_in(28,*,2)/yptsi_in(28,*,1),name='$O^{+}$ 3729Å to 3726Å',/overplot,color = 'firebrick',yrange=[0.4,1.6],/xlog,margin=marginz0)

;leg = LEGEND(TARGET=[p1,p2], POSITION=[300,1.4], $
;  /DATA, /AUTO_TEXT_COLOR)

;leg.save,'S+_andO+_doublet_chianti10.1_vs_Tec_at_nec=2000_per_cc.png',resolution=300


;p3  = plot(nel,yptsi_in(*,17,6)/yptsi_in(*,17,7),xtitle='$n_{ec}$ ($cm^{-3}$)',title='$T_{ec}$ = 5 eV',name='$S^{+}$ 6716Å to 6731Å ',/xlog,layout=[2,1,2],/current,yrange=[0.4,1.6],margin=marginz0)
  ;p4  = plot(nel,yptsi_in(*,17,2)/yptsi_in(*,17,1),name='$O^{+}$ 3729Å to 3726Å',/overplot,color = 'firebrick',/xlog,yrange=[0.4,1.6],margin=marginz0)


  p1  = plot(tec,yptsi_in(28,*,6)/nel(28),xtitle='$T_{ec}$ (eV)',ytitle='Emissivity',title='$n_{ec}$ = 2000 $cm^{-3}$',/xlog,/ylog)

    p2  = plot(tec,yptsi_in_hote_vstec_slice(*,0,0,6)/nel(28),/xlog,/ylog,/overplot,color='red')
  
  p3  = plot(nel,yptsi_in(*,17,6)/nel,xtitle='$n_{el}$ ($cm^{-3}$)',ytitle='Emissivity',title='$T_{ec}$ = 5 eV',/xlog,/ylog)

    p4  = plot(nel,yptsi_in_hote_vsnel_slice(*,0,0,6)/nel,/xlog,/ylog,/overplot,color='red')



stop
p1  = plot(tec,yptsi_in(28,*,7)/nel(28),xtitle='$T_{ec}$ (eV)',ytitle='Emissivity',title='$n_{ec}$ = 2000 $cm^{-3}$',/xlog,/ylog)

  p2  = plot(tec,yptsi_in_hote_vstec_slice(*,0,0,7)/nel(28),/xlog,/ylog,/overplot,color='red')


  p3  = plot(nel,yptsi_in(*,17,7)/nel,xtitle='$n_{el}$ ($cm^{-3}$)',ytitle='Emissivity',title='$n_{ec}$ =5 eV',/xlog,/ylog)

    p4  = plot(nel,yptsi_in_hote_vsnel_slice(*,0,0,7)/nel,/xlog,/ylog,/overplot,color='red')


stop

  p1  = plot(tec,reform(yptsi_in(28,*,5)/yptsi_in(28,*,0)),xtitle='$T_{ec}$ (eV)',ytitle='Emissivity Ratio',title='$n_{ec}$ = 2000 $cm^{-3}$, $S^{++}$ 6312Å to 3722Å ',name='$S^{++}$ 6312Å to 3722Å ',/xlog,layout=[2,1,1],yrange=[2,4],margin=marginz0)



  p333  = plot(nel,reform(yptsi_in(*,17,5)/yptsi_in(*,17,0)),xtitle='$n_{ec}$ ($cm^{-3}$)',title='$T_{ec}$ = 5 eV, $S^{++}$ 6312Å to 3722Å ',name='$S^{++}$ 6312Å to 3722Å ',/xlog,layout=[2,1,2],/current,yrange=[2,4],margin=marginz0)

;stop

  p11111  = plot(tec,reform(yptsi_in(28,*,7)/yptsi_in(28,*,4)),xtitle='$T_{ec}$ (eV)',ytitle='Emissivity Ratio',title='$n_{ec}$ = 2000 $cm^{-3}$, $S^{+}$ 6731Å to 4076Å ',name='$S^{+}$ 6312Å to 4076Å ',/xlog,layout=[2,1,1],yrange=[1,50],margin=marginz0)



  p3333  = plot(nel,reform(yptsi_in(*,17,7)/yptsi_in(*,17,4)),xtitle='$n_{ec}$ ($cm^{-3}$)',title='$T_{ec}$ = 5 eV, $S^{+}$ 6731Å to 4076Å ',name='$S^{+}$ 6731Å to 4076Å ',/xlog,layout=[2,1,2],/current,margin=marginz0)

;stop

  p01110  = plot(tec,reform(yptsi_in(28,*,0)),xtitle='$T_{ec}$ (eV)',ytitle='Emission Rate Coefficients (Photons/s)',title='$n_{ec}$ = 2000 $cm^{-3}$',name='$S^{++}$ 3722 Å ',color='pink',/xlog,layout=[2,1,1],yrange=[1e-20,6e-5],/ylog,margin=marginz0)

  p01111  = plot(tec,reform(yptsi_in(28,*,1)),name='$O^{+}$ 3726 Å ',color='purple',/overplot)
  p01112  = plot(tec,reform(yptsi_in(28,*,2)),name='$O^{+}$ 3729 Å ',color='blue',/overplot)
  p01113  = plot(tec,reform(yptsi_in(28,*,3)),name='$S^{+}$ 4069 Å ',color='light Green',/overplot)
  p01114  = plot(tec,reform(yptsi_in(28,*,4)),name='$S^{+}$ 4076 Å ',color='dark Green',/overplot)
  p01115  = plot(tec,reform(yptsi_in(28,*,5)),name='$S^{++}$ 6312 Å ',color='orange',/overplot)
  p01116  = plot(tec,reform(yptsi_in(28,*,6)),name='$S^{+}$ 6716 Å ',color='red',/overplot)
  p01117  = plot(tec,reform(yptsi_in(28,*,7)),name='$S^{+}$ 6731 Å ',color='brown',/overplot)
  
  
  
  p11110  = plot(nel,reform(yptsi_in(*,17,0)),xtitle='$n_{ec}$ (eV)',ytitle='Emission Rate Coefficients (Photons/s)',title='$T_{ec}$ = 5 eV',name='$S^{++}$ 3722 Å ',color='pink',/xlog,/current,layout=[2,1,2],yrange=[6e-10,1e-4],/ylog,margin=marginz0)

    p11111  = plot(nel,reform(yptsi_in(*,17,1)),name='$O^{+}$ 3726 Å ',color='purple',/overplot)
  p11112  = plot(nel,reform(yptsi_in(*,17,2)),name='$O^{+}$ 3729 Å ',color='blue',/overplot)
  p11113  = plot(nel,reform(yptsi_in(*,17,3)),name='$S^{+}$ 4069 Å ',color='light Green',/overplot)
  p11114  = plot(nel,reform(yptsi_in(*,17,4)),name='$S^{+}$ 4076 Å ',color='dark Green',/overplot)
  p11115  = plot(nel,reform(yptsi_in(*,17,5)),name='$S^{++}$ 6312 Å ',color='orange',/overplot)
  p11116  = plot(nel,reform(yptsi_in(*,17,6)),name='$S^{+}$ 6716 Å ',color='red',/overplot)
  p11117  = plot(nel,reform(yptsi_in(*,17,7)),name='$S^{+}$ 6731 Å ',color='brown',/overplot)

  leg12345 = LEGEND(TARGET=[p01117,p01116,p01115,p01114,p01113,p01112,p01111,p01110], POSITION=[300,1e-10], $
    /DATA, /AUTO_TEXT_COLOR)



eff_emiss_coeffs_vs_tec=dblarr(57,8)

eff_emiss_coeffs_vs_nec=dblarr(61,8)


;yptsi_in=dblarr(61,57,8) 

for iiii = 0, 7 do begin
  eff_emiss_coeffs_vs_tec[*,iiii]= reform(yptsi_in(28,*,iiii)/nel(28))
  eff_emiss_coeffs_vs_nec[*,iiii] = reform(yptsi_in(*,17,iiii))/nel
endfor

eff_emiss_coeffs_vs_tec_fixed_feh_TEH=dblarr(57,8)

eff_emiss_coeffs_vs_nel_fixed_feh_TEH=dblarr(61,8)


;yptsi_in_hote_vsnel_slice=dblarr(57,2,2,8)
;yptsi_in_hote_vstec_slice=dblarr(61,2,2,8)


for iiii = 0, 7 do begin
  eff_emiss_coeffs_vs_tec_fixed_feh_TEH[*,iiii]= reform(yptsi_in_hote_vstec_slice(*,0,0,iiii)/nel(28))
  eff_emiss_coeffs_vs_nel_fixed_feh_TEH[*,iiii] = reform(yptsi_in_hote_vsnel_slice(*,0,0,iiii))/nel
endfor



;write_csv,'APO_eff_emiss_coeffs_vs_tec_57x8.csv',eff_emiss_coeffs_vs_tec

;write_csv,'APO_eff_emiss_coeffs_vs_nec_61x8.csv',eff_emiss_coeffs_vs_nec

write_csv,'APO_eff_emiss_coeffs_vs_tec_57x8_feh=0.0025_Teh=45eV.csv',eff_emiss_coeffs_vs_tec_fixed_feh_TEH

write_csv,'APO_eff_emiss_coeffs_vs_nel_61x8_feh=0.0025_Teh=45eV.csv',eff_emiss_coeffs_vs_nel_fixed_feh_TEH

eff_emiss_coeffs_vs_tec_fixed_feh_TEH=dblarr(57,8)

eff_emiss_coeffs_vs_nel_fixed_feh_TEH=dblarr(61,8)


;yptsi_in_hote_vsnel_slice=dblarr(57,2,2,8)
;yptsi_in_hote_vstec_slice=dblarr(61,2,2,8)


for iiii = 0, 7 do begin
  eff_emiss_coeffs_vs_tec_fixed_feh_TEH[*,iiii]= reform(yptsi_in_hote_vstec_slice(*,1,1,iiii)/nel(28))
  eff_emiss_coeffs_vs_nel_fixed_feh_TEH[*,iiii] = reform(yptsi_in_hote_vsnel_slice(*,1,1,iiii))/nel
endfor

write_csv,'APO_eff_emiss_coeffs_vs_tec_57x8_feh=0.01_Teh=300eV.csv',eff_emiss_coeffs_vs_tec_fixed_feh_TEH

write_csv,'APO_eff_emiss_coeffs_vs_nel_61x8_feh=0.01_Teh=300eV.csv',eff_emiss_coeffs_vs_nel_fixed_feh_TEH



;  leg12345.save,'all_apo_lines_Emission_Rate_Coefficients_Photons_per_s_vs_nec_and_tec_at_constantothers_chianti10.1_singlemaxwellians_nohote.png',resolution=300

;print,min(reform(yptsi_in(28,*,*))),max(reform(yptsi_in(28,*,*)))

;print,min(reform(yptsi_in(*,17,*))),max(reform(yptsi_in(*,17,*)))

;p4.save,'S+_andO+_doublet_chianti10.1_vs_nec_at_Tec=5eV_and_vs_Tec_at_nec=2000_percc.png',resolution=300


;write_csv,'tec_for_ratio_plots_57_elements.csv',tec
;write_csv,'chianti10.1_APO_sp_emiss_ratio_vs_tec_at_nec=2000_percc_57_elements.csv',reform(yptsi_in(28,*,6)/yptsi_in(28,*,7))
;write_csv,'chianti10.1_APO_op_emiss_ratio_vs_tec_at_nec=2000_percc_57_elements.csv',reform(yptsi_in(28,*,2)/yptsi_in(28,*,1))
;write_csv,'nec_for_ratio_plots_61_elements.csv',nel
;write_csv,'chianti10.1_APO_sp_emiss_ratio_vs_nec_at_Tec=5eV_61_elements.csv',reform(yptsi_in(*,17,6)/yptsi_in(*,17,7))
;write_csv,'chianti10.1_APO_op_emiss_ratio_vs_nec_at_Tec=5eV_61_elements.csv',reform(yptsi_in(*,17,2)/yptsi_in(*,17,1))
;write_csv,'chianti10.1_APO_sp_6716to6731ang_emiss_ratio_vs_nec_and_Tec_61x57_elements.csv',reform(yptsi_in(*,*,6)/yptsi_in(*,*,7))
;write_csv,'chianti10.1_APO_op_3729to3726ang_emiss_ratio_vs_nec_and_Tec_61x57_elements.csv',reform(yptsi_in(*,*,2)/yptsi_in(*,*,1))
;write_csv,'chianti10.1_APO_sp_4069to6731ang_emiss_ratio_vs_nec_and_Tec_61x57_elements.csv',reform(yptsi_in(*,*,3)/yptsi_in(*,*,7))
;write_csv,'chianti10.1_APO_sp_emissline2_ratio_vs_tec_at_nec=2000_percc_57_elements.csv',reform(yptsi_in(28,*,3)/yptsi_in(28,*,7))
;write_csv,'chianti10.1_APO_sp_emissline2_ratio_vs_nec_at_Tec=5eV_61_elements.csv',reform(yptsi_in(*,17,3)/yptsi_in(*,17,7))

;write_csv,'chianti10.1_APO_sp_emissline3_ratio_vs_tec_at_nec=2000_percc_57_elements.csv',reform(yptsi_in(28,*,4)/yptsi_in(28,*,7))
;write_csv,'chianti10.1_APO_sp_emissline3_ratio_vs_nec_at_Tec=5eV_61_elements.csv',reform(yptsi_in(*,17,4)/yptsi_in(*,17,7))
;write_csv,'chianti10.1_APO_sp_4076to6731ang_emiss_ratio_vs_nec_and_Tec_61x57_elements.csv',reform(yptsi_in(*,*,4)/yptsi_in(*,*,7))


;write_csv,'chianti10.1_APO_s2p_emiss_ratio6312to3722ang_vs_tec_at_nec=2000_percc_57_elements.csv',reform(yptsi_in(28,*,5)/yptsi_in(28,*,0))
;write_csv,'chianti10.1_APO_s2p_emiss_ratio6312to3722ang_vs_nec_at_Tec=5eV_61_elements.csv',reform(yptsi_in(*,17,5)/yptsi_in(*,17,0))
;write_csv,'chianti10.1_APO_s2p_6312to3722ang_emiss_ratio_vs_nec_and_Tec_61x57_elements.csv',reform(yptsi_in(*,*,5)/yptsi_in(*,*,0))


;write_csv,'chianti10.1_APO_sp_6716to6731_emiss_ratio_vs_tec_at_nel=2000_percc_feh=0.0025_Teh=45eV_lowres.csv',reform(yptsi_in_hote(5,*,2,2,6)/yptsi_in_hote(5,*,2,2,7))
;write_csv,'chianti10.1_APO_op_3729to3726_emiss_ratio_vs_tec_at_nel=2000_percc_feh=0.0025_Teh=45eV_lowres.csv',reform(yptsi_in_hote(5,*,2,2,2)/yptsi_in_hote(5,*,2,2,1))
;write_csv,'chianti10.1_APO_sp_4069to6731_emiss_ratio_vs_tec_at_nel=2000_percc_feh=0.0025_Teh=45eV_lowres.csv',reform(yptsi_in_hote(5,*,2,2,3)/yptsi_in_hote(5,*,2,2,7))

;write_csv,'chianti10.1_APO_sp_6716to6731_emiss_ratio_vs_nel_at_Tec=5eV_feh=0.0025_Teh=45eV_lowres.csv',reform(yptsi_in_hote(*,6,2,2,6)/yptsi_in_hote(*,6,2,2,7))
;write_csv,'chianti10.1_APO_op_3729to3726_emiss_ratio_vs_nel_at_Tec=5eV_feh=0.0025_Teh=45eV_lowres.csv',reform(yptsi_in_hote(*,6,2,2,2)/yptsi_in_hote(*,6,2,2,1))
;write_csv,'chianti10.1_APO_sp_4069to6731_emiss_ratio_vs_nel_at_Tec=5eV_feh=0.0025_Teh=45eV_lowres.csv',reform(yptsi_in_hote(*,6,2,2,3)/yptsi_in_hote(*,6,2,2,7))

;write_csv,'tec_for_ratio_plots_10_elements.csv',tec2

;write_csv,'nel_for_ratio_plots_13_elements.csv',nel2

;write_csv,'chianti10.1_APO_sp_6716to6731_emiss_ratio_vs_tec_at_nel=2000_percc_feh=0.1_Teh=300eV_lowres.csv',reform(yptsi_in_hote(5,*,6,8,6)/yptsi_in_hote(5,*,6,8,7))
;write_csv,'chianti10.1_APO_op_3729to3726_emiss_ratio_vs_tec_at_nel=2000_percc_feh=0.1_Teh=300eV_lowres.csv',reform(yptsi_in_hote(5,*,6,8,2)/yptsi_in_hote(5,*,6,8,1))
;write_csv,'chianti10.1_APO_sp_4069to6731_emiss_ratio_vs_tec_at_nel=2000_percc_feh=0.1_Teh=300eV_lowres.csv',reform(yptsi_in_hote(5,*,6,8,3)/yptsi_in_hote(5,*,6,8,7))

;write_csv,'chianti10.1_APO_sp_6716to6731_emiss_ratio_vs_nel_at_Tec=5eV_feh=0.1_Teh=300eV_lowres.csv',reform(yptsi_in_hote(*,6,6,8,6)/yptsi_in_hote(*,6,6,8,7))
;write_csv,'chianti10.1_APO_op_3729to3726_emiss_ratio_vs_nel_at_Tec=5eV_feh=0.1_Teh=300eV_lowres.csv',reform(yptsi_in_hote(*,6,6,8,2)/yptsi_in_hote(*,6,6,8,1))
;write_csv,'chianti10.1_APO_sp_4069to6731_emiss_ratio_vs_nel_at_Tec=5eV_feh=0.1_Teh=300eV_lowres.csv',reform(yptsi_in_hote(*,6,6,8,3)/yptsi_in_hote(*,6,6,8,7))


;write_csv,'chianti10.1_APO_s2p_6312to3722ang_emiss_ratio_vs_tec_at_nel=2000_percc_feh=0.0025_Teh=45eV_lowres.csv',reform(yptsi_in_hote(5,*,2,2,5)/yptsi_in_hote(5,*,2,2,0))
;write_csv,'chianti10.1_APO_s2p_6312to3722ang_emiss_ratio_vs_nel_at_Tec=5eV_feh=0.0025_Teh=45eV_lowres.csv',reform(yptsi_in_hote(*,6,2,2,5)/yptsi_in_hote(*,6,2,2,0))

;write_csv,'chianti10.1_APO_s2p_6312to3722ang_emiss_ratio_vs_tec_at_nel=2000_percc_feh=0.1_Teh=300eV_lowres.csv',reform(yptsi_in_hote(5,*,6,8,5)/yptsi_in_hote(5,*,6,8,0))
;write_csv,'chianti10.1_APO_s2p_6312to3722ang_emiss_ratio_vs_nel_at_Tec=5eV_feh=0.1_Teh=300eV_lowres.csv',reform(yptsi_in_hote(*,6,6,8,5)/yptsi_in_hote(*,6,6,8,0))



;write_csv,'chianti10.1_APO_sp_6716to6731_emiss_ratio_vs_tec_at_nel=2000_percc_feh=0.0025_Teh=45eV_highres61elements.csv',reform(yptsi_in_hote_vstec_slice(*,0,0,6)/yptsi_in_hote_vstec_slice(*,0,0,7))
;write_csv,'chianti10.1_APO_op_3729to3726_emiss_ratio_vs_tec_at_nel=2000_percc_feh=0.0025_Teh=45eV_highres61elements.csv',reform(yptsi_in_hote_vstec_slice(*,0,0,2)/yptsi_in_hote_vstec_slice(*,0,0,1))
;write_csv,'chianti10.1_APO_sp_4069to6731_emiss_ratio_vs_tec_at_nel=2000_percc_feh=0.0025_Teh=45eV_highres61elements.csv',reform(yptsi_in_hote_vstec_slice(*,0,0,3)/yptsi_in_hote_vstec_slice(*,0,0,7))
;write_csv,'chianti10.1_APO_s2p_6312to3722_emiss_ratio_vs_tec_at_nel=2000_percc_feh=0.0025_Teh=45eV_highres61elements.csv',reform(yptsi_in_hote_vstec_slice(*,0,0,5)/yptsi_in_hote_vstec_slice(*,0,0,0))


;write_csv,'chianti10.1_APO_sp_6716to6731_emiss_ratio_vs_nel_at_Tec=5eV_feh=0.0025_Teh=45eV_highres61elements.csv',reform(yptsi_in_hote_vsnel_slice(*,0,0,6)/yptsi_in_hote_vsnel_slice(*,0,0,7))
;write_csv,'chianti10.1_APO_op_3729to3726_emiss_ratio_vs_nel_at_Tec=5eV_feh=0.0025_Teh=45eV_highres61elements.csv',reform(yptsi_in_hote_vsnel_slice(*,0,0,2)/yptsi_in_hote_vsnel_slice(*,0,0,1))
;write_csv,'chianti10.1_APO_sp_4069to6731_emiss_ratio_vs_nel_at_Tec=5eV_feh=0.0025_Teh=45eV_highres61elements.csv',reform(yptsi_in_hote_vsnel_slice(*,0,0,3)/yptsi_in_hote_vsnel_slice(*,0,0,7))
;write_csv,'chianti10.1_APO_s2p_6312to3722_emiss_ratio_vs_nel_at_Tec=5eV_feh=0.0025_Teh=45eV_highres61elements.csv',reform(yptsi_in_hote_vsnel_slice(*,0,0,5)/yptsi_in_hote_vsnel_slice(*,0,0,0))



;write_csv,'chianti10.1_APO_sp_6716to6731_emiss_ratio_vs_tec_at_nel=2000_percc_feh=0.01_Teh=300eV_highres61elements.csv',reform(yptsi_in_hote_vstec_slice(*,1,1,6)/yptsi_in_hote_vstec_slice(*,1,1,7))
;write_csv,'chianti10.1_APO_op_3729to3726_emiss_ratio_vs_tec_at_nel=2000_percc_feh=0.01_Teh=300eV_highres61elements.csv',reform(yptsi_in_hote_vstec_slice(*,1,1,2)/yptsi_in_hote_vstec_slice(*,1,1,1))
;write_csv,'chianti10.1_APO_sp_4069to6731_emiss_ratio_vs_tec_at_nel=2000_percc_feh=0.01_Teh=300eV_highres61elements.csv',reform(yptsi_in_hote_vstec_slice(*,1,1,3)/yptsi_in_hote_vstec_slice(*,1,1,7))
;write_csv,'chianti10.1_APO_s2p_6312to3722_emiss_ratio_vs_tec_at_nel=2000_percc_feh=0.01_Teh=300eV_highres61elements.csv',reform(yptsi_in_hote_vstec_slice(*,1,1,5)/yptsi_in_hote_vstec_slice(*,1,1,0))


;write_csv,'chianti10.1_APO_sp_6716to6731_emiss_ratio_vs_nel_at_Tec=5eV_feh=0.01_Teh=300eV_highres61elements.csv',reform(yptsi_in_hote_vsnel_slice(*,1,1,6)/yptsi_in_hote_vsnel_slice(*,1,1,7))
;write_csv,'chianti10.1_APO_op_3729to3726_emiss_ratio_vs_nel_at_Tec=5eV_feh=0.01_Teh=300eV_highres61elements.csv',reform(yptsi_in_hote_vsnel_slice(*,1,1,2)/yptsi_in_hote_vsnel_slice(*,1,1,1))
;write_csv,'chianti10.1_APO_sp_4069to6731_emiss_ratio_vs_nel_at_Tec=5eV_feh=0.01_Teh=300eV_highres61elements.csv',reform(yptsi_in_hote_vsnel_slice(*,1,1,3)/yptsi_in_hote_vsnel_slice(*,1,1,7))
;write_csv,'chianti10.1_APO_s2p_6312to3722_emiss_ratio_vs_nel_at_Tec=5eV_feh=0.01_Teh=300eV_highres61elements.csv',reform(yptsi_in_hote_vsnel_slice(*,1,1,5)/yptsi_in_hote_vsnel_slice(*,1,1,0))




stop
end