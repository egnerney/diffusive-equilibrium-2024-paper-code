Pro plot_emission_rate_ratios_ch10_noIONRECandnorrec_for_1356oxygen_emissivity_vs_tec


  openr,1,'neutrals_yptsi_550-2100_Angstroms_vary_ne_te_lookuptable_61x57x7_CHIANTI10_1_noprotons_NOIONREC_NO_RREC.dat'
  yptsi_in=dblarr(61,57,7) ; ne x tec x discrete wavelength centers of each emission line
  readf,1,yptsi_in
  close,1



  openr,1,'neutrals_xwavi_550-2100_Angstroms_7elements_CHIANTI10_1_noprotons_NOIONREC_NO_RREC.dat'
  xwavi_in=dblarr(7)
  readf,1,xwavi_in
  close,1



  openr,1,'neutrals_yname_550-2100_Angstroms_7elements_CHIANTI10_1_noprotons_NOIONREC_NO_RREC.dat'
  yname_in=strarr(7)
  readf,1,yname_in
  close,1

  s_idxs = where(yname_in eq '"S I"')
  o_idxs = where(yname_in eq '"O I"')

wavemax_want = 1314d; 

wavemin_want = 1294d;



  idxxx = where((xwavi_in(o_idxs)  le wavemax_want) and (xwavi_in(o_idxs)  ge wavemin_want) )
  xwavi_in =xwavi_in(o_idxs)

  yptsi_in_temp = reform(yptsi_in(28,17,o_idxs))
  print,xwavi_in(idxxx)
  print, yptsi_in_temp(idxxx);/max(yptsi_in_temp(idxxx))
  print, yptsi_in_temp(idxxx)/max(yptsi_in_temp(idxxx))
  
  print, yptsi_in_temp(idxxx)/total(yptsi_in_temp(idxxx))
  print,total(xwavi_in(idxxx)*(yptsi_in_temp(idxxx)/total(yptsi_in_temp(idxxx))))

  
  
  
  wavemax_want = 1366d;

  wavemin_want = 1346d;



  idxxx = where((xwavi_in(o_idxs)  le wavemax_want) and (xwavi_in(o_idxs)  ge wavemin_want) )
  xwavi_in =xwavi_in(o_idxs)

  yptsi_in_temp = reform(yptsi_in(28,17,o_idxs))
  print,xwavi_in(idxxx)
  print, yptsi_in_temp(idxxx);/max(yptsi_in_temp(idxxx))
  print, yptsi_in_temp(idxxx)/max(yptsi_in_temp(idxxx))

  print, yptsi_in_temp(idxxx)/total(yptsi_in_temp(idxxx))
 print,total(xwavi_in(idxxx)*(yptsi_in_temp(idxxx)/total(yptsi_in_temp(idxxx))))
 
  nel = [dindgen(9) + 1d, 10d*dindgen(9) + 10d, 100d*dindgen(4) + 100d, 250d*dindgen(39) + 500d]

  Tec = [0.1d*dindgen(9) + 0.1d, 0.5d*dindgen(17)+1d,dindgen(10)+10d,5d*dindgen(8) + 20d, 20d*dindgen(7) + 60d, 80d*dindgen(6) + 200d]

  n_ne = n_elements(nel)

  n_te = n_elements(tec)
  
  stop
  
  ;680ish feature Csssini UVIS
  ;wavemax_want = 689d

  ;wavemin_want = 673d;


  ;1198 is a combo of the 1191 and 1201 features... I am using UVIS FUV channels 98:115, so 1191.65-1204.89A.
  ;1191.65-1204.89A.
  ;wavemax_want = 1204.89d;

  ;wavemin_want = 1191.65d;

  ;What I am calling 1725 is a combo of the 1713 and 1729 features, UVIS channels 760:792, 1707.97-1732.92A.

  ;wavemax_want = 1732.92d;

  ;wavemin_want = 1707.97d;

  wavemax_want = 1314d

  wavemin_want = 1294d;

  
  idxxx_1304 = where((xwavi_in(o_idxs)  le wavemax_want) and (xwavi_in(o_idxs)  ge wavemin_want) )
  
  wavemax_want = 1366d;

  wavemin_want = 1346d;

  idxxx_1356 = where((xwavi_in(o_idxs)  le wavemax_want) and (xwavi_in(o_idxs)  ge wavemin_want) )

  

ratio_vs_tec_1 = dblarr(n_te) ;1304Å/1356Å
emiss_vs_tec_1 = dblarr(n_te) ;1304
emiss_vs_tec_2 = dblarr(n_te) ;1356
print,  idxxx_1304
print,  idxxx_1356
stop
  
  for iii = 0, n_te - 1 do begin
    print,iii
    
    yptsi_in_temp = reform(yptsi_in(28,iii,o_idxs))
    
    ratio_vs_tec_1(iii) =     Total(yptsi_in_temp(idxxx_1304),/double)/total(yptsi_in_temp(idxxx_1356),/double)
    emiss_vs_tec_1(iii) =     Total(yptsi_in_temp(idxxx_1304),/double)
 emiss_vs_tec_2(iii) =     Total(yptsi_in_temp(idxxx_1356),/double)
    
  endfor
  
  p1=plot(Tec[9:26],ratio_vs_tec_1[9:26],xtitle='$T_{ec}$ (eV)',ytitle='130.4nm/135.6nm Emission rate ratio',title='Single Maxwellian, CHIANTI 10.1, Only electron Impact Excitation',xlog=0,ylog=1,xrange=[1,10])
  p2=plot(Tec[9:26],emiss_vs_tec_1[9:26],xtitle='$T_{ec}$ (eV)',ytitle='130.4 nm Volume Emission Rate/Oxygen Density (Photons/s)',title='Single Maxwellian, CHIANTI 10.1, Only electron Impact Excitation',xlog=0,ylog=1,xrange=[1,10])
  p3=plot(Tec[9:26],emiss_vs_tec_2[9:26],xtitle='$T_{ec}$ (eV)',ytitle='135.6 nm Volume Emission Rate/Oxygen Density (Photons/s)',title='Single Maxwellian, CHIANTI 10.1, Only electron Impact Excitation',xlog=0,ylog=1,xrange=[1,10])

p1.save,'noionrecandnorrec1304Å_to_1356Å_Emission_rate_ratio_single_maxwellian_CHIANTI10.1_only_electron_impact_excitation_nec=2000cm-3_tec_1-10eV.png',resolution=300
p2.save,'noionrecandnorrec1304Å_volume_emiss_per_density_single_maxwellian_CHIANTI10.1_only_electron_impact_excitation_nec=2000cm-3_tec_1-10eV.png',resolution=300
p3.save,'noionrecandnorrec1356Å_volume_emiss_per_density_single_maxwellian_CHIANTI10.1_only_electron_impact_excitation_nec=2000cm-3_tec_1-10eV.png',resolution=300

  write_csv,'noionrecandnorrec1304Å_to_1356Å_ratio_vs_tec_CHIANTI10.1_maxwellian_nec=2000cm-3.csv',ratio_vs_tec_1
  write_csv,'noionrecandnorrec1304Å_volume_emiss_per_density_vs_tec_CHIANTI10.1_maxwellian_nec=2000cm-3.csv',emiss_vs_tec_1
  write_csv,'noionrecandnorrec1356Å_volume_emiss_per_density_vs_tec_CHIANTI10.1_maxwellian_nec=2000cm-3.csv',emiss_vs_tec_2
  ;write_csv,'noionrecandnorrec1198_to_1725ANG_ratio_vs_tec_CHIANTI10.1_maxwellian_nec=2000cm-3.csv',ratio_vs_tec_2
  
stop
end