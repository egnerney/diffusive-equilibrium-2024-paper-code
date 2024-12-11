Pro plot_emission_rate_ratios_ch10_noIONRECandnorrec


  openr,1,'yptsi_550-2100_Angstroms_vary_ne_te_lookuptable_61x57x4153_CHIANTI10_1_noprotons_NOIONREC_NO_RREC.dat'
  yptsi_in=dblarr(61,57,4153) ; ne x tec x discrete wavelength centers of each emission line
  readf,1,yptsi_in
  close,1



  openr,1,'xwavi_550-2100_Angstroms_CHIANTI10_1_noprotons_NOIONREC_NO_RREC.dat'
  xwavi_in=dblarr(4153)
  readf,1,xwavi_in
  close,1



  openr,1,'yname_550-2100_Angstroms_CHIANTI10_1_noprotons_NOIONREC_NO_RREC.dat'
  yname_in=strarr(4153)
  readf,1,yname_in
  close,1

  sp_idxs = where(yname_in eq '"S II"')
  s2p_idxs = where(yname_in eq '"S III"')
  s3p_idxs = where(yname_in eq '"S IV"')
  s4p_idxs = where(yname_in eq '"S V"')
  op_idxs = where(yname_in eq '"O II"')
  o2p_idxs = where(yname_in eq '"O III"')

;680ish feature Csssini UVIS
;wavemax_want = 689d

;wavemin_want = 673d;

;1198 is a combo of the 1191 and 1201 features... I am using UVIS FUV channels 98:115, so 1191.65-1204.89A.
;1191.65-1204.89A.
;wavemax_want = 1204.89d; 

;wavemin_want = 1191.65d;

;What I am calling 1725 is a combo of the 1713 and 1729 features, UVIS channels 760:792, 1707.97-1732.92A.

wavemax_want = 1732.92d;

wavemin_want = 1707.97d;



  ;idxxx = where((xwavi_in(s2p_idxs)  le wavemax_want) and (xwavi_in(s2p_idxs)  ge wavemin_want) )
  ;xwavi_in =xwavi_in(s2p_idxs)

 ; yptsi_in_temp = reform(yptsi_in(28,17,s2p_idxs))
 ; print,xwavi_in(idxxx)
 ; print, yptsi_in_temp(idxxx);/max(yptsi_in_temp(idxxx))
  
 
  nel = [dindgen(9) + 1d, 10d*dindgen(9) + 10d, 100d*dindgen(4) + 100d, 250d*dindgen(39) + 500d]

  Tec = [0.1d*dindgen(9) + 0.1d, 0.5d*dindgen(17)+1d,dindgen(10)+10d,5d*dindgen(8) + 20d, 20d*dindgen(7) + 60d, 80d*dindgen(6) + 200d]

  n_ne = n_elements(nel)

  n_te = n_elements(tec)
  
  
  
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

  wavemax_want = 689d

  wavemin_want = 673d;

  
  idxxx_680 = where((xwavi_in(s2p_idxs)  le wavemax_want) and (xwavi_in(s2p_idxs)  ge wavemin_want) )
  
  wavemax_want = 1204.89d;

  wavemin_want = 1191.65d;


  idxxx_1198 = where((xwavi_in(s2p_idxs)  le wavemax_want) and (xwavi_in(s2p_idxs)  ge wavemin_want) )

   wavemax_want = 1732.92d;

  wavemin_want = 1707.97d;

  idxxx_1725 = where((xwavi_in(s2p_idxs)  le wavemax_want) and (xwavi_in(s2p_idxs)  ge wavemin_want) )

ratio_vs_tec_1 = dblarr(n_te) ;680Å/1725Å
ratio_vs_tec_2 = dblarr(n_te) ;1198Å/1725Å 
print,  idxxx_680
print,  idxxx_1198
print,  idxxx_1725
stop
  
  for iii = 0, n_te - 1 do begin
    print,iii
    
    yptsi_in_temp = reform(yptsi_in(28,iii,s2p_idxs))
    
    ratio_vs_tec_1(iii) =     Total(yptsi_in_temp(idxxx_680),/double)/total(yptsi_in_temp(idxxx_1725),/double)
    
    ratio_vs_tec_2(iii) =     total(yptsi_in_temp(idxxx_1198),/double)/total(yptsi_in_temp(idxxx_1725),/double)
    
  endfor
  
  p1=plot(Tec,ratio_vs_tec_1,xtitle='$T_{ec}$ (eV)',ytitle='680Å/1725Å Emission rate ratio',title='Single Maxwellian, CHIANTI 10.1, Only electron Impact Excitation',xlog=0,ylog=1,xrange=[1,10],yrange=[0.01,4])
  p2=plot(Tec,ratio_vs_tec_2,xtitle='$T_{ec}$ (eV)',ytitle='1198Å/1725Å Emission rate ratio',title='Single Maxwellian, CHIANTI 10.1, Only electron Impact Excitation',xlog=0,ylog=1,xrange=[1,10],yrange=[0.01,2])

p1.save,'noionrecandnorrec680Å_to_1725Å_Emission_rate_ratio_for_cassini_UVIS_single_maxwellian_CHIANTI10.1_only_electron_impact_excitation_nec=2000cm-3_tec_1-10eV.png',resolution=300
p2.save,'noionrecandnorrec1198Å_to_1725Å_Emission_rate_ratio_for_cassini_UVIS_single_maxwellian_CHIANTI10.1_only_electron_impact_excitation_nec=2000cm-3_tec_1-10eV.png',resolution=300

  write_csv,'noionrecandnorrec680_to_1725ANG_ratio_vs_tec_CHIANTI10.1_maxwellian_nec=2000cm-3.csv',ratio_vs_tec_1
  write_csv,'noionrecandnorrec1198_to_1725ANG_ratio_vs_tec_CHIANTI10.1_maxwellian_nec=2000cm-3.csv',ratio_vs_tec_2
  write_csv,'tec_vals_for_Cassini_UVIS_line_ratios.csv',Tec

stop
end