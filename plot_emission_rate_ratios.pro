Pro plot_emission_rate_ratios


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



  idxxx = where((xwavi_in(s2p_idxs)  le wavemax_want) and (xwavi_in(s2p_idxs)  ge wavemin_want) )
  xwavi_in =xwavi_in(s2p_idxs)

  yptsi_in_temp = reform(yptsi_in(28,17,s2p_idxs))
  print,xwavi_in(idxxx)
  print, yptsi_in_temp(idxxx);/max(yptsi_in_temp(idxxx))
  
 
  nel = [dindgen(9) + 1d, 10d*dindgen(9) + 10d, 100d*dindgen(4) + 100d, 250d*dindgen(39) + 500d]

  Tec = [0.1d*dindgen(9) + 0.1d, 0.5d*dindgen(17)+1d,dindgen(10)+10d,5d*dindgen(8) + 20d, 20d*dindgen(7) + 60d, 80d*dindgen(6) + 200d]

  n_ne = n_elements(nel)

  n_te = n_elements(tec)
;stop

  openr,1,'yptsi_550-2100_Angstroms_vary_ne_te_lookuptable_61x57x4153_CHIANTI10_1_noprotons_withIONRECandrrec.dat'
  yptsi_in_chi10_np_withIONRECandrrec=dblarr(61,57,4153) ; ne x tec x discrete wavelength centers of each emission line
  readf,1,yptsi_in_chi10_np_withIONRECandrrec
  close,1



  openr,1,'xwavi_550-2100_Angstroms_CHIANTI10_1_noprotons_withIONRECandrrec.dat'
  xwavi_in_chi10_np_withIONRECandrrec=dblarr(4153)
  readf,1,xwavi_in_chi10_np_withIONRECandrrec
  close,1



  openr,1,'yname_550-2100_Angstroms_CHIANTI10_1_noprotons_withIONRECandrrec.dat'
  yname_in_chi10_np_withIONRECandrrec=strarr(4153)
  readf,1,yname_in_chi10_np_withIONRECandrrec
  close,1

  sp_idxs_chi10_np_withIONRECandrrec = where(yname_in_chi10_np_withIONRECandrrec eq '"S II"')
  s2p_idxs_chi10_np_withIONRECandrrec = where(yname_in_chi10_np_withIONRECandrrec eq '"S III"')
  s3p_idxs_chi10_np_withIONRECandrrec = where(yname_in_chi10_np_withIONRECandrrec eq '"S IV"')
  s4p_idxs_chi10_np_withIONRECandrrec = where(yname_in_chi10_np_withIONRECandrrec eq '"S V"')
  op_idxs_chi10_np_withIONRECandrrec = where(yname_in_chi10_np_withIONRECandrrec eq '"O II"')
  o2p_idxs_chi10_np_withIONRECandrrec = where(yname_in_chi10_np_withIONRECandrrec eq '"O III"')

;stop


  idxxx_chi10_np_withIONRECandrrec = where((xwavi_in_chi10_np_withIONRECandrrec(s2p_idxs_chi10_np_withIONRECandrrec)  le wavemax_want) and (xwavi_in_chi10_np_withIONRECandrrec(s2p_idxs_chi10_np_withIONRECandrrec)  ge wavemin_want) )
  xwavi_in_chi10_np_withIONRECandrrec =xwavi_in_chi10_np_withIONRECandrrec(s2p_idxs_chi10_np_withIONRECandrrec)

  yptsi_in_temp_chi10_np_withIONRECandrrec = reform(yptsi_in_chi10_np_withIONRECandrrec(28,17,s2p_idxs_chi10_np_withIONRECandrrec))
  print,xwavi_in_chi10_np_withIONRECandrrec(idxxx_chi10_np_withIONRECandrrec)
  print, yptsi_in_temp_chi10_np_withIONRECandrrec(idxxx_chi10_np_withIONRECandrrec);/max(yptsi_in_temp_chi10_np_withIONRECandrrec(idxxx_chi10_np_withIONRECandrrec))



print, (yptsi_in_temp_chi10_np_withIONRECandrrec(idxxx_chi10_np_withIONRECandrrec) - yptsi_in_temp(idxxx))/yptsi_in_temp_chi10_np_withIONRECandrrec(idxxx_chi10_np_withIONRECandrrec)

print,xwavi_in(idxxx)-xwavi_in_chi10_np_withIONRECandrrec(idxxx_chi10_np_withIONRECandrrec)

stop
end