Pro plot_emission_rate_ratios_ch10_noIONRECandnorrec_APO


  openr,1,'yptsi_3000-7000_Angstroms_vary_ne_te_lookuptable_61x57x2503_CHIANTI10_1_noprotons_NOIONREC_NO_RREC.dat'
  yptsi_in=dblarr(61,57,2503) ; ne x tec x discrete wavelength centers of each emission line
  readf,1,yptsi_in
  close,1



  openr,1,'xwavi_3000-7000_Angstroms_CHIANTI10_1_noprotons_NOIONREC_NO_RREC.dat'
  xwavi_in=dblarr(2503)
  readf,1,xwavi_in
  close,1



  openr,1,'yname_3000-7000_Angstroms_CHIANTI10_1_noprotons_NOIONREC_NO_RREC.dat'
  yname_in=strarr(2503)
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

wavemax_want = 6320d;1732.92d;

wavemin_want = 6300d ;1707.97d;



  idxxx = where((xwavi_in(s2p_idxs)  le wavemax_want) and (xwavi_in(s2p_idxs)  ge wavemin_want) )
 xwavi_inn = xwavi_in(s2p_idxs)

  yptsi_in_temp = reform(yptsi_in(28,17,s2p_idxs))
  print,xwavi_inn(idxxx)
  print, yptsi_in_temp(idxxx)/max(yptsi_in_temp(idxxx))
  
  STOP
  
  
  openr,1,'yptsi_8_APO_LINES_optical_vary_ne_te_lookuptable_61x57x8_CHIANTI801.dat'
  yptsi_in_APO_chianti8=dblarr(61,57,8) ; ne x tec x discrete wavelength centers of each emission line
  readf,1,yptsi_in_APO_chianti8
  close,1
  
   ;yptsi_in_temp_APO_chianti8 = reform(yptsi_in_APO_chianti8(28,17,6:7))
  
 ; print, yptsi_in_temp_APO_chianti8
  
  

  

  nel = [dindgen(9) + 1d, 10d*dindgen(9) + 10d, 100d*dindgen(4) + 100d, 250d*dindgen(39) + 500d]

  Tec = [0.1d*dindgen(9) + 0.1d, 0.5d*dindgen(17)+1d,dindgen(10)+10d,5d*dindgen(8) + 20d, 20d*dindgen(7) + 60d, 80d*dindgen(6) + 200d]

nel_halfwidth = [replicate(1d,9)/2d,replicate(10d,9)/2d,replicate(100d,4)/2d,replicate(250d,39)/2d]
tec_halfwidth = [replicate(0.1d,9)/2d,replicate(0.5d,17)/2d,replicate(1d,10)/2d,replicate(5d,8)/2d,replicate(20d,7)/2d,replicate(80d,6)/2d]


  n_ne = n_elements(nel)

  n_te = n_elements(tec)
  
diff_6731= dblarr(n_ne,n_te)
diff_6716= dblarr(n_ne,n_te)

  
  for iii = 0, n_ne - 1 do begin
      print,iii
    for jjj = 0, n_te - 1 do begin
  
    
    yptsi_in_temp = reform(yptsi_in(iii,jjj,s2p_idxs))
    
    yptsi_in_temp = yptsi_in_temp(idxxx)
    
   diff_6731(iii,jjj) = 100d*abs(( yptsi_in_temp(0) -  reform(yptsi_in_APO_chianti8(iii,jjj,5)) ) /yptsi_in_temp(0))
 ;  diff_6716(iii,jjj) = 100d*abs(( yptsi_in_temp(1) - reform(yptsi_in_APO_chianti8(iii,jjj,6))  ) /yptsi_in_temp(1))
    
  endfor
  endfor
  
  p1= jade_spectrogram(   diff_6731,nel,nel_halfwidth,tec,tec_halfwidth,logdata=1,yrange=[0.1,100]);,ylog=1,xlog=1,
  
  stop
  
 ; p1= jade_spectrogram(   diff_6716,nel,nel_halfwidth,tec,tec_halfwidth,logdata=1,yrange=[0.1,100])

  
  
stop

end