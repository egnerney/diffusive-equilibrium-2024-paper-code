Pro make_oxygen_emission_rate_tables_ch10_noIONRECandnorrec_noprotons_for_vincent


  openr,1,'neutrals_yptsi_550-2100_Angstroms_vary_ne_te_lookuptable_61x57x7_CHIANTI10_1_noprotons_NOIONREC_NO_RREC.dat'
  yptsi_in=dblarr(61,57,7) ; ne x tec x discrete wavelength centers of each emission line
  readf,1,yptsi_in
  close,1
  
  openr,1,'neutral_onlyO_yptsi_6300_Angstroms_vary_ne_te_lookuptable_61x57_CHIANTI10_1_noprotons_NOIONREC_NO_RREC.dat'
  yptsi_in_6300=dblarr(61,57) ; ne x tec x discrete wavelength centers of each emission line
  readf,1,yptsi_in_6300
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
 
  nec = [dindgen(9) + 1d, 10d*dindgen(9) + 10d, 100d*dindgen(4) + 100d, 250d*dindgen(39) + 500d]

  tec = [0.1d*dindgen(9) + 0.1d, 0.5d*dindgen(17)+1d,dindgen(10)+10d,5d*dindgen(8) + 20d, 20d*dindgen(7) + 60d, 80d*dindgen(6) + 200d]

  n_ne = n_elements(nec)

  n_te = n_elements(tec)
  
  ;stop
  
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

  

emiss_table_oxygen = dblarr(n_ne,n_te,3) ;1304Å, 1356Å, 6300 ang with dimensions 61x57x3

wavelengths = [1304,1356,6300]
;stop
  
  for iii = 0, n_ne - 1 do begin
    print,iii
    for jjj = 0, n_te - 1 do begin
    
    
    yptsi_in_temp = reform(yptsi_in(iii,jjj,o_idxs))
    yptsi_in_temp3 = reform(yptsi_in_6300(iii,jjj))
    
    
 emiss_table_oxygen(iii,jjj,*) = [Total(yptsi_in_temp(idxxx_1304),/double),Total(yptsi_in_temp(idxxx_1356),/double),yptsi_in_temp3]
    
  endfor
  endfor
  
  SAVE, emiss_table_oxygen,nec,tec,wavelengths, FILENAME = 'CHIANTI_10.1_volume_emission_rates_per_Oxygen_density_for_Vincent.sav'

  
stop
end