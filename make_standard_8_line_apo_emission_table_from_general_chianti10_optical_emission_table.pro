Pro make_standard_8_line_APO_emission_table_from_general_CHIANTI10_optical_emission_table


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



APO_emiss_table = dblarr(61,57,8)


wavemax_want = 3730d;1732.92d;

wavemin_want = 3710d ;1707.97d;


species_want = s2p_idxs

  idxxx = where((xwavi_in(species_want)  le wavemax_want) and (xwavi_in(species_want)  ge wavemin_want) )
 xwavi_inn = xwavi_in(species_want)

  yptsi_in_temp = reform(yptsi_in(28,17,species_want))
  print,xwavi_inn(idxxx)
  print, yptsi_in_temp(idxxx)/max(yptsi_in_temp(idxxx))
  
  yptsi_in_temp = yptsi_in(*,*,species_want)
    yptsi_in_temp =   yptsi_in_temp(*,*,idxxx)
  
  APO_emiss_table(*,*,0) =   yptsi_in_temp(*,*,3) ; S2+ 3722 ANGSTROMS
  
  
  
  ;3rd index in final index
  print,species_want(idxxx(3))


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  wavemax_want = 3740d;1732.92d;

  wavemin_want = 3720d ;1707.97d;

  
  species_want = op_idxs

  idxxx = where((xwavi_in(species_want)  le wavemax_want) and (xwavi_in(species_want)  ge wavemin_want) )
  xwavi_inn = xwavi_in(species_want)

  yptsi_in_temp = reform(yptsi_in(28,17,species_want))
  print,xwavi_inn(idxxx)
  print, yptsi_in_temp(idxxx)/max(yptsi_in_temp(idxxx))

  yptsi_in_temp = yptsi_in(*,*,species_want)
    yptsi_in_temp =   yptsi_in_temp(*,*,idxxx)

  APO_emiss_table(*,*,1) =   yptsi_in_temp(*,*,0) ; O+ 3726 ANGSTROMS
  APO_emiss_table(*,*,2) =   yptsi_in_temp(*,*,3) ; O+ 3729 ANGSTROMS
  
  print,species_want(idxxx(0))
   print,species_want(idxxx(3))

  
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  wavemax_want = 4085d;1732.92d;

  wavemin_want = 4060d ;1707.97d;


  species_want = sp_idxs

  idxxx = where((xwavi_in(species_want)  le wavemax_want) and (xwavi_in(species_want)  ge wavemin_want) )
  xwavi_inn = xwavi_in(species_want)

  yptsi_in_temp = reform(yptsi_in(28,17,species_want))
  print,xwavi_inn(idxxx)
  print, yptsi_in_temp(idxxx)/max(yptsi_in_temp(idxxx))


  yptsi_in_temp = yptsi_in(*,*,species_want)
  yptsi_in_temp =   yptsi_in_temp(*,*,idxxx)

  APO_emiss_table(*,*,3) =   yptsi_in_temp(*,*,1) ; S+ 4069 ANGSTROMS
  APO_emiss_table(*,*,4) =   yptsi_in_temp(*,*,3) ; S+ 4076 ANGSTROMS


  print,species_want(idxxx(1))
  print,species_want(idxxx(3))

 


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  wavemax_want = 6320d;1732.92d;

  wavemin_want = 6305d ;1707.97d;


  species_want = s2p_idxs

  idxxx = where((xwavi_in(species_want)  le wavemax_want) and (xwavi_in(species_want)  ge wavemin_want) )
  xwavi_inn = xwavi_in(species_want)

  yptsi_in_temp = reform(yptsi_in(28,17,species_want))
  print,xwavi_inn(idxxx)
  print, yptsi_in_temp(idxxx)/max(yptsi_in_temp(idxxx))

  yptsi_in_temp = yptsi_in(*,*,species_want)
  
    yptsi_in_temp =   yptsi_in_temp(*,*,idxxx)

  APO_emiss_table(*,*,5) =   yptsi_in_temp(*,*,0) ; S++ 6312 ANGSTROMS
  
  print,species_want(idxxx(0))
  
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  wavemax_want = 6740d;1732.92d;

  wavemin_want = 6705d ;1707.97d;


  species_want = sp_idxs

  idxxx = where((xwavi_in(species_want)  le wavemax_want) and (xwavi_in(species_want)  ge wavemin_want) )
  xwavi_inn = xwavi_in(species_want)

  yptsi_in_temp = reform(yptsi_in(28,17,species_want))
  print,xwavi_inn(idxxx)
  print, yptsi_in_temp(idxxx)/max(yptsi_in_temp(idxxx))

  yptsi_in_temp = yptsi_in(*,*,species_want)

  yptsi_in_temp =   yptsi_in_temp(*,*,idxxx)

  APO_emiss_table(*,*,6) =   yptsi_in_temp(*,*,1) ; S+ 6716 ANGSTROMS
  APO_emiss_table(*,*,7) =   yptsi_in_temp(*,*,2) ; S++ 6731 ANGSTROMS


  print,species_want(idxxx(1))
  print,species_want(idxxx(2))

  stop

write_csv,'yptsi_8_APO_LINES_optical_vary_ne_te_lookuptable_61x57x8_CHIANTI_10.1.dat',APO_emiss_table
  


  
  STOP
  
  

  

  
  
  
stop

end