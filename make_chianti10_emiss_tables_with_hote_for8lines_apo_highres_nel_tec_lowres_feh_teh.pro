
function calculate_IPT_emiss_citep2_output_yname_emiss_calc_given_feh_teh, Tec, nel,feh,Teh, Nsp, Ns2p, Ns3p,Ns4p, Nop, No2p, $
  min = min, max = max, xwavi=xwavi, yname = yname
  ;Tel in eV, nel in #/cm^3, Nsp etc in #/cm^2, max/min/fwhm in angstroms
  ; sample spectra calculate_IPT_emiss_citep2_output_yname
  if keyword_set(max) then maxwav = max else maxwav = 1800d
  if keyword_set(min) then minwav = min else minwav = 550d

  ;p_ratio=replicate(0.1d,n_elements(tel))
  conversion=11604.5221d ; conversion between eV and Kelvin from NIST
  log10tec=alog10(Tec*conversion)
  log10teh=alog10(Teh*conversion)
  log10nel=alog10(nel)
  
  fec = 1.d - feh

  ;tic
  ; emission line calculation: Sulfur ions
  ;s2em = emiss_calc(16, 2, temp = log10tel, dens = log10nel, /no_de, radt = 1.d, /quiet,/NOPROT,/NOIONREC,/NO_RREC );,proton_ratio=0.1d)
  ;s3em = emiss_calc(16, 3, temp = log10tel, dens = log10nel, /no_de, radt = 1.d, /quiet,/NOPROT,/NOIONREC,/NO_RREC );,proton_ratio=0.1d)
  ;s4em = emiss_calc(16, 4, temp = log10tel, dens = log10nel, /no_de, radt = 1.d, /quiet,/NOPROT,/NOIONREC,/NO_RREC );,proton_ratio=0.1d)
  ;s5em = emiss_calc(16, 5, temp = log10tel, dens = log10nel, /no_de, radt = 1.d, /quiet,/NOPROT,/NOIONREC,/NO_RREC );,proton_ratio=0.1d)


  ; emission line calculation: oxygen ions
  ;o2em = emiss_calc(8, 2, temp = log10tel, dens = log10nel, /no_de, radt = 1.d, /quiet,/NOPROT,/NOIONREC,/NO_RREC );,proton_ratio=0.1d)
  ;o3em = emiss_calc(8, 3, temp = log10tel, dens = log10nel, /no_de, radt = 1.d, /quiet,/NOPROT,/NOIONREC,/NO_RREC );,proton_ratio=0.1d)
  
  
  ;s1em = emiss_calc(16, 1,  temp = [log10tec,log10teh],sum_mwl_coeff=[fec,feh], dens = log10nel, /no_de, radt = 1.d, /quiet,/NOPROT,/NOIONREC,/NO_RREC)
  s2em = emiss_calc(16, 2,  temp = [ log10tec,log10teh],sum_mwl_coeff=[fec,feh], dens = log10nel, /no_de, radt = 1.d, /quiet,/NOPROT,/NOIONREC,/NO_RREC)
  s3em = emiss_calc(16, 3,  temp = [log10tec,log10teh],sum_mwl_coeff=[fec,feh], dens = log10nel, /no_de, radt = 1.d, /quiet,/NOPROT,/NOIONREC,/NO_RREC)
  s4em = emiss_calc(16, 4, temp = [log10tec,log10teh],sum_mwl_coeff=[fec,feh], dens = log10nel,  /no_de, radt = 1.d, /quiet,/NOPROT,/NOIONREC,/NO_RREC)
  s5em = emiss_calc(16,5, temp = [log10tec,log10teh],sum_mwl_coeff=[fec,feh], dens = log10nel, /no_de, radt = 1.d, /quiet,/NOPROT,/NOIONREC,/NO_RREC)
  ; emission line calculation: oxygen ions
  ;o1em = emiss_calc(8, 1,  temp = [log10tec,log10teh],sum_mwl_coeff=[fec,feh], dens = log10nel,/no_de, radt = 1.d, /quiet,/NOPROT,/NOIONREC,/NO_RREC)
  o2em = emiss_calc(8, 2,  temp = [log10tec,log10teh],sum_mwl_coeff=[fec,feh], dens = log10nel, /no_de, radt = 1.d, /quiet,/NOPROT,/NOIONREC,/NO_RREC)
  o3em = emiss_calc(8, 3,  temp = [log10tec,log10teh],sum_mwl_coeff=[fec,feh], dens = log10nel, /no_de, radt = 1.d, /quiet,/NOPROT,/NOIONREC,/NO_RREC)

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

function calculate_IPT_emiss_citep2_output_yname_emiss_calc_given_singlemax, Tel, nel, Nsp, Ns2p, Ns3p,Ns4p, Nop, No2p, $
  min = min, max = max, xwavi=xwavi, yname = yname
  ;Tel in eV, nel in #/cm^3, Nsp etc in #/cm^2, max/min/fwhm in angstroms
  ; sample spectra calculate_IPT_emiss_citep2_output_yname
  if keyword_set(max) then maxwav = max else maxwav = 1800d
  if keyword_set(min) then minwav = min else minwav = 550d

  ;p_ratio=replicate(0.1d,n_elements(tel))
  conversion=11604.5221d ; conversion between eV and Kelvin from NIST

  log10tel=alog10(Tel*conversion)
  log10nel=alog10(nel)



  ;tic
  ; emission line calculation: Sulfur ions
  s2em = emiss_calc(16, 2, temp = log10tel, dens = log10nel, /no_de, radt = 1.d, /quiet,/NOPROT,/NOIONREC,/NO_RREC );,proton_ratio=0.1d)
  s3em = emiss_calc(16, 3, temp = log10tel, dens = log10nel, /no_de, radt = 1.d, /quiet,/NOPROT,/NOIONREC,/NO_RREC );,proton_ratio=0.1d)
  s4em = emiss_calc(16, 4, temp = log10tel, dens = log10nel, /no_de, radt = 1.d, /quiet,/NOPROT,/NOIONREC,/NO_RREC );,proton_ratio=0.1d)
  s5em = emiss_calc(16, 5, temp = log10tel, dens = log10nel, /no_de, radt = 1.d, /quiet,/NOPROT,/NOIONREC,/NO_RREC );,proton_ratio=0.1d)


  ; emission line calculation: oxygen ions
  o2em = emiss_calc(8, 2, temp = log10tel, dens = log10nel, /no_de, radt = 1.d, /quiet,/NOPROT,/NOIONREC,/NO_RREC );,proton_ratio=0.1d)
  o3em = emiss_calc(8, 3, temp = log10tel, dens = log10nel, /no_de, radt = 1.d, /quiet,/NOPROT,/NOIONREC,/NO_RREC );,proton_ratio=0.1d)


  
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


Pro make_CHIANTI10_emiss_tables_with_hote_for8lines_apo_highres_nel_tec_lowres_feh_teh

tic
 

 ;stop
;xmin=550d;3000d;
;xmax=2100d;7000d;

xmin = 3000d;
xmax = 7000d;


nel = [dindgen(9) + 1d, 10d*dindgen(9) + 10d, 100d*dindgen(4) + 100d, 250d*dindgen(39) + 500d]

Tec = [0.1d*dindgen(9) + 0.1d, 0.5d*dindgen(17)+1d,dindgen(10)+10d,5d*dindgen(8) + 20d, 20d*dindgen(7) + 60d, 80d*dindgen(6) + 200d]



feh = [0.0025,0.01];[0.0001d,0.001d,0.0025d,0.005d,0.01d,0.05d,0.1d]
Teh = [0.0025,0.01];[5d*dindgen(4) + 35d, 20d*dindgen(3) + 60d, 100d*dindgen(3) + 200d]

n_ne = n_elements(nel)

n_te = n_elements(tec)
n_feh = n_elements(feh)
n_teh = n_elements(teh)

print,n_ne,n_te,n_feh,n_teh

;nxwavi = 884 ; for 550 - 2100 CHIANTI 8.0.1
;nxwavi = 4153 ; for 550 - 2100 CHIANTI 10.1

nxwavi = 2503 ; for 3000 - 7000 CHIANTI 10.1

yptsi_vary_nel_tec_feh_teh = dblarr(n_ne, n_te, n_feh, n_teh, nxwavi)
for i=0,n_ne -1 do begin
  print,'i = ',i
for j=0, n_te -1 do begin
  ;print,'j = ',j
  for k=0, n_feh -1 do begin  
;print,'k = ',k

  for l=0, n_teh -1 do begin 
  
  neli = nel(i)
  tecj = tec(j)
  
  fehk = feh(k)
  tehl = teh(l)
  

  
    yptsi=calculate_IPT_emiss_citep2_output_yname_emiss_calc_given_feh_teh( Tecj, neli,fehk,tehl, 1d, 1d, 1d,1d, 1d, 1d, $
    min = xmin, max = xmax, xwavi=xwavi, yname = yname)
    
    
    
    
    

   ; yptsi=calculate_IPT_emiss_citep2_output_yname_emiss_calc_given_singlemax( Tecj, neli, 1d, 1d, 1d,1d, 1d, 1d, $
   ;   min = xmin, max = xmax, xwavi=xwavi, yname = yname)

 
    
  
   ;print,n_elements(  yptsi)
   ; stop
   yptsi_vary_nel_tec_feh_teh(i,j,k,l,*) = yptsi;reform(yptsi,1,1,n_elements(  yptsi))

   
   
  ; print,xwavi
  ; print,yptsi
  ; print,n_elements(xwavi)
 ;stop

Endfor
endfor
Endfor
endfor
;,/NOIONREC,/NO_RREC
write_csv,'yptsi_3000-7000_Angstroms_vary_nec_tec_feh_teh_lookuptable_61x57x9x21x2503_CHIANTI10_1_noprotons_NOIONREC_NO_RREC.dat',yptsi_vary_nel_tec
;write_csv,'xwavi_3000-7000_Angstroms_CHIANTI10_1_noprotons_NOIONREC_NO_RREC.dat',xwavi
;write_csv,'yname_3000-7000_Angstroms_CHIANTI10_1_noprotons_NOIONREC_NO_RREC.dat',yname



openr,1,'yname_3000-7000_Angstroms_CHIANTI10_1_noprotons_NOIONREC_NO_RREC.dat'
yname_in=strarr(2503)
readf,1,yname_in
close,1

 yptsi_in = yptsi_vary_nel_tec_feh_teh
 
  xwavi_in = xwavi
  
   
   sp_idxs = where(yname_in eq '"S II"')
   s2p_idxs = where(yname_in eq '"S III"')
   s3p_idxs = where(yname_in eq '"S IV"')
   s4p_idxs = where(yname_in eq '"S V"')
   op_idxs = where(yname_in eq '"O II"')
   o2p_idxs = where(yname_in eq '"O III"')


APO_emiss_table = dblarr(n_ne, n_te, n_feh, n_teh, 8)


wavemax_want = 3730d;1732.92d;

wavemin_want = 3710d ;1707.97d;


species_want = s2p_idxs

idxxx = where((xwavi_in(species_want)  le wavemax_want) and (xwavi_in(species_want)  ge wavemin_want) )
xwavi_inn = xwavi_in(species_want)

yptsi_in_temp = reform(yptsi_in(28,17,0,17,species_want))
print,xwavi_inn(idxxx)
print, yptsi_in_temp(idxxx)/max(yptsi_in_temp(idxxx))

yptsi_in_temp = yptsi_in(*,*,*,*,species_want)
yptsi_in_temp =   yptsi_in_temp(*,*,*,*,idxxx)

APO_emiss_table(*,*,*,*,0) =   yptsi_in_temp(*,*,*,*,3) ; S2+ 3722 ANGSTROMS



;3rd index in final index
print,species_want(idxxx(3))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

wavemax_want = 3740d;1732.92d;

wavemin_want = 3720d ;1707.97d;


species_want = op_idxs

idxxx = where((xwavi_in(species_want)  le wavemax_want) and (xwavi_in(species_want)  ge wavemin_want) )
xwavi_inn = xwavi_in(species_want)

yptsi_in_temp = reform(yptsi_in(28,17,0,17,species_want))
print,xwavi_inn(idxxx)
print, yptsi_in_temp(idxxx)/max(yptsi_in_temp(idxxx))

yptsi_in_temp = yptsi_in(*,*,*,*,species_want)
yptsi_in_temp =   yptsi_in_temp(*,*,*,*,idxxx)

APO_emiss_table(*,*,*,*,1) =   yptsi_in_temp(*,*,*,*,0) ; O+ 3726 ANGSTROMS
APO_emiss_table(*,*,*,*,2) =   yptsi_in_temp(*,*,*,*,3) ; O+ 3729 ANGSTROMS

print,species_want(idxxx(0))
print,species_want(idxxx(3))



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

wavemax_want = 4085d;1732.92d;

wavemin_want = 4060d ;1707.97d;


species_want = sp_idxs

idxxx = where((xwavi_in(species_want)  le wavemax_want) and (xwavi_in(species_want)  ge wavemin_want) )
xwavi_inn = xwavi_in(species_want)

yptsi_in_temp = reform(yptsi_in(28,17,0,17,species_want))
print,xwavi_inn(idxxx)
print, yptsi_in_temp(idxxx)/max(yptsi_in_temp(idxxx))


yptsi_in_temp = yptsi_in(*,*,*,*,species_want)
yptsi_in_temp =   yptsi_in_temp(*,*,*,*,idxxx)

APO_emiss_table(*,*,*,*,3) =   yptsi_in_temp(*,*,*,*,1) ; S+ 4069 ANGSTROMS
APO_emiss_table(*,*,*,*,4) =   yptsi_in_temp(*,*,*,*,3) ; S+ 4076 ANGSTROMS


print,species_want(idxxx(1))
print,species_want(idxxx(3))




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

wavemax_want = 6320d;1732.92d;

wavemin_want = 6305d ;1707.97d;


species_want = s2p_idxs

idxxx = where((xwavi_in(species_want)  le wavemax_want) and (xwavi_in(species_want)  ge wavemin_want) )
xwavi_inn = xwavi_in(species_want)

yptsi_in_temp = reform(yptsi_in(28,17,0,17,species_want))
print,xwavi_inn(idxxx)
print, yptsi_in_temp(idxxx)/max(yptsi_in_temp(idxxx))

yptsi_in_temp = yptsi_in(*,*,*,*,species_want)

yptsi_in_temp =   yptsi_in_temp(*,*,*,*,idxxx)

APO_emiss_table(*,*,*,*,5) =   yptsi_in_temp(*,*,*,*,0) ; S++ 6312 ANGSTROMS

print,species_want(idxxx(0))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

wavemax_want = 6740d;1732.92d;

wavemin_want = 6705d ;1707.97d;


species_want = sp_idxs

idxxx = where((xwavi_in(species_want)  le wavemax_want) and (xwavi_in(species_want)  ge wavemin_want) )
xwavi_inn = xwavi_in(species_want)

yptsi_in_temp = reform(yptsi_in(28,17,0,17,species_want))
print,xwavi_inn(idxxx)
print, yptsi_in_temp(idxxx)/max(yptsi_in_temp(idxxx))

yptsi_in_temp = yptsi_in(*,*,*,*,species_want)

yptsi_in_temp =   yptsi_in_temp(*,*,*,*,idxxx)

APO_emiss_table(*,*,*,*,6) =   yptsi_in_temp(*,*,*,*,1) ; S+ 6716 ANGSTROMS
APO_emiss_table(*,*,*,*,7) =   yptsi_in_temp(*,*,*,*,2) ; S++ 6731 ANGSTROMS


print,species_want(idxxx(1))
print,species_want(idxxx(2))



write_csv,'yptsi_8_APO_LINES_optical_vary_nec_tec_feh_teh_lookuptable_61x57x9x57x8_CHIANTI_10.1.dat',APO_emiss_table


 
  
  
  
  toc
  stop


end




