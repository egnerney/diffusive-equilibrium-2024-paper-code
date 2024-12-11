
function calculate_IPT_emiss_citep2_output_yname_emiss_calc_given_neutrals_onlyO, Tel, nel, No,  $
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
  ; emission line calculation: Sulfur neutral
  

  ; emission line calculation: oxygen neutral
  o1em = emiss_calc(8, 1, temp = log10tel, dens = log10nel, /no_de, radt = 1.d, /quiet,/NOPROT,/NOIONREC,/NO_RREC );,proton_ratio=0.1d)
  ;toc
  ;help,s3em.em
  ;stop


  ;puts in rayleighs, nsp etc. in local density units #/cm^3

  o1emiss = reform(o1em.em)*No

  ; discrete
  xwavi = [ o1em.lambda]
  wsort = sort(xwavi) ; sort by wavelength
  xwavi = xwavi[wsort]

  yptsi = [ o1emiss]
  yptsi = yptsi[wsort]

  yname = [ o1em.ion_name]
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


Pro make_CHIANTI10_emiss_tables_neutral_only6300_O

tic
 

 ;stop
;xmin=550d;3000d;
;xmax=2100d;7000d;

xmin = 6295d;
xmax = 6305d;


nel = [dindgen(9) + 1d, 10d*dindgen(9) + 10d, 100d*dindgen(4) + 100d, 250d*dindgen(39) + 500d]

Tec = [0.1d*dindgen(9) + 0.1d, 0.5d*dindgen(17)+1d,dindgen(10)+10d,5d*dindgen(8) + 20d, 20d*dindgen(7) + 60d, 80d*dindgen(6) + 200d]

n_ne = n_elements(nel)

n_te = n_elements(tec)

;nxwavi = 884 ; for 550 - 2100 CHIANTI 8.0.1
;nxwavi = 4153 ; for 550 - 2100 CHIANTI 10.1

 ; for 550 - 2100 CHIANTI 10.1

yptsi_vary_nel_tec = dblarr(n_ne, n_te)
for i=0,n_ne -1 do begin
  print,'i = ',i
for j=0, n_te -1 do begin  
  
  neli = nel(i)
  tecj = tec(j)
  
    yptsi=calculate_IPT_emiss_citep2_output_yname_emiss_calc_given_neutrals_onlyO( Tecj, neli, 1d, $
    min = xmin, max = xmax, xwavi=xwavi, yname = yname)
   
    ;print,n_elements(  yptsi)

   yptsi_vary_nel_tec(i,j) = yptsi;reform(yptsi,1,1,n_elements(  yptsi))
   print,xwavi
   print,yptsi
   print,n_elements(xwavi)
 stop

Endfor
endfor
;,/NOIONREC,/NO_RREC
write_csv,'neutral_onlyO_yptsi_6300_Angstroms_vary_ne_te_lookuptable_61x57_CHIANTI10_1_noprotons_NOIONREC_NO_RREC.dat',yptsi_vary_nel_tec

  toc
  stop


end




