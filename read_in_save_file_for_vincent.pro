pro read_in_save_file_for_vincent

  Restore,  FILENAME = 'CHIANTI_10.1_volume_emission_rates_per_Oxygen_density_for_Vincent.sav',/verbose
; EMISS_TABLE_OXYGEN   this is 61x57x3 element double array in Photon/s or volume emission rate/oxygen number density from CHIANTI 10.1 
; for electron impact excitation at 61x57 electron density and electron temeprature values
; Rayleighs in 1304 angstrom emission line = 10^-6 Integral of (volume emission rate/n_O)*n_O integrated over line of sight
;  given electron temp and density over line of sight assuming single maxwellian electron impact excitation of neutral O 
;  with no ionization/recombination or radiative recombination or proton impact excitation corrections
; NEC this is 61 element double array in cm^-3 from 
 ;TEC this is 57 element double array in eV
 ; WAVELENGTHS 3 element double array with wavelength labels for each 'feature' in Angstroms
;1304 Angstroms is actually made up of 3 discrete lines according to CHIANTI 10.1 at 1302.1680       1304.8580       1306.0291 Angstroms
;At 5 eV and 2000 cm^-3 these make up 56%, 33%, and 11% of total emission with a weighted wavelength center of 1303.4886 Angstroms
;1356 Angstroms is actually made up of 2 discrete lines according to CHIANTI 10.1 at 1355.5980       1358.5120 Angstroms
;At 5 eV and 2000 cm^-3 these make up 78%, and 22% of total emission with a weighted wavelength center of 1356.2449 Angstroms

;6302.05 is 6300 according to CHIANTI 10.1

;CHIANTI 10.1 dbase file for OI or neutral oxygen for info on cross sections/ocillator strengths used:
;%filename:  o_1.splups
;%oscillator strengths: Froese Fischer, C., & Tachiev, G. 2004, ADNDT, 87, 1
;%effective collision strengths: Zatsarinny, O., & Tayal, S.S. 2003, ApJS, 148, 575
;%effective collision strengths (levels 1-3): Bell, K.L., Berrington, K.A., & Thomas, M.R.J., 1998, MNRAS, 293, L83
;% produced as part of the Arcetri/Cambridge/NRL 'CHIANTI' atomic data base collaboration



p1=plot(tec[9:26],EMISS_TABLE_OXYGEN[28,9:26,0],xtitle='$T_{e}$ (eV)',ytitle='1304Å (Photons/s)',title='Volume Emission Rate/$n_O$ (Photons/s)',layout=[1,3,1],xlog=0,ylog=1,xrange=[1,10])
p2=plot(tec[9:26],EMISS_TABLE_OXYGEN[28,9:26,1],xtitle='$T_{e}$ (eV)',ytitle='1356Å (Photons/s)',layout=[1,3,2],/current,xlog=0,ylog=1,xrange=[1,10])
p3=plot(tec[9:26],EMISS_TABLE_OXYGEN[28,9:26,2],xtitle='$T_{e}$ (eV)',ytitle='6300Å (Photons/s)',layout=[1,3,3],/current,xlog=0,ylog=1,xrange=[1,10])


stop
end