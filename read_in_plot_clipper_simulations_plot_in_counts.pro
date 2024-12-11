pro read_in_plot_CLIPPER_simulations_plot_in_counts
n_intervals = 45

dt = 1d ; seconds
n_et_each_interval = floor(30d*60d/dt) + 1;1801 ;30 minutes for each interval at 1 second resolution is (30*60/0.01) + 1 = 1801, duration/dt + 1
;stop
xmin=550d;3000d;
xmax=2100d;7000d;
xstepsize=1d ; 0.6049d  ;4.47d
nx=Round((xmax-xmin)/xstepsize) + 1d
x_wave=dindgen(nx)*xstepsize + xmin
;print,x



;euvbin=0.6049d ;euv binsize cassini UVIS
spec_binsize = xstepsize;0.6049d ; cassini res ; 1d
fwhm = 6d; JUICE/CLIPPER UVS FWHM best case scenario (Angstrom)... ;4.47 ang was  UVIS res



; Read the CSV file into a 2D array
;data = read_csv('CLIPPER_UVS_Rayleighs_per_ang_num_et_each_interval=1801_num_intervals=45_num_xwavebins=1501-001.csv')

; Reshape the 2D array back to the original dimensions
;original_dimensions = [n_et_each_interval,n_intervals,n_elements(x_wave)];[1801, 45, 1551]
;ypts = reshape(data, original_dimensions)

restore,'clipper_simulation_ypts_out.sav',/verbose

restore, 'positions_and_densities_to_plot_for_3_et_intervals_for_each_of_45_Clipper_observations_v1.sav',/verbose
;save,s_outz,x_outz,y_outz,z_outz,utc_outz,tec_outz,nec_outz,nsp_outz,ns2p_outz,ns3p_outz,nop_outz,no2p_outz,filename='positions_and_densities_to_plot_for_3_et_intervals_for_each_of_45_Clipper_observations_v1.sav',/verbose


  ;openr,1,'CLIPPER_UVS_Rayleighs_per_ang_num_et_each_interval=1801_num_intervals=45_num_xwavebins=1501-001.csv'
  ;ypts = dblarr(n_et_each_interval,n_intervals,n_elements(x_wave)) ; ne x tec x discrete wavelength centers of each emission line
  ;readf,1,ypts
 ; close,1
 left = 0.2
 right = 0.2
 MARGINz0 =[left,0.1,right,0.2]
 MARGINz1 =[left,0.1,right,0.0]
 MARGINz2 =[left,0.25,right,0.0]

 restore,'Europa_UVS_EffArea_20220331.sav'

 
 p_out = OBJARR(n_intervals)
 lll=-1
 
 indicies_for_jjj_proxy = [0,900,1800] ; Beginning of et_interval, middle, and end. 
 
 strings_for_et_intervals = [', Beginning of 30 min Observation',', Middle of 30 min Observation',', End of 30 min Observation' ]
 strings_for_et_intervals_less = [', Beginning of 30 min',', Middle of 30 min',', End of 30 min' ]
n_tot_smooth = 5
for iii = 0, n_intervals -1 do begin ; 45 intervals 
;for jjj_proxy = 0 , 2 do begin ; only showing 9 frames per 30 minute interval 
 
 jjj = indicies_for_jjj_proxy[jjj_proxy]
 
 current_string =  strings_for_et_intervals[jjj_proxy]
 
 current_string_less = strings_for_et_intervals_less [jjj_proxy]
 ((3d/2.2d)^2d)*ypts_out(jjj,iii,*)


counts = 
  p1=plot(x_wave,counts,xtitle='Wavelength (Å)',YTITLE='Counts',title='Torus Observation ' + string(iii + 1) )
 
  lll = lll + 1
 p_out(lll) = p1
 

;endfor
endfor

for q=0, n_elements(p_out) -2 do p_out[q].Save, 'Clipper_UVS_plots_1s_resolution_only_3_plots_per_observation_with_densities_along_LOS.pdf', /APPEND

p_out[n_elements(p_out) -1].Save, 'Clipper_UVS_plots_1s_resolution_only_3_plots_per_observation_with_densities_along_LOS.pdf', /APPEND, /CLOSE


stop

end