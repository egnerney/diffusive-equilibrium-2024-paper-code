pro read_in_plot_CLIPPER_simulations
n_intervals = 45

dt = 1d ; seconds


n_et_each_interval = floor(30d*60d/dt) + 1;1801 ;30 minutes for each interval at 1 second resolution is (30*60/0.01) + 1 = 1801, duration/dt + 1
;stop
t = dt*dindgen(n_et_each_interval)
xmin=550d;3000d;
xmax=2100d;7000d;
xstepsize=1d ; 0.6049d  ;4.47d
nx=Round((xmax-xmin)/xstepsize) + 1d
x_wave=dindgen(nx)*xstepsize + xmin
;print,x
nspec_bins = nx
restore,'Europa_UVS_EffArea_20220331.sav'
effective_area_clipper = interpol(effarea_ap,w*10d,x_wave)

;help,x_wave
;help,w
print,w*10d
;w_diff = dblarr(n_elements(w) -1)
stop
;for i=0,n_elements(w) -2 do w_diff(i) = 10d*(w(i+1) - w(i))
;print,w_diff
;stop


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
;save,slit_pos_vec_outz,norm_vec_outz,s_outz,x_outz,y_outz,z_outz,utc_outz,tec_outz,nec_outz,nsp_outz,ns2p_outz,ns3p_outz,nop_outz,no2p_outz,filename='positions_and_densities_to_plot_for_3_et_intervals_for_each_of_45_Clipper_observations_v1.sav',/verbose


r_outz = Sqrt(x_outz^2d + y_outz^2d + z_outz^2d)

lat_outz = (180d/!dpi)*asin(z_outz/r_outz)

phi_outz = (180d/!dpi)*atan(y_outz,x_outz)

  ;openr,1,'CLIPPER_UVS_Rayleighs_per_ang_num_et_each_interval=1801_num_intervals=45_num_xwavebins=1501-001.csv'
  ;ypts = dblarr(n_et_each_interval,n_intervals,n_elements(x_wave)) ; ne x tec x discrete wavelength centers of each emission line
  ;readf,1,ypts
 ; close,1
 left = 0.2
 right = 0.2
 MARGINz0 =[left,0.1,right,0.2]
 MARGINz1 =[left,0.1,right,0.0]
 MARGINz2 =[left,0.25,right,0.0]

 
 p_out = OBJARR(n_intervals*3*3 + n_intervals)
 lll=-1
 
 indicies_for_jjj_proxy = [0,900,1800] ; Beginning of et_interval, middle, and end. 
 
 strings_for_et_intervals = [', Beginning of 30 min Observation',', Middle of 30 min Observation',', End of 30 min Observation' ]
 strings_for_et_intervals_less = [', Beginning of 30 min',', Middle of 30 min',', End of 30 min' ]
n_tot_smooth = 5
strds_per_pixel = 0.3d*0.1d*((!dpi/180d)^2d)
factor = ((10d^6d)/(4d*!dpi))*strds_per_pixel*effective_area_clipper

 wave_idx_want = where(x_wave le 2060d)
for iii = 0, n_intervals -1 do begin ; 45 intervals 
for jjj_proxy = 0 , 2 do begin ; only showing 9 frames per 30 minute interval 
  print,'iii = ', iii, ' of 44'
  print,'jjjproxy = ', jjj_proxy, ' of 2'
 jjj = indicies_for_jjj_proxy[jjj_proxy]
 
 current_string =  strings_for_et_intervals[jjj_proxy]
 
 current_string_less = strings_for_et_intervals_less [jjj_proxy]
 
 sample = utc_outz[iii,jjj_proxy,*] + ', $Pointing_{III}$ = [' + string(sigfig(norm_vec_outz[0,iii,jjj_proxy],3)) + ','+ string(sigfig(norm_vec_outz[1,iii,jjj_proxy],3)) + ',' + string(sigfig(norm_vec_outz[2,iii,jjj_proxy],3)) + ']' + ', $[x,y,z]_{III}$ = [' + string(sigfig(slit_pos_vec_outz[0,iii,jjj_proxy],3)) + ',' + string(sigfig(slit_pos_vec_outz[1,iii,jjj_proxy],3)) + ',' + string(sigfig(slit_pos_vec_outz[2,iii,jjj_proxy],3)) + ']'
 
 s_old = s_outz[iii,jjj_proxy,*]
;stop
 if ((finite(s_outz[iii,jjj_proxy,0]) eq 0) or (finite(s_outz[iii,jjj_proxy,1]) eq 0)) then begin
  sss = 0.025d*dindgen(floor(Sqrt(40d^2d + 40d^2d + 5d^2d)/0.025d) + 1)
  xxx = slit_pos_vec_outz[0,iii,jjj_proxy] + norm_vec_outz[0,iii,jjj_proxy]*sss
  yyy = slit_pos_vec_outz[1,iii,jjj_proxy] + norm_vec_outz[1,iii,jjj_proxy]*sss
  zzz = slit_pos_vec_outz[2,iii,jjj_proxy] + norm_vec_outz[2,iii,jjj_proxy]*sss
  s_outz[iii,jjj_proxy,*] = sss
  x_outz[iii,jjj_proxy,*] = xxx
  y_outz[iii,jjj_proxy,*] = yyy
  z_outz[iii,jjj_proxy,*] = zzz
  r_outz[iii,jjj_proxy,*] = sqrt(xxx^2d + yyy^2d + zzz^2d)
  lat_outz[iii,jjj_proxy,*] = (180d/!dpi)*asin(zzz/r_outz[iii,jjj_proxy,*])
  phi_outz[iii,jjj_proxy,*] = (180d/!dpi)*atan(yyy,xxx )
 endif else begin
  s_begin_have = s_outz[iii,jjj_proxy,0]
  n_need_s_more = s_begin_have/0.025d + 1
  s_more = 0.025d*dindgen(n_need_s_more)
  x_more = slit_pos_vec_outz[0,iii,jjj_proxy] + norm_vec_outz[0,iii,jjj_proxy]*s_more
  y_more = slit_pos_vec_outz[1,iii,jjj_proxy] + norm_vec_outz[1,iii,jjj_proxy]*s_more
  z_more = slit_pos_vec_outz[2,iii,jjj_proxy] + norm_vec_outz[2,iii,jjj_proxy]*s_more
  r_more = sqrt(x_more^2d + y_more^2d + z_more^2d)
  lat_more = (180d/!dpi)*asin(z_more/r_more)
  phi_more= (180d/!dpi)*atan(y_more,x_more )
  
  idx_want = where(finite(s_outz[iii,jjj_proxy,*]) eq 1)
  s_outz[iii,jjj_proxy,0:(n_elements(idx_want) + n_need_s_more - 1)] = [s_more,reform(s_outz[iii,jjj_proxy,idx_want])]
  x_outz[iii,jjj_proxy,0:(n_elements(idx_want) + n_need_s_more - 1)] = [x_more,reform(x_outz[iii,jjj_proxy,idx_want])]
  y_outz[iii,jjj_proxy,0:(n_elements(idx_want) + n_need_s_more - 1)] = [y_more,reform(y_outz[iii,jjj_proxy,idx_want])]
  z_outz[iii,jjj_proxy,0:(n_elements(idx_want) + n_need_s_more - 1)] = [z_more,reform(z_outz[iii,jjj_proxy,idx_want])]
  r_outz[iii,jjj_proxy,0:(n_elements(idx_want) + n_need_s_more - 1)] = [r_more,reform(r_outz[iii,jjj_proxy,idx_want])]
  lat_outz[iii,jjj_proxy,0:(n_elements(idx_want) + n_need_s_more - 1)] = [lat_more,reform(lat_outz[iii,jjj_proxy,idx_want])]
  phi_outz[iii,jjj_proxy,0:(n_elements(idx_want) + n_need_s_more - 1)] = [phi_more,reform(phi_outz[iii,jjj_proxy,idx_want])]

 endelse
;stop

   p0001=plot(s_outz[iii,jjj_proxy,*],r_outz[iii,jjj_proxy,*],layout=[2,3,1],ytitle='$r_{III}$ along LOS $(R_J)$',XTICKFORMAT="(A1)",title='Observation ' + string(iii + 1) +  current_string_less,margin=marginz0)
 p0002=plot(s_outz[iii,jjj_proxy,*],lat_outz[iii,jjj_proxy,*],layout=[2,3,3],/current,XTICKFORMAT="(A1)",margin=marginz1,ytitle='$Lat_{III}$ along LOS (Degrees)')
 p0003=plot(s_outz[iii,jjj_proxy,*],phi_outz[iii,jjj_proxy,*],layout=[2,3,5],/current,margin=marginz2,xtitle='s (LOS Distance, $R_J$)',ytitle='$\phi_{III}$ along LOS (Degrees)')
 p0001=plot(s_outz[iii,jjj_proxy,*],x_outz[iii,jjj_proxy,*],layout=[2,3,2],/current,ytitle='$x_{III}$ along LOS $(R_J)$',XTICKFORMAT="(A1)",margin=marginz0)
   p0002=plot(s_outz[iii,jjj_proxy,*],y_outz[iii,jjj_proxy,*],layout=[2,3,4],/current,XTICKFORMAT="(A1)",margin=marginz1,ytitle='$y_{III}$ along LOS $(R_J)$')
 p0003=plot(s_outz[iii,jjj_proxy,*],z_outz[iii,jjj_proxy,*],layout=[2,3,6],/current,margin=marginz2,xtitle='s (LOS Distance, $R_J$)',ytitle='$z_{III}$ along LOS $(R_J)$')

 t3 = TEXT( 0.02, 0.98, $
   sample,/NORMAL, /DEVICE)
    lll = lll + 1
 p_out(lll) = t3
 
 
 lll = lll + 1
 p01=plot( s_old,Ts_smooth(reform(tec_outz[iii,jjj_proxy,*]),n_tot_smooth),layout=[2,4,1],ytitle='$T_{ec}$ (eV)',XTICKFORMAT="(A1)",title='Observation ' + string(iii + 1) +  current_string_less,margin=marginz0)
 p02=plot( s_old,Ts_smooth(reform((3d/2.2d)*nec_outz[iii,jjj_proxy,*]),n_tot_smooth),layout=[2,4,2],/current,XTICKFORMAT="(A1)",margin=marginz0,ytitle='$n_{ec}$ ($cm^{-3}$)')
 p03=plot( s_old,Ts_smooth(reform((3d/2.2d)*nsp_outz[iii,jjj_proxy,*]),n_tot_smooth),layout=[2,4,3],/current,XTICKFORMAT="(A1)",margin=marginz1,ytitle='$n_{S^{+}}$ ($cm^{-3}$)')
 p04=plot( s_old,Ts_smooth(reform((3d/2.2d)*ns2p_outz[iii,jjj_proxy,*]),n_tot_smooth),layout=[2,4,4],/current,XTICKFORMAT="(A1)",margin=marginz1,ytitle='$n_{S^{++}}$ ($cm^{-3}$)')
 p05=plot( s_old,Ts_smooth(reform((3d/2.2d)*ns3p_outz[iii,jjj_proxy,*]),n_tot_smooth),layout=[2,4,5],/current,XTICKFORMAT="(A1)",margin=marginz1,ytitle='$n_{S^{+++}}$ ($cm^{-3}$)')
 p06=plot( s_old,Ts_smooth(reform((3d/2.2d)*nop_outz[iii,jjj_proxy,*]),n_tot_smooth),layout=[2,4,6],/current,xtitle='s (LOS Distance, $R_J$)',margin=marginz2,ytitle='$n_{O^{+}}$ ($cm^{-3}$)')
 p07=plot( s_old,Ts_smooth(reform((3d/2.2d)*no2p_outz[iii,jjj_proxy,*]),n_tot_smooth),layout=[2,4,7],/current,xtitle='s (LOS Distance, $R_J$)',margin=marginz2,ytitle='$n_{O^{++}}$ ($cm^{-3}$)')
 t03 = TEXT( 0.02, 0.98, $
 sample,/NORMAL, /DEVICE)
 ;stop
 p_out(lll) = t03
 

  p1=plot(x_wave(wave_idx_want)/10d,((3d/2.2d)^2d)*ypts_out(jjj,iii,wave_idx_want)*10d,xtitle='Wavelength (nm)',YTITLE='Rayleighs/nm',title='Torus Observation ' + string(iii + 1) +  current_string,xrange=[50,210] )
  t003 = TEXT( 0.02, 0.98, $
    sample,/NORMAL, /DEVICE)
  lll = lll + 1
 p_out(lll) =t003
 

endfor
counts = dblarr(nspec_bins)
for mmm=0, nspec_bins - 1 do begin
  counts(mmm) = int_tabulated(t,factor[mmm]*((3d/2.2d)^2d)*reform(ypts_out(*,iii,mmm)),/double)
endfor

sample2 = '30 Min Torus Stare Starting at UTC = ' + utc_outz[iii,0,*]


p001=plot(x_wave(wave_idx_want)/10d,counts(wave_idx_want),xtitle='Wavelength (nm)',YTITLE='Counts',title='Torus Observation ' + string(iii + 1) + ', 1 second resolution',/ylog,xrange=[50,210] )

t0003 = TEXT( 0.02, 0.98, $
  sample2,/NORMAL, /DEVICE)

lll = lll + 1
p_out(lll) = t0003
;stop
endfor

for q=0, n_elements(p_out) -2 do p_out[q].Save, 'Clipper_UVS_AP_torus_stare_plots_1s_resolution_10_plots_per_observation_with_densities_locations_along_LOS.pdf', /APPEND

p_out[n_elements(p_out) -1].Save, 'Clipper_UVS_AP_torus_stare_plots_1s_resolution_10_plots_per_observation_with_densities_locations_along_LOS.pdf', /APPEND, /CLOSE


stop

end