pro read_in_ceq_radial_o4_APO_simulations_make_pdf_book_linear_interped
nframes= 90 
p_out = OBJARR(nframes)

openr,1,'linear_interpedhigher_res_simulated_B6731A_in_each_bin_o4_first_try_90framesx376_ceq_radialbins.csv'
B=fltarr(nframes,376)
readf,1,B
close,1

nframes = double(nframes)


dPhi = 360d/nframes

CMLdeg = dphi*findgen(nframes) + 270d

CMLdeg =  CMLdeg Mod 360d

xRight = 7.5d
xLeft = -7.5d
spacing = 0.04d

n_ceq_radial_bins = Floor((xRight - xLeft)/spacing) + 1;

rho_ceq = spacing*findgen(n_ceq_radial_bins) +  xLeft

p_mean = plot(rho_ceq,mean(B,dimension = 1),xtitle='$\rho_{ceq}$',ytitle='Rayleighs',title='Linear Interped each Frame, CML every 4$\deg$',name='Mean over CML',xrange=[4.5,7.0],yrange=[0,400])
p_median = plot(rho_ceq,median(B,dimension = 1),xtitle='$\rho_{ceq}$',ytitle='Rayleighs',title='Linear Interped each Frame, CML every 4$\deg$',name='Median over CML',/overplot,color='firebrick',xrange=[4.5,7.0],yrange=[0,400])
leg = LEGEND(TARGET=[p_mean ,p_median], POSITION=[7.5,600], $
  /DATA, /AUTO_TEXT_COLOR)
leg.save,'0to400_Dusk_no_offset_prediction_nominal_model_mean_and_median_interped_each_frame_then_foundo4_Ceq_pixels_6731Ang.png',resolution=300

p_avg_only_dusk = plot(rho_ceq,mean(B,dimension = 1),xtitle='$\rho_{ceq}$',ytitle='Rayleighs',title='Avg given CML every 4$\deg$',xrange=[4.5,7.5],name='Linear Interped then Avg over Longitude',yrange=[0,600])
p_avg_only_dusk_smooth3 = plot(rho_ceq,ts_smooth(mean(B,dimension = 1),3),xtitle='$\rho_{ceq}$',ytitle='Rayleighs',title='Avg given CML every 4$\deg$',xrange=[4.5,7.5],/overplot,color='red',name='3 point running avg of black',yrange=[0,600])
p_avg_only_dusk_smooth5 = plot(rho_ceq,ts_smooth(mean(B,dimension = 1),5),xtitle='$\rho_{ceq}$',ytitle='Rayleighs',title='Avg given CML every 4$\deg$',xrange=[4.5,7.5],/overplot,color='blue',name='5 point running avg of black',yrange=[0,600])

leg = LEGEND(TARGET=[p_avg_only_dusk ,p_avg_only_dusk_smooth3,p_avg_only_dusk_smooth5 ], POSITION=[7.5,600], $
  /DATA, /AUTO_TEXT_COLOR)

stop

for k= 0, nframes -1 do begin
  p_out(k) = plot(rho_ceq,reform(B(k,*)),xtitle='$\rho_{ceq}$',ytitle='Rayleighs',title='Linear Interp., LH Dusk $\lambda_{III}$ = ' + string(sigfig((360d - CMLdeg(k)) - 90d,4)) + '$\deg$',xrange=[4.5,7.7])
endfor

for q=0, nframes -2 do p_out[q].Save, 'LH_labeled_dusk_only_4.5-7.7_linear_interped_APO_6731_simulation_convolved_for_seeing_and_avged_over_each_bin_using_o4_constant_ceq.pdf', /APPEND

p_out[nframes -1].Save, 'LH_labeled_dusk_only_4.5-7.7_linear_interped_APO_6731_simulation_convolved_for_seeing_and_avged_over_each_bin_using_o4_constant_ceq.pdf', /APPEND, /CLOSE

stop

end