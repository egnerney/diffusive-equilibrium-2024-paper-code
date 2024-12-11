FUNCTION gaussian2D, x, z, sigma
  RETURN, EXP(-(x^2 + z^2) / (2 * sigma^2))
End

pro read_in_plot_frames

  


  openr,1,'ypts_vary_x0_and_z0_emission_table_45x279x129_new_way_UV_only680s2plines.txt'
  ypts=dblarr(45,279,129)
  readf,1,ypts
  close,1
  
  
  


  dx1 = 0.1
  nx1 = 40

  dx2 = 0.05;0.025
  nx2 = 78;390 is for 0.01 ; 156 is 0.025 ; 78 for 0.05 ;
  ;429 for x and y for 0.025, 281 for 0.05
  nx3 = 42 ; middle bits at 0.1, was 44 before with mistake


  xgrid1d = [ dx1*findgen(nx1) - 10. , dx2*findgen(nx2) - 6. , dx1*findgen(nx3) - 2.1 , dx2*findgen(nx2) + 2.1 , dx1*findgen(nx1+1) + 6. ]

  dz1 = 0.1
  nz1 = 12

  dz2 = 0.025
  nz2 = 104;104 for 0.025 for -1.3 to 1.3 so 2.6RJ total

  zgrid1d = [dz1*findgen(nz1) - 2.5, dz2*findgen(nz2) - 1.3,dz1*findgen(nz1 + 1) + 1.3 ];z_step_size*findgen(n_z) + z_min

  z_half_step_sizegrid1d =  [ replicate(dz1/2.,nz1) , replicate(dz2/2.,nz2) , replicate(dz1/2.,nz1 +1) ]
  x_half_step_sizegrid1d =  [ replicate(dx1/2.,nx1) , replicate(dx2/2.,nx2) , replicate(dx1/2.,nx3) , replicate(dx2/2.,nx2) , replicate(dx1/2.,nx1+1) ]
  ygrid1d = xgrid1d
  y_half_step_sizegrid1d  = x_half_step_sizegrid1d

  i_z0 = 64 ; z =0 index
  i_x0 = 139 ; x=0 index
  i_y0 = i_x0 ; y=0 index
  
  k=0
  
  p1=jade_spectrogram(reform(ypts(k,*,*)),xgrid1d,x_half_step_sizegrid1d,zgrid1d,z_half_step_sizegrid1d)

  openr,1,'ypts_8apolinesa_vary_x0_and_z0_emission_table_90x279x129x8_movieframes_trydiff3.txt'
  ypts=dblarr(90,279,129,8)
  readf,1,ypts
  close,1

ypts = reform(ypts[*, *, *,7])
nframes = 90

sigma = 0.065 / (2 * SQRT(2 * ALOG(2)))

; Create 2D Gaussian kernel
kernel2D = FLTARR(279, 129)
FOR i=0, 278 DO BEGIN
  print,i
    FOR j=0, 128 DO BEGIN
        kernel2D[i, j] = gaussian2D(xgrid1d[i], zgrid1d[j], sigma)
    ENDFOR
ENDFOR

kernel2D = kernel2D / TOTAL(kernel2D)



ypts_conv = FLTARR(90, 279, 129)
FOR k=0, 89 DO BEGIN
  print,k
    ypts_conv[k, *, *] = CONVOL(reform(ypts[k, *, *]), kernel2D, /EDGE_TRUNCATE)
ENDFOR

  data2D = REFORM(ypts_conv, nframes, 279L*129L)
  ;OPENW, lun, 'data2D.csv', /GET_LUN
  ;FOR i=0,nframes -1 DO PRINTF, lun, data2D[i,*]
  ;FREE_LUN, lun

write_csv,'convolved_6731sp_90frames_for_movie_2D_mapped_from_3D_90x279x129.csv', data2D


;write_csv,'xgrid_for_movie.csv',xgrid1d
;write_csv,'zgrid_for_movie.csv',zgrid1d

p1=jade_spectrogram(reform(ypts_conv(0,*,*)) - reform(ypts(0,*,*)),xgrid1d,x_half_step_sizegrid1d,zgrid1d,z_half_step_sizegrid1d)



stop
end