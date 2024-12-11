FUNCTION gaussian2D, x, z, sigma
  RETURN, EXP(-(x^2 + z^2) / (2 * sigma^2))
End


pro convolve_frames_then_sample_o4_eq_for_APO_like_carl

nframes = 90

  openr,1,'ypts_8apolinesa_vary_x0_and_z0_emission_table_90x279x129x8_movieframes_trydiff3_chianti10_duskfit1.2ions.txt'
  ypts_frames=dblarr(nframes,279,129,8)
  readf,1,ypts_frames
  close,1
  
 ypts_conv_out = dblarr(nframes, 279L*129L,8)
 
 


 dx1 = 0.1d
 nx1 = 40

 dx2 = 0.05d;0.025
 nx2 = 78;390 is for 0.01 ; 156 is 0.025 ; 78 for 0.05 ;
 ;429 for x and y for 0.025, 281 for 0.05
 nx3 = 42 ; middle bits at 0.1, was 44 before with mistake


 xgrid1d = [ dx1*dindgen(nx1) - 10.d , dx2*dindgen(nx2) - 6.d , dx1*dindgen(nx3) - 2.1d , dx2*dindgen(nx2) + 2.1d , dx1*dindgen(nx1+1) + 6.d ]

 dz1 = 0.1d
 nz1 = 12

 dz2 = 0.025d
 nz2 = 104;104 for 0.025 for -1.3 to 1.3 so 2.6RJ total

 zgrid1d = [dz1*dindgen(nz1) - 2.5d, dz2*dindgen(nz2) - 1.3d,dz1*dindgen(nz1 + 1) + 1.3d ];z_step_size*dindgen(n_z) + z_min

 z_half_step_sizegrid1d =  [ replicate(dz1/2.d,nz1) , replicate(dz2/2.d,nz2) , replicate(dz1/2.d,nz1 + 1) ]
 x_half_step_sizegrid1d =  [ replicate(dx1/2.d,nx1) , replicate(dx2/2.d,nx2) , replicate(dx1/2.d,nx3) , replicate(dx2/2.d,nx2) , replicate(dx1/2.d,nx1 + 1) ]
 ygrid1d = xgrid1d
 y_half_step_sizegrid1d  = x_half_step_sizegrid1d

 i_z0 = 64 ; z =0 index
 i_x0 = 139 ; x=0 index
 i_y0 = i_x0 ; y=0 index

 
 


  sigma = 0.065d / (2d * SQRT(2d * ALOG(2d)))

  ; Create 2D Gaussian kernel
  kernel2D = dblARR(279, 129)
  FOR i=0, 278 DO BEGIN
    print,i
    FOR j=0, 128 DO BEGIN
      kernel2D[i, j] = gaussian2D(xgrid1d[i], zgrid1d[j], sigma)
    ENDFOR
  ENDFOR

  kernel2D = kernel2D / TOTAL(kernel2D,/double)


  for jj_jj =0, 7 do begin
print,'convolution step emission line index = ', jj_jj


    ypts = reform(ypts_frames[*, *, *,jj_jj])


  ypts_conv = dblARR(nframes, 279, 129)
  FOR k=0, nframes - 1 DO BEGIN
    print,k
    ypts_conv[k, *, *] = CONVOL(reform(ypts[k, *, *]), kernel2D, /EDGE_TRUNCATE)
  ENDFOR

  data2D = REFORM(ypts_conv, nframes, 279L*129L)
  ;OPENW, lun, 'data2D.csv', /GET_LUN
  ;FOR i=0,nframes -1 DO PRINTF, lun, data2D[i,*]
  ;FREE_LUN, lun


ypts_conv_out(*, *,jj_jj) = data2D
     ;ypts_conv = data2D


endfor



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

 ; openr,1,'convolved_6731sp_90frames_for_movie_2D_mapped_from_3D_90x279x129.csv'
 ; ypts_conv=fltarr(nframes, 279L*129L)
  ;readf,1,ypts_conv
  ;close,1
  
  ;data2D = REFORM(ypts_conv, nframes, 279L*129L)
;ypts_conv = FLTARR(90, 279, 129)
;for jj_jj=0,7 do begin
  

ypts_conv = ypts_conv_out;(*, *,jj_jj)
ypts_conv = reform(ypts_conv,90, 279, 129,8) ; 90 frames, 279 rhovalues, 129 z values , 8 APO emission lines




;d\[Phi] = (360/nframes)
;CMLdeg =
;Flatten[Table[{Mod[270 + i*d\[Phi], 360]}, {i, 0, nframes - 1, 1}],
;1]
nframes = double(nframes);90d


dPhi = 360d/nframes

CMLdeg = dphi*dindgen(nframes) + 270d

xRight = 7.5d
xLeft = -7.5d
spacing = 0.04d

n_ceq_radial_bins = Floor((xRight - xLeft)/spacing) + 1;
print,n_ceq_radial_bins
bn = dblarr(n_ceq_radial_bins + 1)

B_sim = dblarr(nframes,n_ceq_radial_bins,8)


;B_frame = dblarr(n_elements(xgrid1d), n_elements(zgrid1d))
d = 0.025d

rho_z = dblarr(n_elements(xgrid1d), n_elements(zgrid1d),2)
c1cond = dblarr(n_elements(xgrid1d), n_elements(zgrid1d))
c2cond = dblarr(n_elements(xgrid1d), n_elements(zgrid1d))
c3cond = dblarr(n_elements(xgrid1d), n_elements(zgrid1d))
c4cond = dblarr(n_elements(xgrid1d), n_elements(zgrid1d))
;stop
rho1d_interpfrom = fltarr(279L,129L)
z1d_interpfrom = fltarr(279L,129L)


for i=0, 279L - 1 do begin
  for j=0, 129L - 1 do begin
    rho1d_interpfrom(i,j) = xgrid1d_old(i)
    z1d_interpfrom(i,j) = zgrid1d_old(j)
  endfor
endfor

rho1d_interpfrom = reform(rho1d_interpfrom,279L*129L)
z1d_interpfrom = reform(z1d_interpfrom,279L*129L)

for k = 0, nframes - 1 do begin

print,'k = ', k, ' of ', nframes - 1
B_frame = reform(ypts_conv(k,*,*,*),279L*129L,8)
;B_frame = reform(B_frame,279L*129L)
;triangulate,nel_,tec_,tr
;yptsi_LOS = griddata(  nel_, tec_,yptsi_,xout = nel_s, yout = tel_s,/linear,triangles = tr  )


triangulate,rho1d_interpfrom, z1d_interpfrom,tr
for jj_jj = 0, 7 do begin
  print,'interping APO emission line index=  ',jj_jj
B_frame(*,jj_jj) = griddata(rho1d_interpfrom, z1d_interpfrom,reform(B_frame(*,jj_jj),279L*129L),xout = rho_higher_res, yout = z_higher_res ,/grid,/linear,triangles = tr )
endfor
;print,'after1'

m = Tan((2d/3d)*ASin(-0.166769d*Sin((!dpi/180d)*(CMLdeg[k] + 21.7d))))

if (abs(m) > 0) then begin
  


perp_m = -1d/m 

b = Sqrt((d*Tan((2d/3d)*ASin(-0.166769d*Sin((!dpi/180d)*(CMLdeg[k] + 21.7d)))))^2d + d^2d)

b2 = xLeft*(m + 1d/m) - b;
b3 = xRight*(m + 1d/m) - b;



if (m > 0 ) then for n=0, n_ceq_radial_bins do bn(n) = b2 + n*spacing*sqrt((1d/m^2d) + 1d)
if (m < 0 ) then for n=0, n_ceq_radial_bins do bn(n) = b2 - n*spacing*sqrt((1d/m^2d) + 1d)


;line_above = m*rho + b
;line_below = m*rho - b
;perplines = -rho/m + bn

;p1=jade_spectrogram(B_frame,xgrid1d,x_half_step_sizegrid1d,zgrid1d,z_half_step_sizegrid1d)
for n = 0, n_ceq_radial_bins -1 do begin
  print,'n = ',n,' of 376'
  c1cond = dblarr(n_elements(xgrid1d), n_elements(zgrid1d))
  c2cond = dblarr(n_elements(xgrid1d), n_elements(zgrid1d))
  c3cond = dblarr(n_elements(xgrid1d), n_elements(zgrid1d))
  c4cond = dblarr(n_elements(xgrid1d), n_elements(zgrid1d))

for i=0, n_elements(xgrid1d) - 1 do begin
  for j=0, n_elements(zgrid1d) - 1 do begin
   ;rho_z(i,j,*) = [xgrid1d(i), zgrid1d(j)]
   c1cond(i,j) = xgrid1d(i) - (bn(n + 1) - zgrid1d(j))*m
   c2cond(i,j) = xgrid1d(i) - (bn(n) - zgrid1d(j))*m
   c3cond(i,j) = zgrid1d(j) - (m*xgrid1d(i) + b)
   c4cond(i,j) = zgrid1d(j) - (m*xgrid1d(i) - b)
    
  endfor
endfor

;B_frame = reform(B_frame,279L*129L)
C1cond = reform(C1cond,long64(n_elements(xgrid1d))*long64(n_elements(zgrid1d)))
C2cond = reform(C2cond,long64(n_elements(xgrid1d))*long64(n_elements(zgrid1d)))
C3cond = reform(C3cond,long64(n_elements(xgrid1d))*long64(n_elements(zgrid1d)))
C4cond = reform(C4cond,long64(n_elements(xgrid1d))*long64(n_elements(zgrid1d)))

idxs = where( (C1cond le 0) and (C2cond ge 0) and (C3cond le 0) and (C4cond ge 0)) ; idxs inside rectangular region of interest to find avg value of brightness within
for jj_jj =0, 7 do begin
  B_sim(k,n,jj_jj) = Total(B_frame(idxs,jj_jj),/double)/N_elements(B_frame(idxs,jj_jj))

endfor


endfor

endif else begin
  
  for n = 0, n_ceq_radial_bins -1 do begin

rho_z = dblarr(n_elements(xgrid1d), n_elements(zgrid1d),2)
    for i=0, n_elements(xgrid1d) - 1 do begin
      for j=0, n_elements(zgrid1d) - 1 do begin
       rho_z(i,j,*) = [xgrid1d(i), zgrid1d(j)]
        

      endfor
    endfor
    rho_z = reform(rho_z, long64(n_elements(xgrid1d))*long64(n_elements(zgrid1d)), 2 )
    ;B_frame = reform(B_frame,279L*129L)
  

    idxs = where( ( rho_z(*,0) le (xLeft + (n + 1d)*spacing)) and (rho_z(*,0) ge (xLeft + n*spacing)) and (rho_z(*,1) le d) and (rho_z(*,1) ge -d)) ; idxs inside rectangular region of interest to find avg value of brightness within

for jj_jj =0, 7 do begin
    B_sim(k,n,jj_jj) = Total(B_frame(idxs,jj_jj),/double)/N_elements(B_frame(idxs,jj_jj))
endfor

  endfor
  
  
  
endelse

endfor

;write_csv,'linear_interpedhigher_res_simulated_B6731A_in_each_bin_o4_first_try_90framesx376_ceq_radialbins.csv',B_sim
save,B_sim,filename='linear_interpedhigher_res_simulated_apo_8lines_dusk_90_frames_given_1.2xionsfit_90framesx376_ceq_radialbins_8lines.sav',/verbose



stop
end