pro read_in_convolved_frames_sample_o4_eq_for_APO_like_carl

nframes= 90 
  openr,1,'convolved_6731sp_90frames_for_movie_2D_mapped_from_3D_90x279x129.csv'
  ypts_conv=fltarr(nframes, 279L*129L)
  readf,1,ypts_conv
  close,1
  
  ;data2D = REFORM(ypts_conv, nframes, 279L*129L)
;ypts_conv = FLTARR(90, 279, 129)
ypts_conv = reform(ypts_conv,90, 279, 129) ; 90 frames, 279 rhovalues, 129 z values 


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

;d\[Phi] = (360/nframes)
;CMLdeg =
;Flatten[Table[{Mod[270 + i*d\[Phi], 360]}, {i, 0, nframes - 1, 1}],
;1]
nframes = 90d


dPhi = 360d/nframes

CMLdeg = dphi*dindgen(nframes) + 270d

xRight = 7.5d
xLeft = -7.5d
spacing = 0.04d

n_ceq_radial_bins = Floor((xRight - xLeft)/spacing) + 1;
print,n_ceq_radial_bins
bn = dblarr(n_ceq_radial_bins + 1)

B_sim = dblarr(nframes,n_ceq_radial_bins)
;B_frame = dblarr(n_elements(xgrid1d), n_elements(zgrid1d))
d = 0.025d

rho_z = dblarr(n_elements(xgrid1d), n_elements(zgrid1d),2)
c1cond = dblarr(n_elements(xgrid1d), n_elements(zgrid1d))
c2cond = dblarr(n_elements(xgrid1d), n_elements(zgrid1d))
c3cond = dblarr(n_elements(xgrid1d), n_elements(zgrid1d))
c4cond = dblarr(n_elements(xgrid1d), n_elements(zgrid1d))
stop
for k = 0, nframes - 1 do begin

print,'k = ', k, ' of 89 '
B_frame = reform(ypts_conv(k,*,*))

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

p1=jade_spectrogram(B_frame,xgrid1d,x_half_step_sizegrid1d,zgrid1d,z_half_step_sizegrid1d)
for n = 0, n_ceq_radial_bins -1 do begin
  
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

B_frame = reform(B_frame,279L*129L)
C1cond = reform(C1cond,279L*129L)
C2cond = reform(C2cond,279L*129L)
C3cond = reform(C3cond,279L*129L)
C4cond = reform(C4cond,279L*129L)

idxs = where( (C1cond le 0) and (C2cond ge 0) and (C3cond le 0) and (C4cond ge 0)) ; idxs inside rectangular region of interest to find avg value of brightness within

B_sim(k,n) = Total(B_frame(idxs),/double)/N_elements(B_frame(idxs))


endfor

endif else begin
  
  for n = 0, n_ceq_radial_bins -1 do begin

rho_z = dblarr(n_elements(xgrid1d), n_elements(zgrid1d),2)
    for i=0, n_elements(xgrid1d) - 1 do begin
      for j=0, n_elements(zgrid1d) - 1 do begin
       rho_z(i,j,*) = [xgrid1d(i), zgrid1d(j)]
        

      endfor
    endfor
    rho_z = reform(rho_z, 279L*129L, 2 )
    B_frame = reform(B_frame,279L*129L)
    C1 = reform(C1,279L*129L)
    C2 = reform(C2,279L*129L)
    C3 = reform(C3,279L*129L)
    C4 = reform(C4,279L*129L)

    idxs = where( ( rho_z(*,0) le (xLeft + (n + 1d)*spacing)) and (rho_z(*,0) ge (xLeft + n*spacing)) and (rho_z(*,1) le d) and (rho_z(*,1) ge -d)) ; idxs inside rectangular region of interest to find avg value of brightness within

    B_sim(k,n) = Total(B_frame(idxs),/double)/N_elements(B_frame(idxs))


  endfor
  
  
  
endelse

endfor

write_csv,'simulated_B6731A_in_each_bin_o4_first_try_90framesx376_ceq_radialbins.csv',B_sim



stop
end