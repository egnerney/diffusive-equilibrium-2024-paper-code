
pro test_new_LOS_algorithm

  ;stop
  dx1 = 0.1
  nx1 = 40

  dx2 = 0.05;0.025
  nx2 = 78;390 is for 0.01 ; 156 is 0.025 ; 78 for 0.05 ;
  ;429 for x and y for 0.025, 281 for 0.05
  nx3 = 42 ; middle bits at 0.1, was 44 with mistake from before


  xgrid1d = [ dx1*findgen(nx1) - 10. , dx2*findgen(nx2) - 6. , dx1*findgen(nx3) - 2.1 , dx2*findgen(nx2) + 2.1 , dx1*findgen(nx1+1) + 6. ]



  ;z_min = -2.5;
  ;z_max = 2.5;
  ;z_step_size = 0.05 ;  RJ
  ;n_z = Round((z_max - z_min)/z_step_size) + 1
  dz1 = 0.1
  nz1 = 12

  dz2 = 0.025
  nz2 = 104;104 for 0.025 ; 80 for 0.05; 260 for 0.01 for -1.3 to 1.3 so 2.6RJ total
  ;285 elements for zgrid if dz=0.01 for z= -1.3 to 1.3, 129 if 0.025
  zgrid1d = [dz1*findgen(nz1) - 2.5, dz2*findgen(nz2) - 1.3,dz1*findgen(nz1 + 1) + 1.3 ];z_step_size*findgen(n_z) + z_min

  z_half_step_sizegrid1d =  [ replicate(dz1/2.,nz1) , replicate(dz2/2.,nz2) , replicate(dz1/2.,nz1 +1) ]
  x_half_step_sizegrid1d =  [ replicate(dx1/2.,nx1) , replicate(dx2/2.,nx2) , replicate(dx1/2.,nx3) , replicate(dx2/2.,nx2) , replicate(dx1/2.,nx1+1) ]
  ygrid1d = xgrid1d
  y_half_step_sizegrid1d  = x_half_step_sizegrid1d

  i_z0 = 64
  i_x0 = 139
  i_y0 = i_x0
  
  
  xgrid = double(xgrid1d)
  ygrid = double(ygrid1d)
  zgrid = double(zgrid1d)
  
  nxgrid = long64(n_elements(xgrid))
  nygrid = long64(n_elements(ygrid))
  nzgrid = long64(n_elements(zgrid))
  nxgridm1 = nxgrid -1
  nygridm1 = nygrid -1
  nzgridm1 = nzgrid -1
  
  n_3dcartgrid  = nxgrid*nygrid*nzgrid

  n_3dcartgridm1 =  n_3dcartgrid - 1
  
  
  xgrid_3d = dblarr(nxgrid,nygrid,nzgrid)
  ygrid_3d = xgrid_3d
  zgrid_3d = xgrid_3d

  for i=0,nxgridm1  do begin
    for j=0, nygridm1 do begin
      for k=0,nzgridm1 do begin
        xgrid_3d(i,j,k) = xgrid(i)
        ygrid_3d(i,j,k) = ygrid(j)
        zgrid_3d(i,j,k) = zgrid(k)
      endfor
    endfor
  endfor
  
  
  norm_vec = [0.d,1.d,0.d] ;looking in yhat direction
  x0 = 0.d
  y0 = -10.d
  z0 = 0.d
  
  slit_pos_vec=[x0,y0,z0]

s_i = 0.1d*dindgen(211) ; 0, 21 in 0.1 RJ increments

x_los_i=slit_pos_vec[0] + s_i*norm_vec[0]; in RJ
y_los_i=slit_pos_vec[1] + s_i*norm_vec[1]; in RJ
z_los_i=slit_pos_vec[2] + s_i*norm_vec[2]; in RJ





x_idxs_in_torus = where((x_los_i ge -10d) and (x_los_i le 10d) )
y_idxs_in_torus = where((y_los_i ge -10d) and (y_los_i le 10d) )
z_idxs_in_torus = where((z_los_i ge -2.5d) and (z_los_i le 2.5d) )

if ((x_idxs_in_torus(0) ge 0) or (y_idxs_in_torus(0) ge 0) or (z_idxs_in_torus(0) ge 0)) then begin ; if less than 0 then no LOS through torus
  


x_los_i_in_torus = x_los_i(x_idxs_in_torus)
y_los_i_in_torus = y_los_i(y_idxs_in_torus)
z_los_i_in_torus = z_los_i(z_idxs_in_torus)

x0 = min(x_los_i_in_torus)
xf = max(x_los_i_in_torus)
y0 = min(y_los_i_in_torus)
yf = max(y_los_i_in_torus)
z0 = min(z_los_i_in_torus)
zf = max(z_los_i_in_torus)


dx = xf - x0
dy = yf - y0
dz = zf - z0




flag = 1 ; proceed



start = [x0, y0, z0]
endd = [xf, yf, zf]

; Define direction of the ray
ray_dir = (endd - start) / sqrt(total((endd - start)^2))

; Define the increment in each direction to move one cell along the ray
delta_dist = abs(1. / ray_dir)

; Determine the sign of the increment in each direction
step = signum(ray_dir)

; Determine how far along the ray we must move in each direction to cross a cell boundary
side_dist = delta_dist * (step * (ceil(step * start) - start))

; Define an array to hold the grid indices that the ray passes through
grid_indices = lonarr(3, n_3dcartgrid)

; Initialize a counter for the number of grid cells crossed
num_cells_crossed = 0

; Begin the loop over the grid cells
for i=0L, nxgridm1 do begin
  for j=0L, nygridm1 do begin
    for k=0L, nzgridm1 do begin
      ; Check if the cell at the current position is inside the region of interest
      if (start[0] ge xgrid_3d[i,j,k] && start[0] le xgrid_3d[i,j,k] + dx1 && $
        start[1] ge ygrid_3d[i,j,k] && start[1] le ygrid_3d[i,j,k] + dx1 && $
        start[2] ge zgrid_3d[i,j,k] && start[2] le zgrid_3d[i,j,k] + dz1) then begin
        ; Record the indices of the grid cell
        grid_indices[*, num_cells_crossed] = [i, j, k]
        num_cells_crossed += 1L
      endif

      ; Move to the next cell along the ray
      if side_dist[0] < side_dist[1] then begin
        if side_dist[0] < side_dist[2] then begin
          side_dist[0] += delta_dist[0]
          start[0] += step[0]
        endif else begin
          side_dist[2] += delta_dist[2]
          start[2] += step[2]
        endelse
      endif else begin
        if side_dist[1] < side_dist[2] then begin
          side_dist[1] += delta_dist[1]
          start[1] += step[1]
        endif else begin
          side_dist[2] += delta_dist[2]
          start[2] += step[2]
        endelse
      endelse
   
    endfor
  endfor
endfor

; Resize the grid_indices array to remove unused elements
grid_indices = reform(grid_indices, 3, num_cells_crossed)

endif else begin
  print,'nothing in LOS goes through torus so emission is 0'
  flag = 0 ; skip next part
endelse

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
help, grid_indices
stop


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;





stop
end