
pro test_new_LOS_algorithm_vary_s_find_closest_grid_point

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
  
  
  norm_vec =[-sqrt(2d)/2d,sqrt(2d)/2d,0d] ;[0.d,1.d,0.d] ;looking in yhat direction
  x0 = 10.d*sqrt(2d)/2d
  y0 = -10.d*sqrt(2d)/2d
  z0 = 0.d
  
  slit_pos_vec=[x0,y0,z0]

s_i = [0.1d*dindgen(40),0.05d*dindgen(30) + 4d,0.1d*dindgen(90) + 5.5d ,0.05d*dindgen(30) + 14.5d ,0.1d*dindgen(41)+16d]
; , so maps to s between 0-4 0.1 RJ res, 4-5.5 0.05 RJ res, 5.5-14.5 0.1 RJ,
; ; 14.5 - 16 RJ 0.05RJ, 16-20 RJ 0.1RJ res


;stop
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


totaldx = xf - x0
totaldy = yf - y0
totaldz = zf - z0




flag = 1 ; proceed



; Define start and end points
start = [x0, y0, z0]
endd = [xf, yf, zf]



x_delta_dist = [ replicate(dx1,nx1) , replicate(dx2,nx2) , replicate(dx1,nx3) , replicate(dx2,nx2) , replicate(dx1,nx1+1) ]
y_delta_dist = x_delta_dist
z_delta_dist = [ replicate(dz1,nz1) , replicate(dz2,nz2) , replicate(dz1,nz1 +1) ]

k  = []
i_k = []

rho_gridsquared = x_los_i_in_torus*x_los_i_in_torus + y_los_i_in_torus*y_los_i_in_torus
;help,rho_gridsquared
;stop
tic
for i=0, n_elements(rho_gridsquared)-1 do begin

  rhogrid = rho_gridsquared(i)


  IF (rhogrid ge 20.25) AND (rhogrid le 100.) and (abs(z_los_i_in_torus(i)) le 2.5) THEN BEGIN ; 4.5^2 = 20.25 ; 4.75^2 = 22.5625
     print,i
    diff_r = ((x_los_i_in_torus(i) - xgrid_3d)*(x_los_i_in_torus(i) - xgrid_3d)) + ((y_los_i_in_torus(i) - ygrid_3d)*(y_los_i_in_torus(i) - ygrid_3d)) + ((z_los_i_in_torus(i) - zgrid_3d )*(z_los_i_in_torus(i) - zgrid_3d ))
    idx_map = where(diff_r eq min(diff_r))
    idx_map = idx_map(0)
    k = [k,idx_map(0)]
    i_k = [i_k,i]
  endif
endfor

k=k[UNIQ(k, SORT(k))]
help,k
print,xgrid_3d(k)
print,ygrid_3d(k)
print,zgrid_3d(k)
toc
;n=replicate({densities_LOS, s: 0.0d, sp: 0.0d, s2p: 0.0d, s3p: 0.0d, s4p:0.0d, $
;  o: 0.0d, op: 0.0d, o2p: 0.0d, el: 0.0d,elh:0.0d},n_elements(k))



;n.sp = n_in[k].sp;griddata(  rho_in, z_in,n_in.sp,xout = rho_los, yout = z_los )
;n.s2p = n_in[k].s2p; griddata( rho_in, z_in,n_in.s2p, xout = rho_los, yout = z_los )
;n.s3p = n_in[k].s3p  ;griddata(rho_in, z_in,n_in.s3p,  xout = rho_los, yout = z_los )
;n.op = n_in[k].op ;griddata( rho_in, z_in,n_in.op,  xout = rho_los, yout = z_los )
;n.o2p = n_in[k].o2p ; griddata( rho_in, z_in,n_in.o2p, xout = rho_los, yout = z_los )
;n.el = n_in[k].el  ;griddata(rho_in, z_in, n_in.el,  xout = rho_los, yout = z_los )
;tel =  T_in[k].el;griddata(rho_in, z_in, T_in.el,  xout = rho_los, yout = z_los )



endif else begin
  print,'nothing in LOS goes through torus so emission is 0'
  flag = 0 ; skip next part
endelse

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

stop


end