
pro find_closest_grid_points, x0, y0, z0, xf, yf, zf, grid_x, grid_y, grid_z
  ; Calculate direction vector of the line of sight
  dx = xf - x0
  dy = yf - y0
  dz = zf - z0

  ; Normalize the direction vector
  norm = sqrt(dx^2 + dy^2 + dz^2)
  dx = dx / norm
  dy = dy / norm
  dz = dz / norm

  ; Determine the maximum number of steps
  num_steps = max([abs(xf - x0), abs(yf - y0), abs(zf - z0)])

  ; Calculate step sizes along each axis
  step_x = (xf - x0) / num_steps
  step_y = (yf - y0) / num_steps
  step_z = (zf - z0) / num_steps

  ; Initialize variables to store closest grid points
  closest_x = fltarr(num_steps)
  closest_y = fltarr(num_steps)
  closest_z = fltarr(num_steps)

  ; Initialize variables to track current position
  current_x = x0
  current_y = y0
  current_z = z0

  ; Perform ray tracing along the line of sight
  for i = 0, num_steps-1 do begin
    ; Store the current grid point as the closest point
    closest_x[i] = current_x
    closest_y[i] = current_y
    closest_z[i] = current_z

    ; Move to the next grid point
    current_x = current_x + step_x
    current_y = current_y + step_y
    current_z = current_z + step_z
  endfor

  ; Find the grid points that are closest to the line of sight
  min_distance = 1e20  ; Initialize with a large value
  for i = 0, num_steps-1 do begin
    distance = sqrt((grid_x - closest_x[i])^2 + (grid_y - closest_y[i])^2 + (grid_z - closest_z[i])^2)
    indices = where(distance eq min(distance), count)
    if min(distance) lt min_distance then begin
      closest_x_final = closest_x[i]
      closest_y_final = closest_y[i]
      closest_z_final = closest_z[i]
      min_distance = min(distance)
    endif
  endfor

  ; Return the closest grid points
  print, 'Closest grid point:'
  print, closest_x_final, closest_y_final, closest_z_final
end

; Example usage
x0 = 1.0
y0 = 2.0
z0 = 3.0
xf = 10.0
yf = 12.0
zf = 14.0

grid_x = [0.0, 2.0, 4.0, 6.0, 8.0, 10.0]
grid_y = [0.0, 3.0, 6.0, 9.0, 12.0]
grid_z = [0.0, 5.0, 10.0, 15.0]

find_closest_grid_points, x0, y0, z0, xf, yf, zf, grid_x, grid_y, grid_z