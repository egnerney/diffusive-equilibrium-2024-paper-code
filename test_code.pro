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
  closest_x = dblarr(num_steps+1)
  closest_y = dblarr(num_steps+1)
  closest_z = dblarr(num_steps+1)

  ; Initialize variables to track current position
  current_x = x0
  current_y = y0
  current_z = z0

  ; Perform ray tracing along the line of sight
  for i = 0, num_steps do begin
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
  indices = fltarr(num_steps+1)
  distance = fltarr(num_steps+1)
  for i = 0, num_steps do begin
    distance = sqrt((grid_x - closest_x[i])^2 + (grid_y - closest_y[i])^2 + (grid_z - closest_z[i])^2)
    indices[i] = where(distance eq min(distance), count)
  endfor

  ; Return the closest grid points
  closest_x = grid_x[indices]
  closest_y = grid_y[indices]
  closest_z = grid_z[indices]

  print, 'Closest grid points:'
  print, closest_x, closest_y, closest_z
  
  help,closest_x
  help,closest_y
  help,closest_z
  
  stop
end

pro test_code

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

stop
end