FUNCTION lineOfSight, x0, y0, z0, xf, yf, zf, x, y, z, x_closest = x_closest, y_closest = y_closest, z_closest = z_closest
  ; Initialize variables
  dx = xf - x0
  dy = yf - y0
  dz = zf - z0

  ; Line of sight in 3D space
  n = MAX([dx, dy, dz]) + 1

  ; Pre-allocate arrays for performance
  xline = DBLARR(n)
  yline = DBLARR(n)
  zline = DBLARR(n)

  ; Initiate variables
  xi = 0.0
  yi = 0.0
  zi = 0.0
  dxi = dx/n
  dyi = dy/n
  dzi = dz/n

  ; The main loop
  FOR i=0L, n-1 DO BEGIN
    xline[i] = x0 + xi
    yline[i] = y0 + yi
    zline[i] = z0 + zi
    IF xi < dx THEN xi += dxi
    IF yi < dy THEN yi += dyi
    IF zi < dz THEN zi += dzi
  ENDFOR

  ; Determine the closest points
  x_closest = DBLARR(n)
  y_closest = DBLARR(n)
  z_closest = DBLARR(n)

  FOR i=0L, n-1 DO BEGIN
    min_dist = !VALUES.F_INFINITY
    FOR j=0L, N_ELEMENTS(x)-1 DO BEGIN
      dist = SQRT((x[j] - xline[i])^2 + (y[j] - yline[i])^2 + (z[j] - zline[i])^2)
      IF dist < min_dist THEN BEGIN
        min_dist = dist
        x_closest[i] = x[j]
        y_closest[i] = y[j]
        z_closest[i] = z[j]
      ENDIF
    ENDFOR
  ENDFOR

  RETURN, x_closest
END

pro test_code_for_my_LOS_try3
  ; Example usage
  x0 = 1.0
  y0 = 2.0
  z0 = 3.0
  xf = 10.0
  yf = 12.0
  zf = 14.0

  x = [0.0, 2.0, 4.0, 6.0, 8.0, 10.0]
  y = [0.0, 3.0, 6.0, 9.0, 12.0]
  z = [0.0, 5.0, 10.0, 15.0]

  f=lineOfSight( x0, y0, z0, xf, yf, zf, x, y, z, x_closest = x_closest, y_closest = y_closest, z_closest = z_closest)



  PRINT, "x values:", x_closest
  PRINT, "y values:", y_closest
  PRINT, "z values:", z_closest

  stop
end