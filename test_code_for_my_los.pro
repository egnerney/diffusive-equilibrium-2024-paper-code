FUNCTION lineOfSight, x0, y0, z0, xf, yf, zf, x, y, z, x_closest = x_closest, y_closest = y_closest, z_closest = z_closest
  ; Initialize variables
  dx = xf - x0
  dy = yf - y0
  dz = zf - z0
  xdir = SIgNum(dx)
  ydir = SIGNum(dy)
  zdir = SIGNum(dz)
  dx = ABS(dx)
  dy = ABS(dy)
  dz = ABS(dz)

  ; Line of sight in 3D space
  n = MAX([dx, dy, dz]) + 1

  ; Pre-allocate arrays for performance
  xline = DBLARR(n)
  yline = DBLARR(n)
  zline = DBLARR(n)

  ; Initiate Bresenham's line algorithm in 3D
  xi = 0d
 yi = xi
  zi = xi
  dxi = ROUND(dx/n)
  dyi = ROUND(dy/n)
  dzi = ROUND(dz/n)

  ; The main loop
  FOR i=0L, n-1 DO BEGIN
    xline[i] = x0 + xi * xdir
    yline[i] = y0 + yi * ydir
    zline[i] = z0 + zi * zdir
    IF xi < dx THEN xi += dxi
    IF yi < dy THEN yi += dyi
    IF zi < dz THEN zi += dzi
  ENDFOR

  ; Determine the closest points
  x_closest = AINTERPOLATE(x, xline, /INTERP)
  y_closest = AINTERPOLATE(y, yline, /INTERP)
  z_closest = AINTERPOLATE(z, zline, /INTERP)

  RETURN, x_closest
END

pro test_code_for_my_LOS
  x0 = 0d
  y0 = 0d
  z0 = 0d

  xf = 10d
  yf = 10d
  zf = 10d

  x = dINDGEN(20)/2.0d ; generates a grid with irregular spacing
  y = dINDGEN(20)/1.5d ; generates a grid with irregular spacing
  z = dINDGEN(20)/2.5d ; generates a grid with irregular spacing

  result = lineOfSight(x0, y0, z0, xf, yf, zf, x, y, z, x_closest = x_closest, y_closest = y_closest, z_closest = z_closest)

  PRINT, "x values:", x_closest
  PRINT, "y values:", y_closest
  PRINT, "z values:", z_closest


stop
end