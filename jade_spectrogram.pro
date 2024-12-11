function jade_spectrogram,im_aaaaaaa,x_aaaaaaa,dx_aaaaaaa,y_aaaaaaa,dy_aaaaaaa, $
  LOGDATA = LOGDATA,  XYPLOT = XYPLOT, $
  XRANGE = XRANGE, YRANGE = YRANGE, CLIM = CLIM, XTIME = XTIME, $
  XTITLE = XTITLE, YTITLE = YTITLE, CTITLE = CTITLE, TITLE = TITLE, $
  XLABEL = XLABEL, YLABEL = YLABEL, CLABEL = CLABEL  ; *LABEL is duplicating *TITLE for Rob's convenince and typos from Matlab
  
  ON_ERROR,2
  
  IF KEYWORD_SET(LOGDATA) THEN LogData = 1 ELSE LogData = 0
  IF KEYWORD_SET(XTIME) THEN XTIME = 1 ELSE XTIME = 0 ; is x a Julian time?
  IF KEYWORD_SET(XYPLOT) THEN XYPLOT = 1 ELSE XYPLOT = 0 ; line plot not spectorgram
  
  IF (N_ELEMENTS(x_aaaaaaa)+N_ELEMENTS(dx_aaaaaaa)+N_ELEMENTS(y_aaaaaaa)+N_ELEMENTS(dy_aaaaaaa) EQ 0) THEN BEGIN
    im_size = size(im_aaaaaaa)
    IF (im_size[0] GT 2) THEN BEGIN ; hope a reform makes it 2D
      im_aaaaaaa = REFORM(im_aaaaaaa) & im_size = size(im_aaaaaaa)
    ENDIF
    IF (im_size[0] GT 2) THEN MESSAGE,'Array must be 2D or 1D'
    x_aaaaaaa = DINDGEN(im_size[1])
    y_aaaaaaa = DINDGEN(im_size[2])
    dx_aaaaaaa  =0.5d
    dy_aaaaaaa = 0.5d
  ENDIF
  
  ; copy & rename variables so I can't accidentally alter input
  im = im_aaaaaaa
  x  = x_aaaaaaa
  dx = dx_aaaaaaa
  y  = y_aaaaaaa
  dy = dy_aaaaaaa
  
  ;im = [[0d,1d,2d,3d],[4d,5d,6d,7d],[8d,9d,10d,11d]] ; 2D data array
  ;x  = [1.5,2.5,3.5 ,4.5]  ; Center of the x-coords   (dim 0)
  ;dx = [0.5,0.5,0.25,0.5] ; Half-width of the x-bins (dim 0)
  ;y  = [1.5,2.5 ,3.5]     ; Center of the y-coords   (dim 1)
  ;dy = [0.5,0.25,0.4]     ; Half-width of the y-bins (dim 1)
  
  nx = N_ELEMENTS(x)
  ny = N_ELEMENTS(y)
  nim = SIZE(im,/DIMENSIONS)
  nxm1 = nx-1L
  nym1 = ny-1L
  IF (XYPLOT EQ 0) THEN BEGIN
    IF (N_ELEMENTS(nim) NE 2) THEN MESSAGE,'ERROR: Spectrogram expects to be fed a 2D array'
    IF (nx NE nim[0]) THEN MESSAGE,'ERROR: dx not the same size as dimension 1 of 2D array'
    IF (ny NE nim[1]) THEN MESSAGE,'ERROR: dy not the same size as dimension 2 of 2D array'
    ; if dx or dy is scalar, make array of right size
    ; check dx and dy size match
    IF (N_ELEMENTS(dx) EQ 1) THEN dx = DBLARR(nx)+dx $
    ELSE IF (N_ELEMENTS(dx) NE nx) THEN MESSAGE,'ERROR: dx not same size as x'
    IF (N_ELEMENTS(dy) EQ 1) THEN dy = DBLARR(ny)+dy $
    ELSE IF (N_ELEMENTS(dy) NE ny) THEN MESSAGE,'ERROR: dy not same size as y'
  ENDIF ELSE BEGIN ; if XYPLOT
    IF (TOTAL(nim) NE 0) THEN MESSAGE,'ERROR: Spectrogram XYPLOT expects to be empty DATA array'
    nxyplot = SIZE(y,/DIMENSIONS)
    IF (nx NE nxyplot[0]) THEN MESSAGE,'ERROR: Spectrogram XYPLOT x and y not the same size'
    ; if dx or dy is scalar, make array of right size
    ; check dx and dy size match
    IF (N_ELEMENTS(dx) EQ 1) THEN dx = DBLARR(nx)+dx $
    ELSE IF (N_ELEMENTS(dx) NE nx) THEN MESSAGE,'ERROR: dx not same size as x'
    IF (N_ELEMENTS(dy) EQ 1) THEN dy = DBLARR(nxyplot[0])+dy $
    ELSE IF ((N_ELEMENTS(dy) NE nxyplot[0]) AND (N_ELEMENTS(nxyplot) EQ 1)) THEN MESSAGE,'ERROR: dy not same size as y'
  ENDELSE
  
  
  
  IF KEYWORD_SET(XRANGE) THEN BEGIN
    IF (N_ELEMENTS(XRANGE) NE 2) THEN MESSAGE,'ERROR: XRANGE must be an array of size 2'
    IF (XRANGE[0] GE XRANGE[1]) THEN MESSAGE,'ERROR: XRANGE first value must be less than second'
    xminmax = XRANGE
  ENDIF ELSE BEGIN
    xminmax = [x[0]-dx[0], x[nxm1]+dx[nxm1]]; presume x is increasing
  ENDELSE
  IF KEYWORD_SET(YRANGE) THEN BEGIN
    IF (N_ELEMENTS(YRANGE) NE 2) THEN MESSAGE,'ERROR: YRANGE must be an array of size 2'
    IF (YRANGE[0] GE YRANGE[1]) THEN MESSAGE,'ERROR: YRANGE first value must be less than second'
    yminmax = YRANGE
  ENDIF ELSE BEGIN
    ymin = MIN(y,yminind,MAX=ymax,SUBSCRIPT_MAX=ymaxind);
    ;    yminmax = [y[0]-dy[0], y[nym1]+dy[nym1]]; presume x is increasing
    IF N_ELEMENTS(dy) EQ 1 THEN BEGIN
      yminmax = [y[yminind]-dy, y[ymaxind]+dy] ; presume x is increasing
    ENDIF ELSE BEGIN
      IF N_ELEMENTS(SIZE(y,/DIMENSIONS)) EQ 1 THEN BEGIN
        yminmax = [y[yminind]-dy[yminind], y[ymaxind]+dy[ymaxind]]; presume x is increasing
      ENDIF ELSE BEGIN
        IF N_ELEMENTS(SIZE(y,/DIMENSIONS)) EQ 2 THEN BEGIN
          dynn = SIZE(y,/DIMENSIONS)
          yminmax = [y[yminind]-dy[yminind MOD dynn[0]], y[ymaxind]+dy[ymaxind MOD dynn[0]]]; presume x is increasing
        ENDIF ELSE MESSAGE,'CONFUSING... see ROB' ; not coded fo this
      ENDELSE
    ENDELSE
    ; IF XYPLOT, IDL has a habbit of hiding lines on top/bottom, so add another delta
    IF ((XYPLOT EQ 1) AND (dy[0] NE 0))THEN yminmax = yminmax+[-dy[yminind],+dy[ymaxind]]
    IF ((XYPLOT EQ 1) AND (dy[0] EQ 0))THEN yminmax = yminmax+(yminmax[1]-yminmax[0])*0.01d*[-1d,1d]
  ENDELSE
  
  ; Log Data if requested, CLIM log done above
  IF (LogData EQ 1) THEN BEGIN
    DATA_min = MIN(im, data_min_index)
    IF DATA_min GT 0 THEN BEGIN
      im = ALOG10(im)
    ENDIF ELSE BEGIN ;  IF (DATA_min LE 0) THEN BEGIN
      zero_ind = WHERE(im EQ 0) ; must be at lest one
      print,'Warning: Logging data where '+STRTRIM(STRING(N_ELEMENTS(zero_ind)),2)+' elements of '$
        +STRTRIM(STRING(N_ELEMENTS(im)),2)+' are zero pre-Log (will be gray)'
      IF (im[data_min_index] LT 0) THEN BEGIN
        zero_neg = WHERE(im EQ 0) ; must be at lest one
        print,'Warning: Logging data where '+STRTRIM(STRING(N_ELEMENTS(zero_neg)),2)+' elements of '$
          +STRTRIM(STRING(N_ELEMENTS(im)),2)+' are negative pre-Log (will be gray)'
        zero_ind = [zero_ind,zero_neg]
      ENDIF
      im[zero_ind]=1e40 ; something positive so ALOG10 doesn't moan
      im = ALOG10(im)
      im[zero_ind]=!Values.D_NAN ; something positive so ALOG10 doesn't moan
    ENDELSE
  ENDIF
  
  
  ; Borrowing the Colormap of the default Matlab colormap:
  ; is Matlab default colormap, from 0-1, but IDL uses 0 to 255, hence the multiply by 255
  colormap = 255.0*[[0.0000,    0.0000,    0.5625],[0.0000,    0.0000,    0.6250],[0.0000,    0.0000,    0.6875],[0.0000,    0.0000,    0.7500],[0.0000,    0.0000,    0.8125],[0.0000,    0.0000,    0.8750],[0.0000,    0.0000,    0.9375],[0.0000,    0.0000,    1.0000],[0.0000,    0.0625,    1.0000],[0.0000,    0.1250,    1.0000],[0.0000,    0.1875,    1.0000],[0.0000,    0.2500,    1.0000],[0.0000,    0.3125,    1.0000],[0.0000,    0.3750,    1.0000],[0.0000,    0.4375,    1.0000],[0.0000,    0.5000,    1.0000],[0.0000,    0.5625,    1.0000],[0.0000,    0.6250,    1.0000],[0.0000,    0.6875,    1.0000],[0.0000,    0.7500,    1.0000],[0.0000,    0.8125,    1.0000],[0.0000,    0.8750,    1.0000],[0.0000,    0.9375,    1.0000],[0.0000,    1.0000,    1.0000],[0.0625,    1.0000,    0.9375],[0.1250,    1.0000,    0.8750],[0.1875,    1.0000,    0.8125],[0.2500,    1.0000,    0.7500],[0.3125,    1.0000,    0.6875],[0.3750,    1.0000,    0.6250],[0.4375,    1.0000,    0.5625],[0.5000,    1.0000,    0.5000],[0.5625,    1.0000,    0.4375],[0.6250,    1.0000,    0.3750],[0.6875,    1.0000,    0.3125],[0.7500,    1.0000,    0.2500],[0.8125,    1.0000,    0.1875],[0.8750,    1.0000,    0.1250],[0.9375,    1.0000,    0.0625],[1.0000,    1.0000,    0.0000],[1.0000,    0.9375,    0.0000],[1.0000,    0.8750,    0.0000],[1.0000,    0.8125,    0.0000],[1.0000,    0.7500,    0.0000],[1.0000,    0.6875,    0.0000],[1.0000,    0.6250,    0.0000],[1.0000,    0.5625,    0.0000],[1.0000,    0.5000,    0.0000],[1.0000,    0.4375,    0.0000],[1.0000,    0.3750,    0.0000],[1.0000,    0.3125,    0.0000],[1.0000,    0.2500,    0.0000],[1.0000,    0.1875,    0.0000],[1.0000,    0.1250,    0.0000],[1.0000,    0.0625,    0.0000],[1.0000,    0.0000,    0.0000],[0.9375,    0.0000,    0.0000],[0.8750,    0.0000,    0.0000],[0.8125,    0.0000,    0.0000],[0.7500,    0.0000,    0.0000],[0.6875,    0.0000,    0.0000],[0.6250,    0.0000,    0.0000],[0.5625,    0.0000,    0.0000],[0.5000,    0.0000,    0.0000]]
  ; Number of different colors
  ncolormap = N_ELEMENTS(colormap)/3L; divide by 3 as n by 3
  ncolormapm1 = ncolormap-1
  ; may need colormap for mupltiple lines
  
  
  IF (XYPLOT EQ 0) THEN BEGIN
  
    ; auto set color limit
    IF KEYWORD_SET(CLIM) THEN BEGIN
      IF (N_ELEMENTS(CLIM) NE 2) THEN MESSAGE,'ERROR: CLIM must be an array of size 2'
      IF (CLIM[0] GE CLIM[1]) THEN MESSAGE,'ERROR: CLIM first value must be less than second'
    ENDIF ELSE BEGIN
      CLIM = [MIN(im, MAX = CLIMz), 0.0]; work out min and max at once
      CLIM[1] = CLIMz
      ;delvar,CLIMz ; not needed anymore
      IF (CLIM[0] EQ CLIM[1]) THEN CLIM[1] = CLIM[0]+1; MESSAGE,'ERROR: Data only has one value!'
    ENDELSE
    ; Log Data if requested
    IF (LogData EQ 1) THEN BEGIN
      IF CLIM[0] EQ 0 THEN CLIM[0] = MIN([1,MIN(im(WHERE(im GT 0)))])
      ;CLIM = ALOG10(CLIM) ; not needed as im is already logged
    ENDIF
    
    ; Work out what index of the colortable each value would have
    cindex = FIX((ncolormap-1L)*(im-CLIM[0])/(CLIM[1]-CLIM[0]))
    ; Now remove ones outside of CLIM range:
    cindex[WHERE((cindex LT 0L) OR (cindex GT ncolormapm1),/NULL)]= -1
    
    
  ENDIF
  
  ; Check that the Log 10 didn't give any NaN's nor Infs
  IF ((LogData EQ 1) AND (XYPLOT EQ 0))  THEN cindex[WHERE(FINITE(im) EQ 0,/NULL)]= -1
  
  ; Set position as [x_lower, y_lower, x_upper, y_upper], ; Use Normalized co-ords, 0-1.
  position_panel    = [0.15  , 0.1  , 0.8  , 0.9  ]
  position_colorbar = [0.85 , 0.1  , 0.9  , 0.9  ]
  
  
  
  ; Check if multiple windows are open:
  w = getwindows() ; if no windows still returns null object of size 1
  IF w[0] THEN BEGIN ; must be at least one, use w[0] in case no windows is null
    ; Now use last Window, and erase it
    nw = N_ELEMENTS(w);
    w[nw-1L].SETCURRENT
    w[nw-1L].ERASE
  ENDIF
  OBJ_DESTROY,w
  
  IF (XYPLOT EQ 1) THEN BEGIN
    ; If nothing to plot, clear figure above and exit
    if (N_ELEMENTS(x) EQ 0) THEN return,0
    ; draw diagonal 'invisible' line to set up figure
    IF (N_ELEMENTS(nxyplot) EQ 1) THEN BEGIN
      p = PLOT(x, y, COLOR = 'b', POSITION = position_panel, /CURRENT, OVERPLOT = 1 )
    ENDIF ELSE BEGIN
      CLIM = [0L,(nxyplot[1]-1L)]
      FOR nxyplotz = CLIM[0],CLIM[1] DO BEGIN
        ;Work out what index of the colortable each value would have
        cindex = FIX((ncolormap-1L)*(nxyplotz-CLIM[0])/(CLIM[1]-CLIM[0]))
        p = PLOT(x, y[*,nxyplotz], COLOR = colormap[*,cindex], POSITION = position_panel, /CURRENT, OVERPLOT = 1 )
      ENDFOR
    ENDELSE
  ENDIF ELSE BEGIN
    ; draw diagonal 'invisible' line to set up figure
    p = PLOT(xminmax, yminmax, LINESTYLE = 6, POSITION = position_panel, /CURRENT, OVERPLOT = 1 ) ; linestyle 6 is no line, so does nto show.
  ENDELSE
  p.XRANGE = xminmax
  p.YRANGE = yminmax
  IF XTIME THEN BEGIN
    ; Only do 3 time stamps
    p.XTICKINTERVAL = (xminmax[1]-xminmax[0])/3d
    IF (xminmax[1]-xminmax[0]) GE 1 THEN $
      DUMMY = LABEL_DATE(DATE_FORMAT=['%Y-%N-%DT%H:%I:%S']) $
    ELSE $
      DUMMY = LABEL_DATE(DATE_FORMAT=['%H:%I:%S'])
    p.XTICKFORMAT = 'LABEL_DATE'
  ENDIF
  
  ; Draw Spectrogram
  IF (XYPLOT EQ 0) THEN BEGIN
  
    ;for zx = 0L,nxm1 DO $
    ;  for zy = 0L,nym1 DO $
    ;    if (cindex[zx,zy] NE -1) THEN $  ; only plot ploygone if in range, if always in range (as above with CLIM set by min/max) then IF statement is not needed.
    ;      poly = POLYGON(x(zx)+dx(zx)*[-1,1,1,-1], y(zy)+dy(zy)*[-1,-1,1,1], /DATA, /FILL_BACKGROUND, FILL_COLOR = colormap[*,cindex[zx,zy]], LINESTYLE = '')
    ; Color can be RGB where it's a vector of size 3, values 0 to 255
    
    polyx = DBLARR(nxm1+1,4,/NOZERO)
    for zx = 0L,nxm1 DO polyx[zx,*] = x(zx)+dx(zx)*[-1d,1d,1d,-1d]
    polyx = TRANSPOSE(polyx)
    polyy = DBLARR(nym1+1,4,/NOZERO)
    for zy = 0L,nym1 DO polyy[zy,*] = y(zy)+dy(zy)*[-1d,-1d,1d,1d]
    polyy = TRANSPOSE(polyy)
    
    FOR cBin =-1L,ncolormapm1 DO BEGIN
      ind = WHERE(cindex EQ cBin, nind)
      IF (ind[0] EQ -1) THEN CONTINUE
      xx = DBLARR(nind*4L)
      yy = xx ; DBLARR(nind*4L)
      cc = LONARR(nind*5L)
      cnt  = 0L
      cntc = 0L
      for zx = 0L,nxm1 DO BEGIN
        for zy = 0L,nym1 DO BEGIN
          if (cindex[zx,zy] EQ cBin) THEN BEGIN  ; only plot ploygone if in range, if always in range (as above with CLIM set by min/max) then IF statement is not needed.
            xx[cnt +[0,1,2,3]]   = polyx[*,zx]
            yy[cnt +[0,1,2,3]]   = polyy[*,zy]
            cc[cntc+[0,1,2,3,4]] = [4,cnt+[0,1,2,3]]
            cnt  = cnt  +4L ; must be after cc line above!
            cntc = cntc +5L
          endif
        endfor
      ENDFOR
      ;  print,cnt
      c_color = colormap[*,cBin]
      if (cBin EQ -1L) THEN c_color = 255*[0.8, 0.8, 0.8] ; grey for fill when not finite or outside of range
      poly = POLYGON(xx, yy, CONNECTIVITY = cc, /DATA, /FILL_BACKGROUND, FILL_COLOR = c_color, LINESTYLE = '')
      
    ENDFOR
  ENDIF ; draw spectrogram
  
  
  ; Draw COLORBAR for spectorgram or many XY plot lines
  is_colorbar = 0
  IF ((XYPLOT EQ 0 ) OR (N_ELEMENTS(nxyplot) GT 1)) THEN BEGIN
    is_colorbar = 1
    ; Add a colorbar
    c = COLORBAR(RGB_TABLE = colormap, RANGE = CLIM, ORIENTATION = 1, POSITION = position_colorbar)
    ; Change some colrobar properties
    ; ORIENTATION = 0 for horizontal, 1 vertical
    ; RANGE is range of colors, so CLIM.
    c.TEXTPOS = 1 ; 0 on left, 1 on right
    c.TICKDIR = 1
    c.BORDER_ON = 1
    ;c.COLOR = 'Black' ; default is black
    c.FONT_STYLE = 'Italic'  ; 'Normal','Bold','Italic' or 'Bold Italic'
    c.FONT_SIZE = 10  ; Font for title and numbers, default is 16
    ;c.TITLE='This is a colorbar title'
    ;IF (LogData EQ 1) THEN BEGIN
    ;    colorTicks = c.TICKVALUES
    ;    c.TICKNAME = STRTRIM(STRING(10^colorTicks))
    ;    stop
    ;ENDIF
    
  ENDIF
  
  ; XLABEL/YLABLE are matlab things, so duplicating XTITLE/YTITLE so Rob doesn't get confused
  IF KEYWORD_SET(XLABEL) THEN XTITLE = XLABEL
  IF KEYWORD_SET(YLABEL) THEN YTITLE = YLABEL
  IF KEYWORD_SET(CLABEL) THEN CTITLE = CLABEL
  
  IF KEYWORD_SET(XTITLE) THEN p.XTITLE = XTITLE ;ELSE p.XTITLE = 'x-axis'
  IF KEYWORD_SET(YTITLE) THEN p.YTITLE = YTITLE ;ELSE p.YTITLE = 'y-axis'
  IF KEYWORD_SET(TITLE)  THEN p.TITLE  = TITLE  ;ELSE p.TITLE = 'Plot Title'
  IF KEYWORD_SET(CTITLE) THEN BEGIN
    IF (N_ELEMENTS(CTITLE) NE 0) AND (is_colorbar EQ 1) THEN BEGIN
      IF (LogData EQ 1) THEN ctitle = 'LOG10[ '+ctitle+' ]
      c.TITLE  = CTITLE ;ELSE c.TITLE = 'Colorbar Title'
    ENDIF
  ENDIF
  
  IF ((XYPLOT EQ 0 ) OR (N_ELEMENTS(nxyplot) GT 1)) THEN BEGIN
    handles = CREATE_STRUCT('Plot', p,'Colorbar',c) ; doesn't seem to matter which save
  ENDIF ELSE BEGIN
    handles = CREATE_STRUCT('Plot', p)
  ENDELSE
  RETURN,handles
END