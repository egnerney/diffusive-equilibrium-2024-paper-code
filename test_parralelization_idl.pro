pro solver_call, i, n, T, result
  ; This is the task we want to run in parallel
  result = n(i)*T(i)
  
  SAVE, result, FILENAME='result_' + STRTRIM(index,2) + '.sav'
END

pro test_parralelization_IDL
tic

n_elements = 501
x = 0.0005*findgen(n_elements ) + 5.

  n = 2000.*Exp(-((x - 5.7)/0.15)^2)
  
  b = 4.32193
  
  T = 1.*(x/5.)^b
  
  result_forloop = fltarr(n_elements )
  for j = 0, n_elements - 1 do begin
  
  result_forloop = n(j)*T(j)
    
  endfor
toc

tic

indices = LINDGEN(n_elements) ; Generate an array [0, 1, 2, ..., 100000]

; Loop and initiate separate IDL sessions
FOR i=0, n_elements-1 DO BEGIN
  command = 'idl -e "solver_call, ' + STRTRIM(indices[i], 2) + ', ' + STRTRIM(n, 2) + ', ' + STRTRIM(T, 2)  + '" &'
  SPAWN, command
ENDFOR

 toc

stop
end
