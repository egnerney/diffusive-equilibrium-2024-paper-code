
pro run_make_line_ist_info_neutral_Oxygen

  ; set variable to path of main chianti directory, assuming chianti is in the current working directory
  cd, current = current
  chiantipath = current + '/chianti/'

  ; or, if chianti is located elsewhere, edit the following variable
  ;chiantipath = '/Users/shinn/pro/chianti/' ; <---- user specific path

  ; use 'exists' keyword to determine whether to execute the following
  defsysv, '!chianti', exists = compile
  if ~compile then defsysv, '!chianti', 1

  if ~compile then defsysv, '!rjkm', 71492d ; RJ in km
  if ~compile then defsysv, '!kev', 8.617385d-5 ; Boltzmann's constant in eV
  if ~compile then defsysv, '!rootpi', sqrt(!dpi)

  ; create a string of all directories under chianti folder separated by ':'s
  if ~compile then chiantipros = ':' + expand_path('+' + chiantipath + 'idl/')

  ; if compile script has not been run before, add chianti path
  if ~compile then !path = !path + chiantipros
  if ~compile then use_chianti, chiantipath + 'dbase'

  ; compile the widget
  resolve_routine, 'make_line_ist_info_neutral_Oxygen'

  ; run the widget
 make_line_ist_info_neutral_Oxygen

end
