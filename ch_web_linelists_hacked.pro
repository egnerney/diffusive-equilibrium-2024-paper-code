
PRO ch_web_linelists_hacked 


;+
; NAME:
;     CH_WEB_LINELISTS
;
; PURPOSE:
;     Generates the set of line list tables that are linked to the CHIANTI
;     website. Only the latex files are created; the user needs to run
;     pdflatex to generate the pdf files.
;
;     The line lists are created for the "extended-flare" DEM for a
;     pressure of 10^16 and for photospheric abundances. Intensity limits
;     are applied to reduce the number of lines in the tables. Population
;     lookup tables are used to reduce computation time.
;
; CATEGORY:
;     CHIANTI; line lists.
;
; CALLING SEQUENCE:
;     CH_WEB_LINELISTS
;
; INPUTS:
;      None.
;
; OUTPUTS:
;      Creates the following six latex files in the working directory.
;       ch_line_list_v10.0.1_1_50.tex
;       ch_line_list_v10.0.1_50_150.tex
;       ch_line_list_v10.0.1_150_912.tex
;       ch_line_list_v10.0.1_912_2000.tex
;       ch_line_list_v10.0.1_2000_10000.tex
;       ch_line_list_v10.0.1_10000_600000.tex
;
; EXAMPLE:
;      IDL> ch_web_linelists
;
; MODIFICATION HISTORY:
;      Ver.1, 08-Jun-2022, Peter Young
;      Ver.2, 31-May-2023, Peter Young
;        Chnaged abundance file to !abund_file.
;-



ch_ver=ch_get_version()
press=1e16

;
; The DEM is hard-coded to the extended flare DEM.
;
dem_name=concat_dir(!xuvtop,'dem')
dem_name=concat_dir(dem_name,'flare_ext.dem')
chck=file_info(dem_name)
IF chck.exists EQ 0 THEN BEGIN
  message,/info,/cont,'The DEM file was not found. Returning...'
  return
ENDIF 

abund_name=!abund_file
chck=file_info(abund_name)
IF chck.exists EQ 0 THEN BEGIN
  message,/info,/cont,'The abundance file was not found. Returning...'
  return
ENDIF 

;6302.05
w0=6302
w1=6303
mini=1.16e6/1e9 ; was over 1e3 before
outfile='egnerneyhack_ch_line_list_v'+ch_ver+'_'+trim(w0)+'_'+trim(w1)+'.tex'
latex_wvl_dem,w0,w1,outfile=outfile,dem_name=dem_name,abund_name=abund_name,/all, $
              mini=mini,pressure=press,/lookup


;1302.1680       1304.8580       1306.0291 Angstroms
w0=1301
w1=1307
mini=2.28e6/1e9
outfile='egnerneyhack_ch_line_list_v'+ch_ver+'_'+trim(w0)+'_'+trim(w1)+'.tex'
latex_wvl_dem,w0,w1,outfile=outfile,dem_name=dem_name,abund_name=abund_name,/all, $
              mini=mini,pressure=press,/lookup


;1355.5980       1358.5120 Angstroms

w0=1354
w1=1359
mini=1.27e6/1e9
outfile='egnerneyhack_ch_line_list_v'+ch_ver+'_'+trim(w0)+'_'+trim(w1)+'.tex'
latex_wvl_dem,w0,w1,outfile=outfile,dem_name=dem_name,abund_name=abund_name,/all, $


END
