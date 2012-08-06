;+
; NAME:
;       FIT_SEDS_DEG2DEC
;
; PURPOSE:
;       Convert degrees to DEC string.
;
; CATEGORY:
;       Bayesian fitting.
;
; CALLING SEQUENCE:
;       dec_string = fit_seds_deg2ra(deg)
;
; INPUTS:
;       deg : dec in degrees (double)
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; PROCEDURE:
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
; 	Written     : Karl Gordon (long, long ago)
;        6 Aug 2012 : documentation (KDG)
;-

function fit_seds_deg2dec,deg,dec_split=dec_split,sep_str=sep_str, $
  trunc_sec=trunc_sec

if (n_elements(sep_str) EQ 0) then sep_str = ':'

dec = deg

if (dec LT 0.0) then begin
   dec_sgn = '-'
endif else begin
   dec_sgn = '+'
endelse
dec = abs(dec) 
decd = fix(dec)
decm = fix(60.0*(dec - decd))
if (keyword_set(trunc_sec)) then begin
    decs = fix(10*60.0*(60.0*(dec - decd) - decm))/10.
endif else begin
    decs = round(10*60.0*(60.0*(dec - decd) - decm))/10.
endelse

dec_split = [decd,decm,decs]

if (decd LE 9) then begin
   dec_str = '0' + strtrim(string(decd),2) + sep_str
endif else begin
   dec_str = strtrim(string(decd),2) + sep_str
endelse
if (decm LE 9) then begin
   dec_str = dec_str + '0' + strtrim(string(decm),2) + sep_str
endif else begin
   dec_str = dec_str + strtrim(string(decm),2) + sep_str
endelse
if (decs LT 10) then begin
   dec_str = dec_str + '0' + strtrim(string(decs,format="(F6.1)"),2)
endif else begin
   dec_str = dec_str + strtrim(string(decs,format="(F6.1)"),2)
endelse

dec_str = dec_sgn + dec_str

return,dec_str

end
