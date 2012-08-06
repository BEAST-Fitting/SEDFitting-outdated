;+
; NAME:
;       FIT_SEDS_DEG2RA
;
; PURPOSE:
;       Convert degrees to RA string.
;
; CATEGORY:
;       Bayesian fitting.
;
; CALLING SEQUENCE:
;       ra_string = fit_seds_deg2ra(ra_deg)
;
; INPUTS:
;       ra_deg : ra in degrees (double)
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

function fit_seds_deg2ra,deg,ra_split=ra_split,sep_str=sep_str, $
                trunc_sec=trunc_sec

if (n_elements(sep_str) EQ 0) then sep_str = ':'

ra = deg/15.0

rah = fix(ra)
ram = fix(60.0*(ra - rah))
if (keyword_set(trunc_sec)) then begin
    ras = fix(100*60.0*(60.0*(ra - rah) - ram))/100.
endif else begin
    ras = round(100*60.0*(60.0*(ra - rah) - ram))/100.
endelse

ra_split = [rah,ram,ras]

if (rah LE 9) then begin
   ra_str = '0' + strtrim(string(rah),2) + sep_str
endif else begin
   ra_str = strtrim(string(rah),2) + sep_str
endelse
if (ram LE 9) then begin
   ra_str = ra_str + '0' + strtrim(string(ram),2) + sep_str
endif else begin
   ra_str = ra_str + strtrim(string(ram),2) + sep_str
endelse
if (ras LT 10) then begin
   ra_str = ra_str + '0' + strtrim(string(ras,format="(F6.2)"),2)
endif else begin
   ra_str = ra_str + strtrim(string(ras,format="(F6.2)"),2)
endelse

return,ra_str

end
