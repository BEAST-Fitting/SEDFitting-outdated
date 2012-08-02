;+
; NAME:
;       CCM
;
; PURPOSE:
;       Procedure to compute the one parameter (R_v) extinction curve.  This
;       is from the work of Cardelli, Clayton, and Mathis 1989, ApJ, 345, 245.
;
; CATEGORY:
;       Bayesian fitting.
;
; CALLING SEQUENCE:
;       result = CCM(rv,x)
;
; INPUTS:
;       rv : R(V) = A(V)/E(B-V) value [rough measure of average grain size]
;       x : wavelength vector in units of 1/micron
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;       result : A(lambda)/A(V) as a function of wavelength
;
; OPTIONAL OUTPUTS:
;       a : zero point coefficient of R(V) dependent relationship
;       b : slope coefficient of R(V) dependent relationship
;
; PROCEDURE:
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
; 	Started     : Karl Gordon (sometime in the distant past)
;       28 Dec 2005 : updated to handle a single x or a vector of xs (KDG)
;       2 Aug 2012  : Cleaned up and full documentation added (KDG)
;-

function ccm,rv,x,a,b

npts = n_elements(x)
a = fltarr(npts)
b = fltarr(npts)
for i = 0L,(npts-1) do begin
    if ((x[i] GE 0.01) AND (x[i] LT 1.1)) then begin
        a[i] = 0.574*(x[i]^1.61)
        b[i] = -0.527*(x[i]^1.61)
    endif else if ((x[i] GE 1.1) AND (x[i] LT 3.3)) then begin
        y = x[i] - 1.82
        a[i] = 1.0 + 0.17699*y - 0.50447*y^2 - 0.02427*y^3 + 0.72085*y^4 $
           + 0.01979*y^5 - 0.77530*y^6 + 0.32999*y^7
        b[i] = 1.41338*y + 2.28305*y^2 + 1.07233*y^3 - 5.38434*y^4 $
           - 0.62251*y^5 + 5.30260*y^6 - 2.09002*y^7
    endif else if ((x[i] GE 3.3) AND (x[i] LE 8.0)) then begin
        if ((x[i] GE 5.9) AND (x[i] LE 8.0)) then begin
            Fa = -0.04473*(x[i] - 5.9)^2 - 0.009779*(x[i] - 5.9)^3
            Fb = 0.2130*(x[i] - 5.9)^2 + 0.1207*(x[i] - 5.9)^3
        endif else begin
            Fa = 0.0
            Fb = 0.0
        endelse
        a[i] = 1.752 - 0.316*x[i] - 0.104/((x[i] - 4.67)^2 + 0.341) + Fa
        b[i] = -3.090 + 1.825*x[i] + 1.206/((x[i] - 4.62)^2 + 0.263) + Fb
    endif else if ((x[i] GT 8.0) AND (x[i] LE 20.0)) then begin
        a[i] = -1.073 - 0.628*(x[i] - 8.0) + 0.137*(x[i] - 8.0)^2 - 0.070*(x[i] - 8.0)^3
        b[i] = 13.670 + 4.257*(x[i] - 8.0) - 0.420*(x[i] - 8.0)^2 + 0.374*(x[i] - 8.0)^3
    endif else begin
        print,'ccm error: outside the bounds of the CCM relation'
        print,x[i]
    endelse
endfor

temp = a + b/Rv

return,temp

end
