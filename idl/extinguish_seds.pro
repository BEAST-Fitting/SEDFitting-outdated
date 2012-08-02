;+
; NAME:
;       EXTINGUISH_SEDS
;
; PURPOSE:
;       Extinguish stellar SEDs for a range of A(V) values assuming a
;       dust extinction curve with a particular value of R(V).
;
; CATEGORY:
;       Bayesian fitting.
;
; CALLING SEQUENCE:
;       EXTINGUISH_SEDS, grid_seds, ext_grid_seds
;
; INPUTS:
;       grid_seds : grid of stellar SEDs (array of structures)
;
; KEYWORD PARAMETERS:
;       av_step : step size in A(V) [default 0.025]
;       av_min : min A(V) value [default 0.0]
;       av_max : max A(V) value [default 1.0]
;       rv : R(V) value [default 3.1]
;       silent : set to turn of screen output
;
; OUTPUTS:
;       ext_grid_seds : grid of stellar SEDs (3D array of structures)
;
; OPTIONAL OUTPUTS:
;
; PROCEDURE:
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
; 	Started     : Karl Gordon (2011)
;       2011-2012   : lots of development (undocumented, KDG)
;       2 Aug 2012  : Cleaned up and full documentation added (KDG)
;-

pro extinguish_seds,grid_seds,ext_grid_seds,rv=rv,max_av=max_av,min_av=min_av,av_step=av_step

; defaults
if (not keyword_set(rv)) then rv = 3.1
if (not keyword_set(max_av)) then max_av = 1.0
if (not keyword_set(min_av)) then min_av = 0.0
if (not keyword_set(av_step)) then av_step = 0.025

; setup the A(V) array
n_av = round((max_av - min_av)/av_step) + 1
many_av = min_av + findgen(n_av)*av_step

; all seds assumed to have the same wavelength scale
; get the A(lambda)/A(V) extinction curve 
; assuming the standard CCM89 R(V) dependent curve
ext_curve = ccm(rv,1.e4/grid_seds[0,0].waves)

; setup the output 3D structure of seds
size_grid_seds = size(grid_seds)
n_wave = n_elements(grid_seds[0,0].waves)
red_sed = {logt: 0.0, logg : 0.0, logl: 0D0, av: 0.0, rv: 0.0, $
           st_logt: 0.0, st_logg : 0.0, st_logl: 0D0, st_radius: 0.0, $
           st_mass: 0.0, st_bmass: 0.0, st_age: 0.0, st_tag: 0, $
           waves: fltarr(n_wave), $
           fluxes: dblarr(n_wave), $
           total_flux: 0D0}
ext_grid_seds = replicate(red_sed,size_grid_seds[1],size_grid_seds[2],n_av)

; loop over the A(V) values and extinguish away!
for i = 0,(n_av-1) do begin
    ext_grid_seds[*,*,i] = grid_seds

    ; save dust extinction values
    ext_grid_seds[*,*,i].av = many_av[i]
    ext_grid_seds[*,*,i].rv = rv

    for j = 0,(size_grid_seds[1]-1) do begin
        for k = 0,(size_grid_seds[2]-1) do begin
            ext_grid_seds[j,k,i].fluxes *= 10^(-0.4*ext_curve*many_av[i])
        endfor
    endfor

endfor

end
