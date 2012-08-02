;+
; NAME:
;       GET_STELLAR_TRACKS
;
; PURPOSE:
;       Read in stellar evolutionary tracks.  Currently a hard coded
;       set of interpolated Geneva based tracks.
;
; CATEGORY:
;       Bayesian fitting.
;
; CALLING SEQUENCE:
;       GET_STELLAR_TRACKS, grid_logl, grid_logt, grid_logg,
;                           grid_radius, grid_mass, grid_bmass, 
;                           grid_tag, grid_age
;
; INPUTS:
;
; KEYWORD PARAMETERS:
;       no_iter : turn off the interpolation by a factor of 25 in mass
;                 and track step
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
; 	Started     : Karl Gordon (23 Jul 2012)
;-

pro get_stellar_tracks,grid_logl,grid_logt,grid_logg,grid_radius,grid_mass,grid_bmass,grid_tag,grid_age, $
  no_inter=no_inter

met_path = '~/Pro/Fit_SEDs/Stellar_Tracks/Geneva/z0.020/'

n_tables = 22
for i = 0,(n_tables-1) do begin
    readcol,met_path+'table'+strtrim(string(i+1),2),step,time,mass,logl,logt,/silent
    ; calculate log(g)
    radius2 = (10^logl)*3.839d26/(4.0*!PI*5.67d-8*((10^logt)^4))
    g = 6.67d-11*(mass*1.989d30)/radius2
    logg = alog10(g*1e2)

    if (i EQ 0) then begin
        n_steps = n_elements(logl)
        grid_logl = fltarr(n_steps,n_tables)
        grid_logt = grid_logl
        grid_logg = grid_logl
        grid_radius = grid_logl
        grid_mass = grid_logl
        grid_bmass = grid_logl
        grid_age = grid_logl
        grid_tag = grid_logl
    endif
    n_cur_steps = n_elements(logl)
    grid_logl[0:n_cur_steps-1,i] = logl
    grid_logt[0:n_cur_steps-1,i] = logt
    grid_logg[0:n_cur_steps-1,i] = logg
    grid_radius[0:n_cur_steps-1,i] = sqrt(radius2)/6.955d8
    grid_mass[0:n_cur_steps-1,i] = mass
    grid_bmass[0:n_cur_steps-1,i] = mass[0]
    grid_age[0:n_cur_steps-1,i] = time
    grid_tag[0:n_cur_steps-1,i] = (i+1)
    if (n_cur_steps LT n_steps) then begin
        grid_logl[n_cur_steps:n_steps-1,i] = replicate(!values.f_nan,n_steps-n_cur_steps)
        grid_logt[n_cur_steps:n_steps-1,i] = replicate(!values.f_nan,n_steps-n_cur_steps)
        grid_logg[n_cur_steps:n_steps-1,i] = replicate(!values.f_nan,n_steps-n_cur_steps)
        grid_radius[n_cur_steps:n_steps-1,i] = replicate(!values.f_nan,n_steps-n_cur_steps)
        grid_mass[n_cur_steps:n_steps-1,i] = replicate(!values.f_nan,n_steps-n_cur_steps)
        grid_bmass[n_cur_steps:n_steps-1,i] = replicate(!values.f_nan,n_steps-n_cur_steps)
        grid_age[n_cur_steps:n_steps-1,i] = replicate(!values.f_nan,n_steps-n_cur_steps)
        grid_tag[n_cur_steps:n_steps-1,i] = replicate(!values.f_nan,n_steps-n_cur_steps)
    endif
endfor

; treat the table as an image and make it bigger to reduce
; interpolation errors when matching to stellar atmosphere models
if (not keyword_set(no_inter)) then begin
    n_inter = 25
    x = findgen((n_inter*n_steps)-(n_inter-1))/n_inter
    y = findgen((n_inter*n_tables)-(n_inter-1))/n_inter
    grid_logl = interpolate(grid_logl,x,y,/grid,missing=!values.f_nan,cubic=-0.5)
    grid_logt = interpolate(grid_logt,x,y,/grid,missing=!values.f_nan,cubic=-0.5)
    grid_logg = interpolate(grid_logg,x,y,/grid,missing=!values.f_nan,cubic=-0.5)
    grid_radius = interpolate(grid_radius,x,y,/grid,missing=!values.f_nan,cubic=-0.5)
    grid_mass = interpolate(grid_mass,x,y,/grid,missing=!values.f_nan,cubic=-0.5)
    grid_bmass = interpolate(grid_bmass,x,y,/grid,missing=!values.f_nan,cubic=-0.5)
    grid_age = interpolate(grid_age,x,y,/grid,missing=!values.f_nan,cubic=-0.5)
    grid_tag = interpolate(grid_tag,x,y,/grid,missing=!values.f_nan,cubic=-0.5)
endif

end
