;+
; NAME:
;       GET_GRID_STELLAR_SEDS
;
; PURPOSE:
;       Creates a grid of single star SEDs at an input distance.  This
;       uses a stellar atmosphere grid with a set of stellar
;       evolutionary tracks.  The currently, these are the Castelli &
;       Kurucz stellar atmospheres and the Genenva stellar
;       evolutionary tracks.
;
; CATEGORY:
;       Bayesian fitting.
;
; CALLING SEQUENCE:
;       GET_GRID_STELLAR_SEDS, distance, grid_seds
;
; INPUTS:
;       distance : distance to galaxy in pc [default is M31, 776e3 pc]
;
; KEYWORD PARAMETERS:
;       no_inter : no interpolation for stellar evolutionary tracks
;                  by default, the tracks are interpolated to higher
;                  mass and time resolution
;       silent : set to turn of screen output
;
; OUTPUTS:
;       grid_seds : grid of stellar SEDs (2D array of structures)
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

pro get_grid_stellar_seds,distance,grid_seds,no_inter=no_inter,silent=silent

; get the stellar evolutionary tracks
get_stellar_tracks,logl,logt,logg,radius,mass,bmass,gtag,age,no_inter=no_inter
size_logl = size(logl)

; get the filenames for the stellar atmospheres
; current uses a custom set of interpolated temperature and gravity
; models (solar metallicity only)
files = file_search('~/Pro/Fit_SEDs/Stellar_Grids/ck04/p00/*_unitemps.moregravs.fits')
n_files = n_elements(files)

; now get the effective temperatures of the models (from filenames)
us_pos = strpos(files,'_')
dot_pos = strpos(files,'.',5)
temp_val = fltarr(n_files)
for i = 0,(n_files-1) do begin
    tus_pos = strpos(files[i],'_',us_pos[i]+1)
    tus_pos = strpos(files[i],'_',tus_pos+1)
    temp_val[i] = float(strmid(files[i],tus_pos+1,dot_pos[i]-tus_pos-10))
endfor
sindxs = sort(temp_val)
files = files[sindxs]
temp_val = temp_val[sindxs]

; final variables
temps = temp_val
n_temps = n_elements(temps)

; gravities are hard coded
gravs = findgen(51)*0.1
n_gravs = n_elements(gravs)

; read in the stellar atmospheres and put them in the output structure
for i = 0,(n_files-1) do begin
    t = mrdfits(files[i],1,/silent)

    ; currently only saving wavelengths that are of interest
    ;  starts at 912 A (may want to change this to provide predictions of ly-continuum photons)
    ;  changed to extend to 30 microns (useful for IRAC/MIPS24 and JWST predictions)
    windxs = where((t.wavelength GT 912.) and (t.wavelength LT 30000.),n_twave)

    tags = tag_names(t) 
    n_tgravs = n_elements(tags) - 1
    if (i EQ 0) then begin
        n_wave = n_twave
        n_gravs = n_tgravs
        waves = fltarr(n_files,n_wave)
        fluxes = dblarr(n_files,n_gravs,n_wave)
        total_fluxes = dblarr(n_files,n_gravs)

        ; setup the grid seds structure
        sed = {logt: 0.0, logg : 0.0, logl: 0D0, av: 0.0, rv: 0.0, $
               st_logt: 0.0, st_logg : 0.0, st_logl: 0D0, st_radius: 0.0, $
               st_mass: 0.0, st_bmass: 0.0, st_age: 0.0, st_tag: 0, $
               waves: fltarr(n_wave), $
               fluxes: dblarr(n_wave), $
               total_flux: 0D0}
        grid_seds = replicate(sed,n_files,n_gravs)

    endif
    if (n_twave NE n_wave) then begin
        print,'***ERROR***'
        print,'current number of wavelengths different from 1st one'
        print,'current # waves = ', n_twave
        print,'    1st # waves = ', n_wave
        stop
    endif
    if (n_tgravs NE n_gravs) then begin
        print,'***ERROR***'
        print,'current number of gravities different from 1st one'
        print,'current # gravs = ', n_tgravs
        print,'    1st # gravs = ', n_grav
        stop
    endif

    for k = 0,(n_gravs-1) do begin
        
        if (total(t[*].(k+1)) GT 0.) then begin

            ; get the closest matching values from the isochrones for this model
            chisqr = ((alog10(temps[i]) - logt)/0.1)^2 + ((gravs[k] - logg)/0.1)^2
            min_chisqr = min(chisqr,indx)
            min_step_indx = indx mod size_logl[1]
            min_mass_indx = (indx/size_logl[1])

            ; by-hand optimized values for the allowed mismatch between
            ; the stellar atmosphere and stellar evolutionary track log(teff) & log(g)
            ; would be better to do a direct interpolation here (KDG 2 Aug 2012)
            if ((abs(logg[min_step_indx,min_mass_indx] - gravs[k]) LT 0.02) AND $
                (abs(logt[min_step_indx,min_mass_indx] - alog10(temps[i])) LT 0.02)) then begin

                ; save the stellar tracks numbers
                grid_seds[i,k].st_logt = logt[min_step_indx,min_mass_indx]
                grid_seds[i,k].st_logg = logg[min_step_indx,min_mass_indx]
                grid_seds[i,k].st_logl = logl[min_step_indx,min_mass_indx]
                grid_seds[i,k].st_radius = radius[min_step_indx,min_mass_indx]
                grid_seds[i,k].st_mass = mass[min_step_indx,min_mass_indx]
                grid_seds[i,k].st_bmass = bmass[min_step_indx,min_mass_indx]
                grid_seds[i,k].st_age = age[min_step_indx,min_mass_indx]
                grid_seds[i,k].st_tag = gtag[min_step_indx,min_mass_indx]

                ; save the luminosity of the source [solar luminosities]
                lum = 4.0*!PI*(grid_seds[i,k].st_radius*6.955d8*1e2)^2*t[*].(k+1)
                grid_seds[i,k].logl = grid_seds[i,k].st_logl
                grid_seds[i,k].total_flux = grid_seds[i,k].logl

                ; put the model at the distance specified and give it the isoc radius
                t[*].(k+1) *= ((grid_seds[i,k].st_radius*6.955d8*1e2)/(distance*3.09e16*1e2))^2
            
                ; save the stellar atmosphere numbers
                grid_seds[i,k].logt = alog10(temps[i])
                grid_seds[i,k].logg = gravs[k]
                grid_seds[i,k].waves = t[windxs].wavelength
                grid_seds[i,k].fluxes = t[windxs].(k+1)
            endif else begin
                grid_seds[i,k].total_flux = !values.f_nan
            endelse
        endif else begin
            grid_seds[i,k].total_flux = !values.f_nan
        endelse
    endfor
endfor

end
