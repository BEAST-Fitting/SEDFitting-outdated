;+
; NAME:
;       MAKE_MODEL_SEDS
;
; PURPOSE:
;       Create an IDL save file with a grid of models of a single star
;       in a galaxy at a known distance.  These single star models
;       include both stellar physics and dust extinction.
;
; CATEGORY:
;       Bayesian fitting.
;
; CALLING SEQUENCE:
;       MAKE_MODEL_SEDS
;
; INPUTS:
;
; KEYWORD PARAMETERS:
;       distance : distance to galaxy in pc [default is M31, 776e3 pc]
;       filter_names : names of the filters
;       filter_files : filenames containing the response functions for
;                      each filter, assumed to be FITS tables
;       filter_ascii : set this keyword if the filter response files
;                      are in ASCII with 2 columns
;                           wavelength [Angstrom] & response
;       av_step : step size in A(V) [default 0.1]
;       av_range : range in A(V) [default 0.0,10.0]
;       rv_step : step size in R(V) [default 0.5
;       rv_range : range in R(V) [default 1.5,6.0]
;
;       modname : name of the saved model file
;                  fit_sed_band_seds_[modname].sav
;       smodname : name of the stellar spectra save file [default smodname=modname]
;                  fit_sed_grid_seds_[smodname].sav
;       use_save_spectra : use fit_sed_grid_seds_[smodname].sav instead of creating the 
;                          stellar models SEDs from stellar atmosphere and evolutionary track models
;       silent : set to turn of screen output
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
; 	Started     : Karl Gordon (2011)
;       2011-2012   : lots of development (undocumented, KDG)
;       31 Jul 2012 : Cleaned up and full documentation added (KDG)
;        7 Aug 2012 : added filter names to the input and output (KDG)
;                     also added filter_files to the output
;-

pro make_model_seds, distance=distance, $
                     filter_names=filter_names,filter_files=filter_files, $
                     filter_ascii=filter_ascii, $
                     av_step=av_step, av_range=av_range, $
                     rv_step=rv_step, rv_range=rv_range, $
                     modname=modname, smodname=smodname, $
                     use_save_spectra=use_save_spectra, $
                     silent=silent

if (not keyword_set(modname)) then modname = 'default'
if (not keyword_set(smodname)) then smodname = modname
if (not keyword_set(distance)) then distance = 776e3 ; M31

if (not keyword_set(av_step)) then av_step = 0.1
if (not keyword_set(av_range)) then av_range = [0.0,10.0]

if (not keyword_set(rv_step)) then rv_step = 0.5
if (not keyword_set(rv_range)) then rv_range = [1.5,6.0]

    ; filter filenames
if (not keyword_set(filter_files)) then begin
    ; defaults for PHAT M31 MCT
    filter_files = '~/Pro/Fit_SEDs/Filters/' + ['wfc3_uvis_f275w_002_syn.fits','wfc3_uvis_f336w_002_syn.fits', $
                                                'acs_f475w_wfc_004_syn.fits','acs_f814w_wfc_005_syn.fits', + $
                                                'wfc3_ir_f110w_002_syn.fits','wfc3_ir_f160w_003_syn.fits']
endif
n_filters = n_elements(filter_files)

if (not keyword_set(filter_names)) then filter_names = repliate('unknown',n_filters)

if (n_elements(filter_names) NE n_filters) then begin
    print,'error: the number of filter_names does not match the number of filter_files'
    print,'correct and try again'
    return
end

; get the stellar atmosphere models, put them at the distance
; specified, and the luminosities (via radii) from stellar
; evolutionary tracks

if (not keyword_set(silent)) then print,'getting the model SEDs (atmosphere and evolutionary tracks)...'
if (not keyword_set(use_save_spectra)) then begin
    get_grid_stellar_seds,distance,grid_seds,silent=silent
    save,grid_seds,filename='fit_sed_grid_seds_'+smodname+'.sav'
endif else begin
    restore,'fit_sed_grid_seds_'+smodname+'.sav'
endelse

; determine the A(V) and R(V) values for the model grid
n_rv = round((rv_range[1] - rv_range[0])/rv_step) + 1
many_rv = rv_range[0] + findgen(n_rv)*rv_step

n_av = round(round(av_range[1] - av_range[0])/av_step) + 1
many_av = av_range[0] + findgen(n_av)*av_step

; setup the final grid variables
size_grid_seds = size(grid_seds)
n_logt_vals = size_grid_seds[1]
n_logg_vals = size_grid_seds[2]

band_grid_seds = fltarr(n_logt_vals,n_logg_vals,n_av,n_rv,n_filters)
grid_total_flux = fltarr(n_logt_vals,n_logg_vals,n_av,n_rv)
grid_bmass = grid_total_flux
grid_mass = grid_total_flux
grid_age = grid_total_flux

; loop over the R(V) values
for i = 0,(n_rv-1) do begin
    rv = many_rv[i]
    if (not keyword_set(silent)) then print,'working on R(V) = ' + strtrim(string(rv,format='(F6.2)'),2)
    
    ; extinguish the SEDs for a range of A(V) values extinction curve
    ; need to do this in chunks to avoid memory problems

    n_av_chunk = round(n_av/40) + 1
    av_chunk_delta = fix(n_av/n_av_chunk)*av_step

    max_chunk_av = av_range[0] - av_step ; setup the start of the 1st chunk
    startk = 0
    for j = 0,(n_av_chunk-1) do begin

        min_chunk_av = max_chunk_av
        max_chunk_av = min_chunk_av + av_chunk_delta
        if (j EQ (n_av_chunk-1)) then max_chunk_av = av_range[1]

        n_av_for_chunk = round((max_chunk_av - (min_chunk_av+av_step))/av_step) + 1
        endk = startk + n_av_for_chunk - 1

        if (not keyword_set(silent)) then $
          print,'   working on A(V) = ' + strtrim(string(min_chunk_av+av_step,format='(F6.2)'),2) + ' to ' + $
          strtrim(string(max_chunk_av,format='(F6.2)'),2) + ' in steps of ' + $
          strtrim(string(av_step,format='(F6.2)'),2)

        ; extinguish the stellar SEDs
        extinguish_seds,grid_seds,ext_grid_seds,rv=rv,min_av=min_chunk_av+av_step,max_av=max_chunk_av,av_step=av_step
        ; multiply the extinguished SEDs by the band response to get band fluxes
        get_sed_band_fluxes,ext_grid_seds,filter_files,band_seds,filter_ascii=filter_ascii

        band_grid_seds[*,*,startk:endk,i,*] = band_seds.band_grid_seds
        grid_total_flux[*,*,startk:endk,i] = band_seds.grid_total_flux
        grid_bmass[*,*,startk:endk,i] = band_seds.grid_bmass
        grid_mass[*,*,startk:endk,i] = band_seds.grid_mass
        grid_age[*,*,startk:endk,i] = band_seds.grid_age
        
        startk = endk + 1
        
        ext_grid_seds = 0.

    endfor

endfor

; setup the final structure to save
band_seds = {band_grid_seds: band_grid_seds, $
             grid_total_flux: grid_total_flux, $
             grid_bmass: grid_bmass, $
             grid_mass: grid_mass, $
             grid_age: grid_age, $
             logt_vals: band_seds.logt_vals, $
             logg_vals: band_seds.logg_vals, $
             av_vals: many_av, $
             rv_vals: many_rv, $
             resp_eff_wave: band_seds.resp_eff_wave, $
             filter_names: filter_names, $
             filter_files: filter_files, $
             distance: distance}

; save the result
save,band_seds,filename='fit_sed_band_seds_'+modname+'.sav'

end
