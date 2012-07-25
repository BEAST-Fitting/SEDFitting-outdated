
pro make_model_seds, modname=modname, smodname=smodname, $
                     distance=distance, $
                     filter_files=filter_files, $
                     av_step=av_step, av_range=av_range, $
                     rv_step=rv_step, rv_range=rv_range, $
                     use_save_spectra=use_save_spectra

if (not keyword_set(modname)) then modname = 'default'
if (not keyword_set(smodname)) then smodname = modname
if (not keyword_set(distance)) then distance = 776e3 ; M31

if (not keyword_set(av_step)) then av_step = 0.1
if (not keyword_set(av_range)) then av_range = [0.0,10.0]

if (not keyword_set(rv_step)) then rv_step = 0.5
if (not keyword_set(rv_range)) then rv_range = [1.5,6.0]

    ; filter filenames
if (not keyword_set(filter_files)) then begin
    filter_files = '~/Pro/Fit_SEDs/Filters/' + ['wfc3_uvis_f275w_002_syn.fits','wfc3_uvis_f336w_002_syn.fits', $
                                                'acs_f475w_wfc_004_syn.fits','acs_f814w_wfc_005_syn.fits', + $
                                                'wfc3_ir_f110w_002_syn.fits','wfc3_ir_f160w_003_syn.fits']
endif
n_filters = n_elements(filter_files)

; get the stellar atmosphere models, put them at the distance
; specified, and the luminosities (via radii) from stellar
; evolutionary tracks (isochrones)

print,'getting the model SEDs (atmosphere and evolutionary tracks)...'
if (not keyword_set(use_save_spectra)) then begin
    get_grid_seds,distance,grid_seds
    save,grid_seds,filename='fit_sed_grid_seds_'+smodname+'.sav'
endif else begin
    restore,'fit_sed_grid_seds_'+smodname+'.sav'
endelse

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
    print,'working on R(V) = ' + strtrim(string(rv,format='(F6.2)'),2)
    
    ; extinguish the SEDs for a range of A(V) values using a standard
    ; extinction curve
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

        print,'   working on A(V) = ' + strtrim(string(min_chunk_av+av_step,format='(F6.2)'),2) + ' to ' + $
              strtrim(string(max_chunk_av,format='(F6.2)'),2) + ' in steps of ' + $
              strtrim(string(av_step,format='(F6.2)'),2)
;        print,j,n_av_for_chunk,startk,endk

        extinguish_seds,grid_seds,ext_grid_seds,rv=rv,min_av=min_chunk_av+av_step,max_av=max_chunk_av,av_step=av_step
        get_sed_band_fluxes,ext_grid_seds,filter_files,band_seds,ascii=filter_ascii

        band_grid_seds[*,*,startk:endk,i,*] = band_seds.band_grid_seds
        grid_total_flux[*,*,startk:endk,i] = band_seds.grid_total_flux
        grid_bmass[*,*,startk:endk,i] = band_seds.grid_bmass
        grid_mass[*,*,startk:endk,i] = band_seds.grid_mass
        grid_age[*,*,startk:endk,i] = band_seds.grid_age
        
        startk = endk + 1
        
        ext_grid_seds = 0.

    endfor

endfor

band_seds = {band_grid_seds: band_grid_seds, $
             grid_total_flux: grid_total_flux, $
             grid_bmass: grid_bmass, $
             grid_mass: grid_mass, $
             grid_age: grid_age, $
             logt_vals: band_seds.logt_vals, $
             logg_vals: band_seds.logg_vals, $
             av_vals: many_av, $
             rv_vals: many_rv, $
             resp_eff_wave: band_seds.resp_eff_wave}

; save the result
save,band_seds,filename='fit_sed_band_seds_'+modname+'.sav'

end
