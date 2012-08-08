;+
; NAME:
;       GET_SED_GRID_FLUXES
;
; PURPOSE:
;       Determine the fluxes in specific bands for a grid of SEDs.
;
; CATEGORY:
;       Bayesian fitting.
;
; CALLING SEQUENCE:
;       GET_SED_GRID_FLUXES, grid_seds, filter_files, band_seds, filter_ascii=filter_ascii
;
; INPUTS:
;       grid_seds : 3D array of structures giving the model SED
;                   3D is logt, logg, A(V)
;       filter_files : filenames containing the filter response functions
;
; KEYWORD PARAMETERS:
;       filter_ascii : set to indicate the filter files are in ascii
;                      can be a vector of boleans for mixed format files
;
; OUTPUTS:
;      band_seds :
;
; OPTIONAL OUTPUTS:
;
; PROCEDURE:
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
; 	Started     : Karl Gordon (2010)
;       2010-2012   : development (undocumented, KDG)
;        7 Aug 2012 : Cleaned up and full documentation added (KDG)
;                     Code fixes to properly(?) compute band fluxes (integration including lambda term)
;-

pro get_sed_band_fluxes,grid_seds,filter_files,band_seds,filter_ascii=filter_ascii

n_filters = n_elements(filter_files)
n_waves = n_elements(grid_seds[0,0,0].waves)

if (not keyword_set(filter_ascii)) then filter_ascii = 0
if (n_elements(filter_ascii) EQ 1) then filter_ascii = replicate(filter_ascii,n_filters)
if (n_elements(filter_ascii) NE n_filters) then begin
    print,'number of filter files and ascii switches is not equal'
    print,'please fix and try again'
    return
endif

; read in the filter response curves
resp_curves = fltarr(n_waves,n_filters)
ave_resp_curves = fltarr(n_waves-1,n_filters)
resp_eff_wave = fltarr(n_filters)
tot_resp_curves = fltarr(n_filters)

; setup the needed arrays for the integration
ave_lambda = 0.5*(grid_seds[0,0,0].waves[1:n_waves-1] + grid_seds[0,0,0].waves[0:n_waves-2])
dlambda = grid_seds[0,0,0].waves[1:n_waves-1] - grid_seds[0,0,0].waves[0:n_waves-2]

; go for it
for k = 0,(n_filters-1) do begin
    if (keyword_set(filter_ascii[k])) then begin
        readcol,filter_files[k],wavelength,throughput,/silent
        print,k,max(wavelength)

;        if (max(wavelength) LT 5.) then stop

        if (max(wavelength) LT 1000.) then wavelength *= 1e4

        linterp,wavelength,throughput,grid_seds[0,0,0].waves,tresp_curve,missing=0.0

        
    endif else begin
        t2 = mrdfits(filter_files[k],1,/silent)

;        if (k EQ 0) then begin
;            indxs = where(t2.wavelength GT 3500.)
;            t2[indxs].throughput = 0.0
;        endif

        linterp,t2.wavelength,t2.throughput,grid_seds[0,0,0].waves,tresp_curve,missing=0.0
    endelse

    resp_curves[*,k] = tresp_curve

    ave_resp_curves[*,k] = 0.5*(resp_curves[1:n_waves-1,k] + resp_curves[0:n_waves-2,k])

    resp_eff_wave[k] = (total(ave_lambda*ave_resp_curves[*,k]*dlambda)/total(ave_resp_curves[*,k]*dlambda))/1e4

    ; multiply the tresp_curve by wavelength to provide the correct intergral in the big loop below
    tot_resp_curves[k] = total(ave_resp_curves[*,k]*dlambda)

endfor

; setup the output information
size_grid_seds = size(grid_seds)
band_grid_seds = fltarr(size_grid_seds[1],size_grid_seds[2],size_grid_seds[3],n_filters)
grid_total_flux = fltarr(size_grid_seds[1],size_grid_seds[2],size_grid_seds[3])
grid_bmass = grid_total_flux
grid_mass = grid_total_flux
grid_age = grid_total_flux

logt_vals = fltarr(size_grid_seds[1])
logg_vals = fltarr(size_grid_seds[2])
av_vals = fltarr(size_grid_seds[3])

ans = ''

; compute the fluxes in each band
for i = 0,(size_grid_seds[1]-1) do begin
    logt_vals[i] = max(grid_seds[i,*,0].logt)
    for j = 0,(size_grid_seds[2]-1) do begin
        logg_vals[j] = max(grid_seds[*,j,0].logg)
        for k = 0,(size_grid_seds[3]-1) do begin
            av_vals[k] = grid_seds[i,j,k].av
            if (grid_seds[i,j,k].logl GT 0.0) then begin
                grid_total_flux[i,j,k] = 10^grid_seds[i,j,k].logl
                grid_bmass[i,j,k] = grid_seds[i,j,k].st_bmass
                grid_mass[i,j,k] = grid_seds[i,j,k].st_mass
                grid_age[i,j,k] = grid_seds[i,j,k].st_age
                for m = 0,(n_filters-1) do begin ; see http://arxiv.org/abs/astro-ph/0210394
                    band_grid_seds[i,j,k,m] = total(0.5*(grid_seds[i,j,k].fluxes[1:n_waves-1] + grid_seds[i,j,k].fluxes[0:n_waves-2])*ave_resp_curves[*,m]*dlambda)/tot_resp_curves[m]
                endfor
                ; keep for red filter leak checks
                print,resp_eff_wave
                print,reform(band_grid_seds[i,j,k,*],n_filters), av_vals[k],logg_vals[j],logt_vals[i]
                kplot,resp_eff_wave,reform(band_grid_seds[i,j,k,*],n_filters),psym=1,kplot_type='oo'
                koplot,grid_seds[i,j,k].waves/1e4,grid_seds[i,j,k].fluxes,psym=100
                read,ans
            endif
        endfor
    endfor
endfor

band_seds = {band_grid_seds: band_grid_seds, $
             grid_total_flux: grid_total_flux, $
             grid_bmass: grid_bmass, $
             grid_mass: grid_mass, $
             grid_age: grid_age, $
             logt_vals: logt_vals, $
             logg_vals: logg_vals, $
             av_vals: av_vals, $
             resp_eff_wave: resp_eff_wave}

end
