; program to get the models fluxes in particular bands

pro get_sed_band_fluxes,grid_seds,filt_files,band_seds,ascii=ascii

n_filt = n_elements(filt_files)
n_wave = n_elements(grid_seds[0,0].waves)

resp_curves = fltarr(n_wave,n_filt)
resp_eff_wave = fltarr(n_filt)
tot_resp_curves = fltarr(n_filt)
for k = 0,(n_filt-1) do begin
    if (keyword_set(ascii)) then begin
        readcol,filt_files[k],wavelength,throughput,/silent
        linterp,wavelength,throughput,grid_seds[0,0,0].waves,tresp_curve,missing=0.0
    endif else begin
        t2 = mrdfits(filt_files[k],1,/silent)
        linterp,t2.wavelength,t2.throughput,grid_seds[0,0,0].waves,tresp_curve,missing=0.0
    endelse
    resp_curves[*,k] = tresp_curve
    tot_resp_curves[k] = total(resp_curves[*,k])
    resp_eff_wave[k] = (total(resp_curves[*,k]*grid_seds[0,0,0].waves)/tot_resp_curves[k])/1e4
endfor

size_grid_seds = size(grid_seds)
band_grid_seds = fltarr(size_grid_seds[1],size_grid_seds[2],size_grid_seds[3],n_filt)
grid_total_flux = fltarr(size_grid_seds[1],size_grid_seds[2],size_grid_seds[3])
grid_bmass = grid_total_flux
grid_mass = grid_total_flux
grid_age = grid_total_flux

logt_vals = fltarr(size_grid_seds[1])
logg_vals = fltarr(size_grid_seds[2])
av_vals = fltarr(size_grid_seds[3])
for i = 0,(size_grid_seds[1]-1) do begin
    logt_vals[i] = max(grid_seds[i,*,0].logt)
    for j = 0,(size_grid_seds[2]-1) do begin
        logg_vals[j] = max(grid_seds[*,j,0].logg)
        for k = 0,(size_grid_seds[3]-1) do begin
            av_vals[k] = grid_seds[i,j,k].av
;            grid_total_flux[i,j,k] = grid_seds[i,j,k].total_flux
            grid_total_flux[i,j,k] = 10^grid_seds[i,j,k].logl
            grid_bmass[i,j,k] = grid_seds[i,j,k].isoc_bmass
            grid_mass[i,j,k] = grid_seds[i,j,k].isoc_mass
            grid_age[i,j,k] = grid_seds[i,j,k].isoc_age
            for m = 0,(n_filt-1) do begin
                band_grid_seds[i,j,k,m] = total(grid_seds[i,j,k].fluxes*resp_curves[*,m])/tot_resp_curves[m]
            endfor
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
