; program to extinguish SEDs

pro extinguish_seds,grid_seds,ext_grid_seds,rv=rv,max_av=max_av,min_av=min_av,av_step=av_step

if (not keyword_set(rv)) then rv = 3.1

if (not keyword_set(max_av)) then max_av = 1.0
if (not keyword_set(min_av)) then min_av = 0.0
if (not keyword_set(av_step)) then av_step = 0.025
n_av = round((max_av - min_av)/av_step) + 1
many_av = min_av + findgen(n_av)*av_step

; assuming all the seds have the same wavelength scale
ext_curve = ccm(rv,1.e4/grid_seds[0,0].waves)

size_grid_seds = size(grid_seds)
n_wave = n_elements(grid_seds[0,0].waves)
red_sed = {logt: 0.0, logg : 0.0, logl: 0D0, av: 0.0, rv: 0.0, $
           isoc_logt: 0.0, isoc_logg : 0.0, isoc_logl: 0D0, isoc_radius: 0.0, $
           isoc_mass: 0.0, isoc_bmass: 0.0, isoc_age: 0.0, isoc_tag: 0, $
           waves: fltarr(n_wave), $
           fluxes: dblarr(n_wave), $
           total_flux: 0D0}
ext_grid_seds = replicate(red_sed,size_grid_seds[1],size_grid_seds[2],n_av)

for i = 0,(n_av-1) do begin
    ext_grid_seds[*,*,i] = grid_seds

;    if (i EQ 0) then begin
;        kplot,ext_grid_seds[10,10,0].waves,ext_grid_seds[10,10,0].fluxes,psym=100,kplot_type='oo'
;    endif

    ext_grid_seds[*,*,i].av = many_av[i]
    ext_grid_seds[*,*,i].rv = rv
    for j = 0,(size_grid_seds[1]-1) do begin
        for k = 0,(size_grid_seds[2]-1) do begin
            ext_grid_seds[j,k,i].fluxes *= 10^(-0.4*ext_curve*many_av[i])
        endfor
    endfor

;    koplot,ext_grid_seds[10,10,i].waves,ext_grid_seds[10,10,i].fluxes,psym=100

endfor

end
