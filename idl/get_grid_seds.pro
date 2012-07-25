; program to get the seds from a particular grid and
; put them at the distance of object of interest (e.g., m31)

pro get_grid_seds,distance,grid_seds,no_inter=no_inter

; get the isochronses first
get_isochrones,logl,logt,logg,radius,mass,bmass,gtag,age,no_inter=no_inter
size_logl = size(logl)

; read in the grid
files = file_search('~/Pro/Fit_SEDs/Stellar_Grids/ck04/p00/*_unitemps.moregravs.fits')
n_files = n_elements(files)
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

; remove temps over 37651 as I'm not seeing good matches with the
; stellar tracks - 24 Mar 2011
;indxs = where(temp_val LT 380000.,n_indxs)
;if (n_indxs GT 0) then begin
;    temp_val = temp_val[indxs]
;    files = files[indxs]
;    n_files = n_indxs
;endif

temps = temp_val
n_temps = n_elements(temps)
gravs = findgen(51)*0.1
n_gravs = n_elements(gravs)

; write the temp file
;openw,unit1,'~/Hubble/PHAT/SEDFitting/ck04/temps.dat',/get_lun
;for i = 0,(n_files-1) do begin
;    printf,unit1,temp_val[i]
;endfor
;free_lun,unit1

;n_files = 5
for i = 0,(n_files-1) do begin
    t = mrdfits(files[i],1,/silent)
    windxs = where((t.wavelength GT 912.) and (t.wavelength LT 30000.),n_twave)
;    windxs = where((t.wavelength NE 0.0),n_twave)
;    n_twave = n_elements(t.wavelength)
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
               isoc_logt: 0.0, isoc_logg : 0.0, isoc_logl: 0D0, isoc_radius: 0.0, $
               isoc_mass: 0.0, isoc_bmass: 0.0, isoc_age: 0.0, isoc_tag: 0, $
               waves: fltarr(n_wave), $
               fluxes: dblarr(n_wave), $
               total_flux: 0D0}
        grid_seds = replicate(sed,n_files,n_gravs)

    endif
    if (n_twave NE n_wave) then begin
        print,'current number of wavelengths different from 1st one'
        print,'current # waves = ', n_twave
        print,'    1st # waves = ', n_wave
        stop
    endif
    if (n_tgravs NE n_gravs) then begin
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

;            print,'***'
;            print,max(logg),max(logt)
;            print,gravs[k],alog10(temps[i])
;            print,logg[min_step_indx,min_mass_indx],logt[min_step_indx,min_mass_indx]

            if ((abs(logg[min_step_indx,min_mass_indx] - gravs[k]) LT 0.02) AND $
                (abs(logt[min_step_indx,min_mass_indx] - alog10(temps[i])) LT 0.02)) then begin
;            if ((abs(logg[min_step_indx,min_mass_indx] - gravs[k]) LT 1000.) AND $
;                (abs(logt[min_step_indx,min_mass_indx] - alog10(temps[i])) LT 1000.)) then begin

;                print,'got one'

                grid_seds[i,k].isoc_logt = logt[min_step_indx,min_mass_indx]
                grid_seds[i,k].isoc_logg = logg[min_step_indx,min_mass_indx]
                grid_seds[i,k].isoc_logl = logl[min_step_indx,min_mass_indx]
                grid_seds[i,k].isoc_radius = radius[min_step_indx,min_mass_indx]
                grid_seds[i,k].isoc_mass = mass[min_step_indx,min_mass_indx]
                grid_seds[i,k].isoc_bmass = bmass[min_step_indx,min_mass_indx]
                grid_seds[i,k].isoc_age = age[min_step_indx,min_mass_indx]
                grid_seds[i,k].isoc_tag = gtag[min_step_indx,min_mass_indx]

            ; save the luminosity of the source [solar luminosities]
                lum = 4.0*!PI*(grid_seds[i,k].isoc_radius*6.955d8*1e2)^2*t[*].(k+1)
;                grid_seds[i,k].total_flux = int_tabulated(t[*].wavelength,lum)/3.839d33
;            grid_seds[i,k].logl = alog10(grid_seds[i,k].total_flux)
                grid_seds[i,k].logl = grid_seds[i,k].isoc_logl
                grid_seds[i,k].total_flux = grid_seds[i,k].logl

            ; put the model at the distance specified and give it the isoc radius
                t[*].(k+1) *= ((grid_seds[i,k].isoc_radius*6.955d8*1e2)/(distance*3.09e16*1e2))^2
            
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
