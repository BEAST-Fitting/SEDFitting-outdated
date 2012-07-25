
pro find_best_fit_sed,obs_sed,obs_sed_unc,band_seds,best_fit_param,chisqr,prob, $
  temp_av_image,temp_grav_image,grav_av_image,av_rv_image, $
  av_prob,temp_prob,grav_prob,rv_prob, $
  bmass_vals,bmass_prob, $
  mass_vals,mass_prob, $
  age_vals,age_prob, $
  best_fit_fluxes,best_unred_fit_fluxes,best_fit_uncs, $
  max_chisqr=max_chisqr,noplot=noplot

if (not keyword_set(max_chisqr)) then max_chisqr = 100.0
if (not keyword_set(noplot)) then begin
    kplot,band_seds.resp_eff_wave,obs_sed,psym=1,symsize=3.0,kplot_type='oo'
endif

; now fit the bugger

size_band_seds = size(band_seds.band_grid_seds)
chisqr = replicate(0d0,size_band_seds[2],size_band_seds[3],size_band_seds[4],size_band_seds[5])

gindxs = where(obs_sed GT 0.0,n_gindxs)

for i = 0,(n_gindxs-1) do begin
    k = gindxs[i]
    chisqr += ((band_seds.band_grid_seds[k,*,*,*,*] - obs_sed[k])/obs_sed_unc[k])^2
endfor

; compute the probability
prob = exp(-1.0d0*chisqr/2.0)

; now take care of models with zero flux
indxs = where(band_seds.band_grid_seds[0,*,*,*,*] LE 0.0,n_indxs)
if (n_indxs GT 0) then begin
;    print,'n_models w/o flux = ', n_indxs, product([size_band_seds[2],size_band_seds[3],size_band_seds[4],size_band_seds[5]])
    prob[indxs] = 0.0
endif

; now find the best fit value
size_prob = size(prob)
max_prob = max(prob,indx)
l = indx/(size_prob[1]*size_prob[2]*size_prob[3])
min_ijk = indx mod (size_prob[1]*size_prob[2]*size_prob[3])
k = min_ijk/(size_prob[1]*size_prob[2])
min_ij = min_ijk mod (size_prob[1]*size_prob[2])
j = min_ij/size_prob[1]
i = min_ij mod size_prob[1]
if (max_prob NE prob[i,j,k,l]) then begin
    print,i,j,k,l,prob[i,j,k,l],max_prob
    stop
endif

best_fit_param = fltarr(12)
best_fit_param[0] = total(((band_seds.band_grid_seds[gindxs,i,j,k,l] - obs_sed[gindxs])/obs_sed_unc[gindxs])^2)
best_fit_param[1] = band_seds.logt_vals[i]
best_fit_param[2] = band_seds.logg_vals[j]
best_fit_param[3] = band_seds.av_vals[k]
best_fit_param[4] = band_seds.rv_vals[l]

best_fit_param[5] = band_seds.grid_total_flux[i,j,k,l]
best_fit_param[6] = prob[i,j,k,l]
best_fit_param[7] = n_gindxs
best_fit_param[8] = band_seds.grid_bmass[i,j,k,l]
best_fit_param[9] = band_seds.grid_mass[i,j,k,l]
best_fit_param[10] = band_seds.grid_age[i,j,k,l]

best_fit_fluxes = band_seds.band_grid_seds[*,i,j,k,l]
best_unred_fit_fluxes = band_seds.band_grid_seds[*,i,j,0,0]

; collapse the 4d cube to 4 images, one for each 2D representation of the
; cube - saves lots of space when writting to disk

prob_norv = total(prob,4)
temp_av_image = total(prob_norv,2)
temp_grav_image = total(prob_norv,3)
grav_av_image = total(prob_norv,1)

; get the av_rv image
prob_no_temp = total(prob,1)
av_rv_image = total(prob_no_temp,1)

; collapse the cube into 3 1D vectors
; one for each parameter

;k = n_elements(band_seds.logt_vals)
;n_avs = n_elements(band_seds.av_vals)
;weights = 10^band_seds.logt_vals[1:k-1] - 10^band_seds.logt_vals[0:k-2]
;weights = weights/total(weights)
;timage = 0.5*(temp_av_image[1:k-1,*] + temp_av_image[0:k-2,*])

;print,weights

;av_prob2 = total(timage*(weights#replicate(1.0,n_avs)),1)


;kplot,band_seds.av_vals,av_prob/max(av_prob),psym=100
;koplot,band_seds.av_vals,av_prob2/max(av_prob2),psym=1

;stop

av_prob = total(temp_av_image,1)
temp_prob = total(temp_av_image,2)
grav_prob = total(temp_grav_image,1)
rv_prob = total(av_rv_image,1)

; get the birth mass probability function
hmin = 0.5
hmax = 130.0
nbins = 50.
hmin_log10 = alog10(hmin)
hmax_log10 = alog10(hmax)
alog10_binsize = (hmax_log10 - hmin_log10)/nbins

histo_vals = hmin_log10 + alog10_binsize*(findgen(nbins)+0.5)
histo_vals = 10^histo_vals

indxs = where(prob GT 0.0,n_indxs)
bmass_vals = histo_vals
bmass_prob = dblarr(nbins)
if (n_indxs GT 0.0) then begin
    ; find all the models that have similar birth masses
    ; the wonder of reverse_indices
    histo = histogram(alog10(band_seds.grid_bmass[indxs]),min=hmin_log10,max=hmax_log10,nbins=nbins,reverse_indices=ri)

    for i = 0,(nbins-1) do begin
        if (ri[i] NE ri[i+1]) then begin
            bmass_prob[i] = total(prob[indxs[ri[ri[i]:ri[i+1]-1]]])
;            print,histo_vals[i],histo[i],min(band_seds.grid_bmass[indxs[ri[ri[i]:ri[i+1]-1]]]),max(band_seds.grid_bmass[indxs[ri[ri[i]:ri[i+1]-1]]]), $
;                  bmass_prob[i]
        endif
    endfor
endif; else begin
;    print,'prob is all zero; wtf?'
;    stop
;endelse

; get the current mass probability function
indxs = where(prob GT 0.0,n_indxs)
mass_vals = histo_vals
mass_prob = dblarr(nbins)
if (n_indxs GT 0.0) then begin
    ; find all the models that have similar birth masses
    ; the wonder of reverse_indices
    histo = histogram(alog10(band_seds.grid_mass[indxs]),min=hmin_log10,max=hmax_log10,nbins=nbins,reverse_indices=ri)

    for i = 0,(nbins-1) do begin
        if (ri[i] NE ri[i+1]) then begin
            mass_prob[i] = total(prob[indxs[ri[ri[i]:ri[i+1]-1]]])
;            print,histo_vals[i],histo[i],min(band_seds.grid_bmass[indxs[ri[ri[i]:ri[i+1]-1]]]),max(band_seds.grid_bmass[indxs[ri[ri[i]:ri[i+1]-1]]]), $
;                  bmass_prob[i]
        endif
    endfor

endif

; get the birth mass probability function
hmin = 1e3
hmax = 1e11
nbins = 100
hmin_log10 = alog10(hmin)
hmax_log10 = alog10(hmax)
alog10_binsize = (hmax_log10 - hmin_log10)/nbins

histo_vals = hmin_log10 + alog10_binsize*(findgen(nbins)+0.5)
histo_vals = 10^histo_vals

; get the current age probability function
indxs = where(prob GT 0.0,n_indxs)
age_vals = histo_vals
age_prob = dblarr(nbins)
if (n_indxs GT 0.0) then begin
    ; find all the models that have similar birth masses
    ; the wonder of reverse_indices
    histo = histogram(alog10(band_seds.grid_age[indxs]),min=hmin_log10,max=hmax_log10,nbins=nbins,reverse_indices=ri)

    for i = 0,(nbins-1) do begin
        if (ri[i] NE ri[i+1]) then begin
            age_prob[i] = total(prob[indxs[ri[ri[i]:ri[i+1]-1]]])
;            print,histo_vals[i],histo[i],min(band_seds.grid_bmass[indxs[ri[ri[i]:ri[i+1]-1]]]),max(band_seds.grid_bmass[indxs[ri[ri[i]:ri[i+1]-1]]]), $
;                  bmass_prob[i]
        endif
    endfor

endif

; now determine the 1-sigma uncertainties on the fit parameters using
; the PDFs

best_fit_uncs = fltarr(7,3)
best_fit_uncs[0,*] = get_uncs_from_pdf(band_seds.logt_vals,temp_prob,best_fit_param[1])
best_fit_uncs[1,*] = get_uncs_from_pdf(band_seds.logg_vals,grav_prob,best_fit_param[2])
best_fit_uncs[2,*] = get_uncs_from_pdf(band_seds.av_vals,av_prob,best_fit_param[3])
best_fit_uncs[3,*] = get_uncs_from_pdf(band_seds.rv_vals,rv_prob,best_fit_param[4])
best_fit_uncs[4,*] = get_uncs_from_pdf(bmass_vals,bmass_prob,best_fit_param[8])
best_fit_uncs[5,*] = get_uncs_from_pdf(mass_vals,mass_prob,best_fit_param[9])
best_fit_uncs[6,*] = get_uncs_from_pdf(age_vals,age_prob,best_fit_param[10])

end
