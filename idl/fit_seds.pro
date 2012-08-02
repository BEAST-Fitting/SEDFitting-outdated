;+
; NAME:
;       FIT_SEDS
;
; PURPOSE:
;       Fit a model of a single star at a specified distance.  The
;       single star model includes both stellar physics and
;       interstellar extinction.
;
; CATEGORY:
;       Bayesian fitting.
;
; CALLING SEQUENCE:
;       FIT_SEDS
;
; INPUTS:
;       cat_fluxes : fluxes of stars [n_bands, n_stars]
;       cat_fluxes_unc : flux uncertainties of stars [n_bands, n_stars]
;       cat_coords : coordinates of stars [2, n_stars]
;       cat_mags : magnitudes of stars [n_bands, n_stars]
;       cat_mags_unc : magnitude uncertainties of stars [n_bands,
;                      n_stars]
;
; KEYWORD PARAMETERS:
;       modname : model name (created from 
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
; 	Started     : Karl Gordon (2010)
;       2010-2012   : lots of development (undocumented, KDG)
;       31 Jul 2012 : Cleaned up and full documentation added (KDG)
;-

pro fit_seds,cat_fluxes,cat_fluxes_unc,cat_coords,cat_mags,cat_mags_unc, $
             modname=modname, $
             use_save_seds=use_save_seds,max_chisqr=max_chisqr,noplot=noplot, $
             min_i=min_i,max_i=max_i,skip_i=skip_i,prompt=prompt,brick=brick, $
             out_path=out_path,distance=distance, $
             output_correlations=output_correlations,output_full_prob=output_full_prob, $
             filter_files=filter_files,filter_ascii=filter_ascii

n_cat = n_elements(cat_fluxes[0,*])

if (not keyword_set(brick)) then brick = 1
if (not keyword_set(distance)) then distance = 776e3 ; M31

if (not keyword_set(modname)) then modname = 'default'

brick_str = '_b' + strtrim(string(brick),2)

if (not keyword_set(out_path)) then outpathsmallbase = 'fit_seds' else $
  outpathsmallbase = out_path
outpathbase = outpathsmallbase + '_fitinfo'
out_path = outpathbase+'_'+modname+'/'

    ; filter filenames
if (not keyword_set(filter_files)) then begin
    filter_files = '~/Hubble/PHAT/SEDFitting/' + ['wfc3_uvis_f275w_002_syn.fits','wfc3_uvis_f336w_002_syn.fits', $
                                                  'acs_f475w_wfc_004_syn.fits','acs_f814w_wfc_005_syn.fits', + $
                                                  'wfc3_ir_f110w_002_syn.fits','wfc3_ir_f160w_003_syn.fits']
endif
 
; get the previously created band SEDs (created using make_model_seds)
restore,'fit_sed_band_seds_'+modname+'.sav'

; now do the fitting
if (not file_test(out_path,/directory)) then begin
    print,'making directory: ' + out_path
    file_mkdir,out_path
endif

ans = ''

if (not keyword_set(min_i)) then min_i = 0
if (not keyword_set(max_i)) then max_i = n_cat - 1
max_i = min([n_cat-1,max_i])

min_i = long(min_i)
save_min_i = min_i
max_i = long(max_i)

if (not keyword_set(skip_i)) then skip_i = 1

;min_i = long(float(min_i)/skip_i)
;max_i = long(float(max_i)/skip_i)
;offset_i = save_min_i mod skip_i

for z = min_i,max_i,skip_i do begin
;for i = min_i,max_i do begin
;    z = offset_i + i*skip_i

;    print,'fitting # = ' + strtrim(string(z+1),2) + ' out of ' + strtrim(string(max_i*skip_i+1),2)
    print,'fitting # = ' + strtrim(string(z+1),2) + ' out of ' + strtrim(string(max_i+1),2)

    find_best_fit_sed,cat_fluxes[*,z],cat_fluxes_unc[*,z],band_seds,best_fit_param,chisqr,prob, $
                      temp_av_image,temp_grav_image,grav_av_image,av_rv_image, $
                      av_prob,temp_prob,grav_prob,rv_prob, $
                      bmass_vals,bmass_prob,mass_vals,mass_prob,age_vals,age_prob, $
                      best_fit_fluxes,best_unred_fit_fluxes,best_fit_uncs, $
                      max_chisqr=max_chisqr,noplot=noplot

;    print,best_fit_param
;    stop

;    if (not keyword_set(noplot)) then begin
;        if (keyword_set(prompt)) then read,ans

;        window,0,xsize=600,ysize=600
;        kplot,10^band_seds.logt_vals,temp_prob,psym=100,kplot_type='oi',xtitle='temp'
;        window,2,xsize=600,ysize=600
;    kplot,band_seds.logg_vals,grav_prob,psym=100,xtitle='grav'
;        kplot,bmass_vals,bmass_prob,psym=100,xtitle='bmass',kplot_type='oi'
;        koplot,mass_vals,mass_prob,psym=100,linestyle=2
;        window,3,xsize=600,ysize=600
;        kplot,band_seds.av_vals,av_prob,psym=100,xtitle='av'
;    endif        

    ; output the images and vectors of probabilities

    out_filename = out_path+outpathsmallbase+brick_str+'_s'+strtrim(string(z+1),2)+'_'+modname+'.fits'
    fits_open,out_filename,ofcb,/write

    ; construct name
    rastr = deg2ra(cat_coords[0,z],sep_str='',/trunc_sec)
    decstr = deg2dec(cat_coords[1,z],sep_str='',/trunc_sec)
    name = 'CATB' + strtrim(string(brick),2) + '_' +  $
           rastr + decstr

    tarray = [[band_seds.resp_eff_wave],[cat_fluxes[*,z]],[cat_fluxes_unc[*,z]], $
              [cat_mags[*,z]],[cat_mags_unc[*,z]],[best_fit_fluxes],[best_unred_fit_fluxes]]
    fxhmake,header,tarray,/initialize
    sxaddpar,header,'NAME',name,'CAT unique name of this star (based on ra & dec)'
    sxaddpar,header,'BRICK',brick,'CAT brick number'
    sxaddpar,header,'SNUM',z+1,'star number'
    sxaddpar,header,'RA',cat_coords[0,z],'RA'
    sxaddpar,header,'DEC',cat_coords[1,z],'DEC'
    sxaddpar,header,'PROB',best_fit_param[6],'unnormalized probability of best fit (use cautiously)'
    sxaddpar,header,'CHISQR',best_fit_param[0],'chisqr of best fit'

    baynorm = product(1.0/sqrt(2.0*!PI*(cat_fluxes_unc[*,z]^2)))
    sxaddpar,header,'BAYNORM',baynorm,'Bayesian weight/normalization (based on observation uncertainties)'

    sxaddpar,header,'TEMP',best_fit_param[1],'log(Teff) of best fit'
    sxaddpar,header,'TEMPUNC',best_fit_uncs[0,0],'uncertainty in log(Teff) of best fit'
    sxaddpar,header,'TEMPMUNC',best_fit_uncs[0,1],'- uncertainty in log(Teff) of best fit'
    sxaddpar,header,'TEMPPUNC',best_fit_uncs[0,2],'+ uncertainty in log(Teff) of best fit'
    sxaddpar,header,'GRAV',best_fit_param[2],'gravity, log(g), of best fit'
    sxaddpar,header,'GRAVUNC',best_fit_uncs[1,0],'uncertainty in log(g) of best fit'
    sxaddpar,header,'GRAVMUNC',best_fit_uncs[1,1],'- uncertainty in log(g) of best fit'
    sxaddpar,header,'GRAVPUNC',best_fit_uncs[1,2],'+ uncertainty in log(g) of best fit'
    sxaddpar,header,'AV',best_fit_param[3],'A(V) of best fit'
    sxaddpar,header,'AVUNC',best_fit_uncs[2,0],'uncertainty in A(V) of best fit'
    sxaddpar,header,'AVMUNC',best_fit_uncs[2,1],'- uncertainty in A(V) of best fit'
    sxaddpar,header,'AVPUNC',best_fit_uncs[2,2],'+ uncertainty in A(V) of best fit'
    sxaddpar,header,'RV',best_fit_param[4],'R(V) of best fit'
    sxaddpar,header,'RVUNC',best_fit_uncs[3,0],'uncertainty in R(V) of best fit'
    sxaddpar,header,'RVMUNC',best_fit_uncs[3,1],'- uncertainty in R(V) of best fit'
    sxaddpar,header,'RVPUNC',best_fit_uncs[3,2],'+ uncertainty in R(V) of best fit'

    sxaddpar,header,'BMASS',best_fit_param[8],'birth mass of best fit'
    sxaddpar,header,'BMASSUNC',best_fit_uncs[4,0],'uncertainty in birth mass of best fit'
    sxaddpar,header,'BMASMUNC',best_fit_uncs[4,1],'- uncertainty in birth mass of best fit'
    sxaddpar,header,'BMASPUNC',best_fit_uncs[4,2],'+ uncertainty in birth mass of best fit'
    sxaddpar,header,'MASS',best_fit_param[9],'current mass of best fit'
    sxaddpar,header,'MASSUNC',best_fit_uncs[5,0],'uncertainty in mass of best fit'
    sxaddpar,header,'MASSMUNC',best_fit_uncs[5,1],'- uncertainty in mass of best fit'
    sxaddpar,header,'MASSPUNC',best_fit_uncs[5,2],'+ uncertainty in mass of best fit'
    sxaddpar,header,'AGE',best_fit_param[10],'current age of best fit'
    sxaddpar,header,'AGEUNC',best_fit_uncs[6,0],'uncertainty in age of best fit'
    sxaddpar,header,'AGEMUNC',best_fit_uncs[6,1],'- uncertainty in age of best fit'
    sxaddpar,header,'AGEPUNC',best_fit_uncs[6,2],'+ uncertainty in age of best fit'
    sxaddpar,header,'TLUM',best_fit_param[5],'best_fit luminosity'
    sxaddpar,header,'GNUM',best_fit_param[7],'number of good photometry points'
    sxaddpar,header,'COMMENT','columns are wavelength, flux, flux_unc, mag, mag_unc, model flux, unred model flux'
    sxaddpar,header,'EXTNAME','OBS_SED','Observed SED'
    fits_write,ofcb,tarray,header

    tarray = [[10^band_seds.logt_vals],[temp_prob]]
    fxhmake,header,tarray,/initialize
    sxaddpar,header,'EXTNAME','TEMP_PROB','probability name'
    fits_write,ofcb,tarray,header

    tarray = [[band_seds.logg_vals],[grav_prob]]
    fxhmake,header,tarray,/initialize
    sxaddpar,header,'EXTNAME','GRAV_PROB','probability name'
    fits_write,ofcb,tarray,header

    tarray = [[band_seds.av_vals],[av_prob]]
    fxhmake,header,tarray,/initialize
    sxaddpar,header,'EXTNAME','AV_PROB','probability name'
    fits_write,ofcb,tarray,header

    tarray = [[band_seds.rv_vals],[rv_prob]]
    fxhmake,header,tarray,/initialize
    sxaddpar,header,'EXTNAME','RV_PROB','probability name'
    fits_write,ofcb,tarray,header

    tarray = [[bmass_vals],[bmass_prob]]
    fxhmake,header,tarray,/initialize
    sxaddpar,header,'EXTNAME','BMASS_PROB','probability name'
    fits_write,ofcb,tarray,header

    tarray = [[mass_vals],[mass_prob]]
    fxhmake,header,tarray,/initialize
    sxaddpar,header,'EXTNAME','MASS_PROB','probability name'
    fits_write,ofcb,tarray,header

    tarray = [[age_vals],[age_prob]]
    fxhmake,header,tarray,/initialize
    sxaddpar,header,'EXTNAME','AGE_PROB','probability name'
    fits_write,ofcb,tarray,header

    if (keyword_set(output_correlations)) then begin
        fxhmake,header,temp_av_image,/initialize
        sxaddpar,header,'EXTNAME','TEMP_AV_IMAGE','image name'
        fits_write,ofcb,temp_av_image,header
        
        fxhmake,header,temp_grav_image,/initialize
        sxaddpar,header,'EXTNAME','TEMP_GRAV_IMAGE','image name'
        fits_write,ofcb,temp_grav_image,header

        fxhmake,header,grav_av_image,/initialize
        sxaddpar,header,'EXTNAME','GRAV_AV_IMAGE','image name'
        fits_write,ofcb,grav_av_image,header

        fxhmake,header,grav_av_image,/initialize
        sxaddpar,header,'EXTNAME','AV_RV_IMAGE','image name'
        fits_write,ofcb,av_rv_image,header
    endif

    if (keyword_set(output_full_prob)) then begin
        fxhmake,header,prob,/initialize
        sxaddpar,header,'EXTNAME','FULL_PROB','full 4D probability grid'
        fits_write,ofcb,prob,header
    endif

    fits_close,ofcb

;    fits_write,out_path+'obs_sed_chisqr'+strtrim(string(z+1),2)+'.fits',chisqr
;    fits_write,out_path+'obs_sed_prob'+strtrim(string(z+1),2)+'.fits',prob

    spawn,'gzip -f ' + out_filename + '&'

    if (keyword_set(prompt)) then read,ans

endfor

end
