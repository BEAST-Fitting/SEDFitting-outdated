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
;       FIT_SEDS,cat_fluxes, cat_fluxes_unc, cat_coords, cat_mags, cag_mags_unc
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
;       modfile : model filename (created from make_model_seds.pro)
;                 necessary for code to run
;       outpath : name of dir in which to output resulting FITS files
;       outbase : output base of filenames
;       output_1d_likelihood : set to output the 1D marginalized likelihoods
;       output_full_likelihood : set to output the full nD likelihoods
;       cat_prefix : catalog name for star names (default = "PHAT")
;       silent : set to turn of screen output
;
;       **parameters to control the number of stars fit**
;       min_i : star number to start fitting
;       max_i : star number to stop fitting
;       skip_i : number of stars to skip between fitting
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
;        2 Aug 2012 : more cleanup (KDG)
;        6 Aug 2012 : changed output to batch output
;                        including binary table of fit parameters (used to be done later)
;        7 Aug 2012 : testing and bugfixes
;-

pro fit_seds,cat_fluxes,cat_fluxes_unc,cat_coords,cat_mags,cat_mags_unc, $
             modfile=modfile,out_path=outpath,outbase=outbase, $
             cat_prefix=cat_prefix,n_out_chunk=n_out_chunk, $
             output_1d_likelihood=output_1d_likelihood, $
             output_full_likelihood=output_full_likelihood, $
             min_i=min_i,max_i=max_i,skip_i=skip_i, $
             silent=silent

; defaults
if (not keyword_set(modfile)) then modfile = 'fits_seds_band_seds_default.sav'
if (not keyword_set(outpath)) then outpath = 'fit_seds_default'
outpath += '/' ; make sure this is a directory
if (not keyword_set(outbase)) then outbase = 'default'
if (not keyword_set(cat_prefix)) then cat_prefix = 'PHAT'
if (not keyword_set(n_out_chunk)) then n_out_chunk = 100

; number of star in catalog
n_cat = n_elements(cat_fluxes[0,*])

if (not keyword_set(min_i)) then min_i = 0
if (not keyword_set(max_i)) then max_i = n_cat - 1
if (not keyword_set(skip_i)) then skip_i = 1
max_i = min([n_cat-1,max_i])
min_i = long(min_i)
save_min_i = min_i
max_i = long(max_i)

; get the previously created band SEDs (created using make_model_seds)
if (file_test(modfile)) then begin
    restore,modfile
endif else begin
    print,'no model IDL save file found'
    print,'modfile = ' + modfile
    return
endelse

; setup the output file
if (not file_test(outpath,/directory)) then begin
    if (not keyword_set(silent)) then print,'making directory: ' + outpath
    file_mkdir,outpath
endif

; size of model grid
size_band_seds = size(band_seds.band_grid_seds)
n_waves = size_band_seds[5]
if (not keyword_set(silent)) then print,'n_bands  ',n_waves

; setup the indexes where the model does not exisit
good_model_indxs = where(band_seds.band_grid_seds[*,*,*,*,0] GT 0.0,n_good_model_indxs)
no_model_indxs = where(band_seds.band_grid_seds[*,*,*,*,0] LE 0.0,n_no_model_indxs)

; setup the likelihood grid - log(likelihood) - only float needed
log_likelihood = replicate(0.0,size_band_seds[1],size_band_seds[2],size_band_seds[3],size_band_seds[4])
log_likelihood[no_model_indxs] = -75.  ; small enough?  (seems safe for float, can go lower if we use double)
if (not keyword_set(silent)) then print,'# of good model points, model points',n_good_model_indxs, product(size_band_seds[1:4])

; ***setup the reverse histogram indices for determining "free" model parameters
; birth and current mass bins setup
mass_hmin = 0.5
mass_hmax = 130.0
mass_nbins = 50.
mass_hmin_log10 = alog10(mass_hmin)
mass_hmax_log10 = alog10(mass_hmax)
mass_alog10_binsize = (mass_hmax_log10 - mass_hmin_log10)/mass_nbins

mass_vals = mass_hmin_log10 + mass_alog10_binsize*(findgen(mass_nbins)+0.5)
mass_vals = 10^mass_vals

; the wonder of reverse_indices
; find all the models that have similar birth masses
histo = histogram(alog10(band_seds.grid_bmass[good_model_indxs]),min=mass_hmin_log10,max=mass_hmax_log10, $
                  nbins=mass_nbins,reverse_indices=bmass_ri)

; find all the models that have similar birth masses
histo = histogram(alog10(band_seds.grid_mass[good_model_indxs]),min=mass_hmin_log10,max=mass_hmax_log10, $
                  nbins=mass_nbins,reverse_indices=mass_ri)

; age bins setup
age_hmin = 1e3
age_hmax = 1e11
age_nbins = 100
age_hmin_log10 = alog10(age_hmin)
age_hmax_log10 = alog10(age_hmax)
age_alog10_binsize = (age_hmax_log10 - age_hmin_log10)/age_nbins

age_vals = age_hmin_log10 + age_alog10_binsize*(findgen(age_nbins)+0.5)
age_vals = age_vals

histo = histogram(alog10(band_seds.grid_age[good_model_indxs]),min=age_hmin_log10,max=age_hmax_log10, $
                  nbins=age_nbins,reverse_indices=age_ri)

bmass_likelihood = fltarr(mass_nbins)
mass_likelihood = fltarr(mass_nbins)
age_likelihood = fltarr(age_nbins)

; setup the output structure (for the FITS binary table)
a = {name: "", number: 0L, ra: double(0.0), dec: 0D0, max_likelihood: 0.0, $
     logt: 0.0, logt_unc: 0.0, logg: 0.0, logg_unc: 0.0, av: 0.0, av_unc: 0.0, rv: 0.0, rv_unc: 0.0, $
     bmass: 0.0, bmass_unc: 0.0, mass: 0.0, mass_unc: 0.0, log_age: 0.0, log_age_unc: 0.0, $
     luminosity: 0.0, radius: 0.0, $
     mags: fltarr(n_waves), mags_unc: fltarr(n_waves), $
;     mod_fluxes: fltarr(n_waves), unred_mod_fluxes: fltarr(n_waves), $
     n_gindxs: 0}
out_table = replicate(a,n_out_chunk)

; setup the output arrays for the 1d likelihoods
;  1st (bottom) row gives the x values
if (keyword_set(output_1d_likelihood)) then begin
    
    out_1d_logt_likelihood = fltarr(size_band_seds[1],n_out_chunk+1)
    out_1d_logt_likelihood[*,0] = band_seds.logt_vals

    out_1d_logg_likelihood = fltarr(size_band_seds[2],n_out_chunk+1)
    out_1d_logg_likelihood[*,0] = band_seds.logg_vals

    out_1d_av_likelihood = fltarr(size_band_seds[3],n_out_chunk+1)
    out_1d_av_likelihood[*,0] = band_seds.av_vals

    out_1d_rv_likelihood = fltarr(size_band_seds[4],n_out_chunk+1)
    out_1d_rv_likelihood[*,0] = band_seds.rv_vals

    out_1d_bmass_likelihood = fltarr(mass_nbins,n_out_chunk+1)
    out_1d_bmass_likelihood[*,0] = mass_vals

    out_1d_mass_likelihood = fltarr(mass_nbins,n_out_chunk+1)
    out_1d_mass_likelihood[*,0] = mass_vals

    out_1d_age_likelihood = fltarr(age_nbins,n_out_chunk+1)
    out_1d_age_likelihood[*,0] = age_vals
endif

cur_tot_out = n_out_chunk + 1
cur_chunk = -1
for z = min_i,max_i,skip_i do begin

    if (not keyword_set(silent)) then $
      print,'fitting # = ' + strtrim(string(z+1),2) + ' out of ' + strtrim(string(max_i+1),2)

    ; need to update this when upper limits are available
    obs_sed = cat_fluxes[*,z]
    obs_sed_unc = cat_fluxes_unc[*,z]

    ; loop over the number of bands with good fluxes
    ;  only compute for the good fluxes *and* good models
    gindxs = where(obs_sed GT 0.0,n_gindxs)
    k = gindxs[0]
    ; initialize
    log_likelihood[good_model_indxs] = (((band_seds.band_grid_seds[*,*,*,*,k])[good_model_indxs] - obs_sed[k])/obs_sed_unc[k])^2
    for i = 1,(n_gindxs-1) do begin
        k = gindxs[i]
        log_likelihood[good_model_indxs] += (((band_seds.band_grid_seds[*,*,*,*,k])[good_model_indxs] - obs_sed[k])/obs_sed_unc[k])^2
    endfor
    ; compute the log_likelihood from the chisqr
    log_likelihood[good_model_indxs] = -0.5*log_likelihood[good_model_indxs]

    ; collapse the 4D grids
    ;   there is probably a better way to do this, maybe with histograms and indices (precompute?)
    max_likelihood = max(log_likelihood)
    likelihood_norv = total(exp(log_likelihood - max_likelihood),4)
    logt_av_image = total(likelihood_norv,2)
    logt_logg_image = total(likelihood_norv,3)
    logg_av_image = total(likelihood_norv,1)
    ; now get the av_rv image
    likelihood_nologt = total(exp(log_likelihood - max_likelihood),1)
    av_rv_image = total(likelihood_nologt,1)

    ; finally, get the 1D likelihoods, marginalized over all other parameters
    logt_likelihood = alog(total(logt_av_image,2)) + max_likelihood
    logg_likelihood = alog(total(logt_logg_image,1)) + max_likelihood
    av_likelihood = alog(total(logt_av_image,1)) + max_likelihood
    rv_likelihood = alog(total(av_rv_image,1)) + max_likelihood

    ; get the 1D likelihoods for the "free" model parameters
    ; birth mass
    for i = 0,(mass_nbins-1) do begin
        if (bmass_ri[i] NE bmass_ri[i+1]) then begin
            bmass_likelihood[i] = alog(total(exp(log_likelihood[good_model_indxs[bmass_ri[bmass_ri[i]:bmass_ri[i+1]-1]]]) $
                                             - max_likelihood)) + max_likelihood
        endif else bmass_likelihood[i] = -75.
    endfor
    ; current mass
    for i = 0,(mass_nbins-1) do begin
        if (mass_ri[i] NE mass_ri[i+1]) then begin
            mass_likelihood[i] = alog(total(exp(log_likelihood[good_model_indxs[mass_ri[mass_ri[i]:mass_ri[i+1]-1]]]) $
                                             - max_likelihood)) + max_likelihood
        endif else mass_likelihood[i] = -75.
    endfor
    ; age
    for i = 0,(age_nbins-1) do begin
        if (age_ri[i] NE age_ri[i+1]) then begin
            age_likelihood[i] = alog(total(exp(log_likelihood[good_model_indxs[age_ri[age_ri[i]:age_ri[i+1]-1]]]) $
                                             - max_likelihood)) + max_likelihood
        endif else age_likelihood[i] = -75.
    endfor

    ; extract the most probable values from the 1D likelihoods and their 1-sigma uncertainties
    mp_logt = get_most_probable(band_seds.logt_vals,logt_likelihood)
    mp_logg = get_most_probable(band_seds.logg_vals,logg_likelihood)
    mp_av = get_most_probable(band_seds.av_vals,av_likelihood)
    mp_rv = get_most_probable(band_seds.rv_vals,rv_likelihood)

    mp_bmass = get_most_probable(mass_vals,bmass_likelihood)
    mp_mass = get_most_probable(mass_vals,mass_likelihood)
    mp_age = get_most_probable(age_vals,age_likelihood)

    ; setup the output
    ;   created so that batch output is the norm - too many files otherwise 

    if ((cur_tot_out+1) GE n_out_chunk) then begin
        
        cur_tot_out = -1 ; reset
        cur_chunk++ ; new output chunk

        ; filename for output
        out_filename = outpath+outbase+'_chunk_'+strtrim(string(cur_chunk+1),2)+'.fits'

        if (keyword_set(output_full_likelihood)) then begin
            fits_open,repstr(out_filename,'.fits','_full_likelihood.fits'),ofcb_full,/write
        endif
        
    endif

    cur_tot_out++  ; add one to keep track of how many have been output in the current file(s)

    ; construct name
    rastr = fit_seds_deg2ra(cat_coords[0,z],sep_str='',/trunc_sec)
    decstr = fit_seds_deg2dec(cat_coords[1,z],sep_str='',/trunc_sec)
    name = cat_prefix + '_' + rastr + decstr

    out_table[cur_tot_out].name = name
    out_table[cur_tot_out].number = z + 1
    out_table[cur_tot_out].ra = cat_coords[0,z]
    out_table[cur_tot_out].dec = cat_coords[1,z]
    out_table[cur_tot_out].max_likelihood = max_likelihood

    ; model parameters
    out_table[cur_tot_out].logt = mp_logt[0]
    out_table[cur_tot_out].logt_unc = mp_logt[1]
    out_table[cur_tot_out].logg = mp_logg[0]
    out_table[cur_tot_out].logg_unc = mp_logg[1]
    out_table[cur_tot_out].av = mp_av[0]
    out_table[cur_tot_out].av_unc = mp_av[1]
    out_table[cur_tot_out].rv = mp_rv[0]
    out_table[cur_tot_out].rv_unc = mp_rv[1]

    ; other model parameters (for free)
    out_table[cur_tot_out].bmass = mp_bmass[0]
    out_table[cur_tot_out].bmass_unc = mp_bmass[1]
    out_table[cur_tot_out].mass = mp_mass[0]
    out_table[cur_tot_out].mass_unc = mp_mass[1]
    out_table[cur_tot_out].log_age = mp_age[0]
    out_table[cur_tot_out].log_age_unc = mp_age[1]

    out_table[cur_tot_out].mags = cat_mags[*,z]
    out_table[cur_tot_out].mags_unc = cat_mags_unc[*,z]
;    out_table[cur_tot_out].mod_fluxes = 
;    out_table[cur_tot_out].unred_mod_fluxes = 
    out_table[cur_tot_out].n_gindxs = n_gindxs

    ; compute the radius
    radius2 = (out_table[cur_tot_out].luminosity)*3.839d26/(4.0*!PI*5.67d-8*((10^out_table[cur_tot_out].logt)^4))
    out_table[cur_tot_out].radius = sqrt(radius2)/6.955d8

    if (keyword_set(output_1d_likelihood)) then begin
        out_1d_logt_likelihood[*,cur_tot_out+1] = logt_likelihood
        out_1d_logg_likelihood[*,cur_tot_out+1] = logg_likelihood
        out_1d_av_likelihood[*,cur_tot_out+1] = av_likelihood
        out_1d_rv_likelihood[*,cur_tot_out+1] = rv_likelihood
        out_1d_bmass_likelihood[*,cur_tot_out+1] = bmass_likelihood
        out_1d_mass_likelihood[*,cur_tot_out+1] = mass_likelihood
        out_1d_age_likelihood[*,cur_tot_out+1] = age_likelihood
    endif

    if (keyword_set(output_full_likelihood)) then begin
        ; each full likelihood gets a new extension (allows for varying sparse matrices)
        out_indxs = where((log_likelihood - max_likelihood) GT -40.0,n_out_indxs)
        if (n_out_indxs GT 0) then begin

;            out_full[0:n_out_indxs-1].indexs = out_indxs
;            out_full[0:n_out_indxs-1].log_likelihood = log_likelihood[out_indxs]

;            mwrfits,out_full,out_filename,/create

        
            tarray = fltarr(n_out_indxs,2)
            tarray[*,0] = out_indxs
            tarray[*,1] = log_likelihood[out_indxs]
            fxhmake,header,tarray,/initialize
            sxaddpar,header,'EXTNAME','FULL_LHD_' + strtrim(string(out_table[cur_tot_out].number),2),'full 4D log likelihood'
            sxaddpar,header,'STARNAME',out_table[cur_tot_out].name,'name of star'
            sxaddpar,header,'STARNUM',out_table[cur_tot_out].number,'number of star'
            fits_write,ofcb_full,tarray,header
        endif

    endif

    ; make sure to close the last file if it is still open
    if (((cur_tot_out+1) EQ n_out_chunk) OR (z GT (max_i-skip_i))) then begin
        if (not keyword_set(silent)) then $
          print,'outputing chuck = ' + strtrim(string(cur_chunk+1),2)

        out_filename = outpath+outbase+'_chunk_'+strtrim(string(cur_chunk+1),2)+'.fits'

        ; trim the unused records (for the last chunk only)
        if (keyword_set(output_1d_likelihood) AND (cur_tot_out LT n_out_chunk)) then begin
            out_table = out_table[0:cur_tot_out] 
            out_1d_logt_likelihood = out_1d_logt_likelihood[*,0:cur_tot_out+1]
            out_1d_logg_likelihood = out_1d_logg_likelihood[*,0:cur_tot_out+1]
            out_1d_av_likelihood = out_1d_av_likelihood[*,0:cur_tot_out+1]
            out_1d_rv_likelihood = out_1d_rv_likelihood[*,0:cur_tot_out+1]
            out_1d_mass_likelihood = out_1d_mass_likelihood[*,0:cur_tot_out+1]
            out_1d_bmass_likelihood = out_1d_bmass_likelihood[*,0:cur_tot_out+1]
            out_1d_age_likelihood = out_1d_age_likelihood[*,0:cur_tot_out+1]
        endif

        ; write the main output file
        mwrfits,out_table,out_filename,/create

        ; write the 1d likelihood file
        if (keyword_set(output_1d_likelihood)) then begin
            fits_open,repstr(out_filename,'.fits','_1d_likelihood.fits'),ofcb_1d,/write
            
            fxhmake,header,out_1d_logt_likelihood,/initialize
            sxaddpar,header,'EXTNAME','LOGT_LHD','Log(Teff) Likelihood'
            sxaddpar,header,'COMMENT','1st row gives the x values'
            fits_write,ofcb_1d,out_1d_logt_likelihood,header

            fxhmake,header,out_1d_logg_likelihood,/initialize
            sxaddpar,header,'EXTNAME','LOGG_LHD','Log(g) Likelihood'
            sxaddpar,header,'COMMENT','1st row gives the x values'
            fits_write,ofcb_1d,out_1d_logg_likelihood,header

            fxhmake,header,out_1d_av_likelihood,/initialize
            sxaddpar,header,'EXTNAME','AV_LHD','A(V) Likelihood'
            sxaddpar,header,'COMMENT','1st row gives the x values'
            fits_write,ofcb_1d,out_1d_av_likelihood,header

            fxhmake,header,out_1d_rv_likelihood,/initialize
            sxaddpar,header,'EXTNAME','RV_LHD','R(V) Likelihood'
            sxaddpar,header,'COMMENT','1st row gives the x values'
            fits_write,ofcb_1d,out_1d_rv_likelihood,header

            fxhmake,header,out_1d_bmass_likelihood,/initialize
            sxaddpar,header,'EXTNAME','BMAS_LHD','Birth Mass Likelihood'
            sxaddpar,header,'COMMENT','1st row gives the x values'
            fits_write,ofcb_1d,out_1d_bmass_likelihood,header

            fxhmake,header,out_1d_mass_likelihood,/initialize
            sxaddpar,header,'EXTNAME','MASS_LHD','Mass Likelihood'
            sxaddpar,header,'COMMENT','1st row gives the x values'
            fits_write,ofcb_1d,out_1d_mass_likelihood,header

            fxhmake,header,out_1d_age_likelihood,/initialize
            sxaddpar,header,'EXTNAME','AGE_LHD','log(Age) Likelihood'
            sxaddpar,header,'COMMENT','1st row gives the x values'
            fits_write,ofcb_1d,out_1d_age_likelihood,header

            fits_close,ofcb_1d
        endif

        ; close the full likelihood file
        if (keyword_set(output_full_prob)) then fits_close,ofcb_full
    endif

endfor

end
