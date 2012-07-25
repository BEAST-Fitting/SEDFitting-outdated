pro read_phat_cat,filename,fluxes,fluxes_unc,coords,all_mags,all_mags_unc,ignoreir=ignoreir,noir=noir, $
                  ignoreuv=ignoreuv,nouv=nouv,use_all=use_all

; old filename = dustin_corr3_all.fits

;table = mrdfits('Catalogs/'+filename,1,range=[1,20])
table = mrdfits(filename,1)

if (keyword_set(use_all)) then begin
    n_indxs = n_elements(table[*].acs_mag1)
    indxs = indgen(n_indxs)

;    if (n_indxs GT 300) then begin
;        nindxs = fix(randomu(seed,300)*n_indxs)
;        indxs = indxs[nindxs]
;        print,nindxs
;        n_indxs = n_elements(indxs)
;    endif
endif else begin
    indxs = where((table[*].acs_nmatch GT 0) AND (table[*].ir_nmatch GT 0) AND (table[*].uv_nmatch GT 0),n_indxs)

    indxs_noir = where((table[*].acs_nmatch GT 0) AND (table[*].uv_nmatch GT 0) AND (table[*].ir_nmatch LE 0),n_indxs_noir)

    indxs_nouv = where((table[*].acs_nmatch GT 0) AND (table[*].uv_nmatch LE 0) AND (table[*].ir_nmatch GT 0),n_indxs_nouv)
    
    indxs_noopt = where((table[*].acs_nmatch LE 0) AND (table[*].uv_nmatch GT 0) AND (table[*].ir_nmatch GT 0),n_indxs_noopt)

    indxs = [indxs,indxs_noir,indxs_nouv]
    n_indxs = n_indxs + n_indxs_noir + n_indxs_nouv
endelse

;indxs = [indxs,indxs_nouv]
;n_indxs = n_indxs + n_indxs_nouv

;print,n_indxs,n_indxs_noir,n_indxs_noopt,n_indxs_nouv
;stop

;indxs41 = where((table[*].acs_index NE -1) AND (table[*].uv_index NE -1),n_indxs41)

;indxs42 = where((table[*].ir_index NE -1) AND (table[*].uv_index NE -1),n_indxs42)

;indxs43 = where((table[*].acs_index NE -1) AND (table[*].ir_index NE -1),n_indxs43)

if (keyword_set(noir)) then begin
    indxs = indxs41
    n_indxs = n_indxs41

;    indxs2 = where(table[indxs].acs_mag2 LT 22.0,n_indxs2)
;    indxs = indxs[indxs2]
;    n_indxs = n_indxs2
endif else if (keyword_set(nouv)) then begin
    indxs = indxs43
    n_indxs = n_indxs43

    indxs2 = where(table[indxs].acs_mag2 LT 20.0,n_indxs2)
    indxs = indxs[indxs2]
    n_indxs = n_indxs2
endif

;print,n_indxs,n_indxs41,n_indxs42,n_indxs43

if (keyword_set(ignoreir) OR keyword_set(noir)) then begin
    es_points = [3.2596e-18,2.5327e-19,1.820979E-19,7.033186E-20]
    vegamag = [22.65,25.79,26.16,25.52]
    waves = [275.,336.,475.,814.]*1e-3
endif else if (keyword_set(ignoreuv) OR keyword_set(nouv)) then begin
    es_points = [1.820979E-19,7.033186E-20,1.5233e-20,1.9106e-20]
    vegamag = [26.16,25.52,26.07,24.70]
    waves = [475.,814.,1100.,1600.]*1e-3
endif else begin
; wrong zero points for F336W (sigh-again)
;    es_points = [3.2596e-18,2.5327e-19,1.820979E-19,7.033186E-20,1.5233e-20,1.9106e-20]
;    vegamag = [22.65,25.79,26.16,25.52,26.07,24.70]
    es_points = [3.2596e-18,1.34e-18,1.820979E-19,7.033186E-20,1.5233e-20,1.9106e-20]
    vegamag = [22.65,23.46,26.16,25.52,26.07,24.70]
    waves = [275.,336.,475.,814.,1100.,1600.]*1e-3
endelse

zero_points = es_points/10^(-0.4*vegamag)
;n_indxs = 10
for i = 0L,(n_indxs-1) do begin
    k = indxs[i]
    if (keyword_set(ignoreir) OR keyword_set(noir)) then begin
        mags = [table[k].uv_mag1,table[k].uv_mag2, $
                table[k].acs_mag1,table[k].acs_mag2]
        mag_uncs = [table[k].uv_mag1err,table[k].uv_mag2err, $
                    table[k].acs_mag1err,table[k].acs_mag2err]
    endif else if (keyword_set(ignoreuv) OR keyword_set(nouv)) then begin
        mags = [table[k].acs_mag1,table[k].acs_mag2, $
                table[k].ir_mag1,table[k].ir_mag2]
        mag_uncs = [table[k].acs_mag1err,table[k].acs_mag2err, $
                    table[k].ir_mag1err,table[k].ir_mag2err]
    endif else begin
        mags = [table[k].uv_mag1,table[k].uv_mag2, $
                table[k].acs_mag1,table[k].acs_mag2, $
                table[k].ir_mag1,table[k].ir_mag2]
        mag_uncs = [table[k].uv_mag1err,table[k].uv_mag2err, $
                    table[k].acs_mag1err,table[k].acs_mag2err, $
                    table[k].ir_mag1err,table[k].ir_mag2err]
    endelse

    one_fluxes = zero_points*10^(-0.4*mags)
    one_fluxes_unc = 0.5*(abs(zero_points*10^(-0.4*(mags-mag_uncs)) - one_fluxes) + $
      abs(zero_points*10^(-0.4*(mags+mag_uncs)) - one_fluxes))

    if (i EQ 0) then begin
        n_waves = n_elements(waves)
        fluxes = fltarr(n_waves,n_indxs)
        fluxes_unc = fluxes
        all_mags = fluxes
        all_mags_unc = fluxes
        coords = dblarr(2,n_indxs)
    endif
    fluxes[*,i] = one_fluxes
    fluxes_unc[*,i] = sqrt(one_fluxes_unc^2 + (replicate(0.1,6)*one_fluxes)^2)
    all_mags[*,i] = mags
    all_mags_unc[*,i] = mag_uncs
    if (keyword_set(use_all)) then begin
        coords[*,i] = [0d0,0d0]
;        fluxes[4:5,i] = -99.
;        fluxes_unc[4:5,i] = -99.
    endif else begin
        coords[*,i] = [table[k].ra,table[k].dec]
    endelse

;    kplot,waves,mags,yerror=mag_uncs,psym=1
;    kplot,waves,fluxes,yerror=fluxes_unc,psym=1

;    ans = ''
;    read,'continue: ',ans
endfor

end
