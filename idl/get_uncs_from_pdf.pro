
function get_uncs_from_pdf, vals, prob, best_val

indxs = where(prob GT 0.0,n_indxs)

if (n_indxs GT 0) then begin
    cum_prob = total(prob[indxs],/CUMULATIVE)
    cum_prob /= max(cum_prob)
    best_val_cum_prob = interpol(cum_prob,vals[indxs],best_val)
    unc_probs = best_val_cum_prob + [-1.0,1.0]*0.67/2.
    unc_probs[0] = max([0.001,unc_probs[0]])
    unc_probs[1] = min([0.999,unc_probs[1]])
    mm_val = interpol(vals[indxs],cum_prob,unc_probs)
    
    uncs = [best_val - mm_val[0],mm_val[1] - best_val]
    uncs[0] = max([uncs[0],0.0])
    uncs[1] = max([uncs[1],0.0])
    
    uncs = [total(uncs)/2.0,uncs]
    
    bindxs = where(uncs LT 0.0,n_bindxs)
    if (n_bindxs GT 0) then stop

endif else begin

    uncs = [0.0,0.0,0.0]

endelse

;print,best_val,uncs,mm_val

return, uncs

end
