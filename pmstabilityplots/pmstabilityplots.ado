
program define pmstabilityplots, rclass

version 16

/* Syntax
	VARLIST = A list of variables in the linear predictor for the model
	 
*/

syntax varlist , [Boot(int 200) ///
				NOSTability ///
				CMDline(string asis) ///
				COLor(string)]


*********************************************** SETUP/CHECKS
* store current frame

local curframe = c(frame)

frame `curframe' {

if "`color'"=="" {
	local color "gs10"
}

frame copy `curframe' pmstabilityplots
frame copy `curframe' boot 

qui `e(cmdline) '
frame pmstabilityplots: predict lp_og, xb
frame pmstabilityplots: predict pr_og, pr

qui forvalues i=1/`boot' {
frame change boot

preserve

* take a bootstrap sample (random sample with replacement) of same size
bsample

* fit model
`e(cmdline) '

* make predictions in OG sample
frame pmstabilityplots: predict lp`i', xb
frame pmstabilityplots: predict pr`i', pr
 
restore

local mata_var_list = "`mata_var_list' pr`i'"

* iteration indicator to help the user know the progress
nois _dots `i' 0
}

frame change pmstabilityplots

// this creates a view onto the data in mata, it will only use the vars listed above
mata: st_view(D=., ., "`mata_var_list'")



// Now we use the program above to calculate the UI based on our data view
// set the centiles we want in a vector p
mata: p = (0.025, 0.975)

// Run the program 
mata: percentiles = row_percentiles(D, p)

// Send matrix to stata
mata: st_matrix("uncertaintyIntervals", percentiles)


// gen variables for UI from the matrix created using mata program 
svmat uncertaintyIntervals 

rename uncertaintyIntervals1 lower
rename uncertaintyIntervals2 upper

* instability plots
local instability_plots 

forvalues rep = 1/`boot' {	
	local instability_plots `instability_plots'  (scatter pr`rep' pr_og, sort jitter(3) msize(vtiny) mcol("`color'")) 
	//nois _dots `rep' 0
}

local n = _N

twoway `instability_plots' (lowess upper pr_og, sort lcol(black) lpattern(dash) bwidth(0.2) plotr(lcol(black)) text(1 0 "N = `n'", size("medsmall") place(se) just(left))) ///
(lowess lower pr_og, sort lcol(black) lpattern(dash) bwidth(0.2)) ///
|| function y = x, clpat(solid) clcol(black) legend(off) ///
xlab(#5, angle(h) grid nogextend format(%3.1f) labsize("medsmall")) ylab(#5, angle(h) grid nogextend format(%3.1f) labsize("medsmall")) xtitle("Estimated risk from developed model", size("medsmall")) ytitle("Estimated risk from bootstrap models", size("medsmall")) aspect(1) graphr(col(white)) name(p_instability_plot, replace)


}

end

mata:
// Assume X is your matrix with rows as observations and columns as variables
function row_percentiles(X, p) {
    n_rows = rows(X)
    n_cols = cols(X)
    result = J(n_rows, 2, .)  // Empty matrix to store 2.5th and 97.5th percentiles
    
    for (i=1; i<=n_rows; i++) {
        x_row = X[i,.]'  // Extract row and transpose to column
        x_sorted = sort(x_row, 1)  // Sort values
        
        // Calculate percentiles
        for (j=1; j<=cols(p); j++) {
            k = n_cols * p[j]
            k1 = floor(k)
            k2 = ceil(k)
            
            if (k1==k2) {
                result[i,j] = x_sorted[k1]
            } else {
                // Linear interpolation
                result[i,j] = x_sorted[k1] + (x_sorted[k2]-x_sorted[k1])*(k-k1)
            }
        }
    }
    return(result)
}
end