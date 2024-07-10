
program define pmcstat, rclass

/* Syntax
	VARLIST = A list of two variables, the linear predictor for the model,
			and the event indicator (observed outcome)
	NOPRINT = suppress the onscreen output of performance stats
	MATRIX = specify the name of a matrix storing the performance stats 
*/

syntax varlist(min=1 max=2 numeric) [if] [in], [noPRINT  ///
				MATrix(name local) HANley FASTER]

*********************************************** SETUP/CHECKS
*SET UP TEMPs
tempvar p rank_disc rank2_disc diff_disc inv_outcome rank_cord rank2_cord diff_cord

// check on the if/in statement 
marksample touse
qui count if `touse'
local samp=r(N)
if `r(N)'==0 { 
	di as err "if statement identifies subgroup with no data?"
	error 2000
	}
	
// parse varlist
tokenize `varlist' , parse(" ", ",")
local lp = `"`1'"'
local outcome = `"`2'"'

// generate probabilities
qui gen `p' = exp(`lp')/(1+exp(`lp'))

// run checks on user input variables in varlist
// check if user has input both LP and obs (for binary outcome)
local varcountcheck: word count `varlist'

if `varcountcheck'!=2 {
	di as err "Varlist must contain two variables. Linear predictor values, followed by observed outcomes (binary variable) are required"
	error 102
	}

// check outcome is binary
cap assert `outcome'==0 | `outcome'==1 if `touse'
        if _rc~=0 {
                noi di as err "Event indicator `outcome' must be coded 0 or 1"
                error 450
        }

// preserve data & keep only the touse sample
preserve
qui keep if `touse'

*********************************************** C-STAT

if "`faster'"=="faster" {
	// check for packages 
	local packs gtools
	foreach pkg of local packs {
		capture which `pkg'
		if _rc==111 {	
			ssc install `pkg'
			}
		}
		
	// discordant pairs
	hashsort `p' `outcome' 		
	qui gen `rank_disc' = _n if `touse'

	hashsort `outcome' `p' `rank_disc'
	qui gen `rank2_disc' = _n if `touse'

	qui gen `diff_disc' = (`rank_disc' - `rank2_disc') if (`outcome'==0) & (`touse')

	// concordant pairs
	qui gen `inv_outcome' = (`outcome'==0)
	hashsort `p' `inv_outcome'
	qui gen `rank_cord' = _n if `touse'

	hashsort `inv_outcome' `p' `rank_cord'
	qui gen `rank2_cord' = _n if `touse'

	qui gen `diff_cord' = (`rank_cord' - `rank2_cord') if `inv_outcome'==0

	// total possible pairs
	qui gstats sum `outcome' if (`outcome'!=.) & (`touse'), meanonly
	local obs = r(N)
	local prev = r(mean)
	local events = r(sum)
	local nonevents = r(N) - r(sum)
	local pairs = `events'*`nonevents'  

	// compute c-stat (allowing for ties)
	qui gstats sum `diff_disc' if `touse'
	local disc = r(sum)
	qui gstats sum `diff_cord' if `touse'
	local cord = r(sum)

	local ties = `pairs'-`disc'-`cord'

	local cstat = (`cord'+(0.5*`ties'))/(`pairs')
	}
	else {
		// discordant pairs
		sort `p' `outcome' 		
		qui gen `rank_disc' = _n if `touse'

		sort `outcome' `p' `rank_disc' 
		qui gen `rank2_disc' = _n if `touse'

		qui gen `diff_disc' = (`rank_disc' - `rank2_disc') if (`outcome'==0) & (`touse')

		// concordant pairs
		qui gen `inv_outcome' = (`outcome'==0) if `touse'
		sort `p' `inv_outcome' 
		qui gen `rank_cord' = _n if `touse'

		sort `inv_outcome' `p' `rank_cord' 
		qui gen `rank2_cord' = _n if `touse'

		qui gen `diff_cord' = (`rank_cord' - `rank2_cord') if (`inv_outcome'==0) & (`touse')

		// total possible pairs
		qui su `outcome' if (`outcome'!=.) & (`touse'), meanonly
		local obs = r(N)
		local prev = r(mean)
		local events = r(sum)
		local nonevents = r(N) - r(sum)
		local pairs = `events'*`nonevents'  

		// compute c-stat (allowing for ties)
		qui su `diff_disc' if `touse'
		local disc = r(sum)
		qui su `diff_cord' if `touse'
		local cord = r(sum)

		local ties = `pairs'-`disc'-`cord'

		local cstat = (`cord'+(0.5*`ties'))/(`pairs')
	}
	
***************************************** CI


if "`hanley'"=="" {
	// default use newcombe SE formula
	local newcombe_c = `cstat'
	local cstat_se = ((`cstat'*(1-`cstat'))*(1+(((`obs'/2)-1)*((1-`cstat')/(2-`cstat'))) ///
	+((((`obs'/2)-1)*`cstat')/(1+`cstat')))/((`obs'^2)*`prev'*(1-`prev')))^.5
	local cstat_lb = `cstat' - (1.96*`cstat_se')
	local cstat_ub = `cstat' + (1.96*`cstat_se')
}
else {
	// if hanley option set then use hanley SE formula 
	local Q1 = `cstat' / (2 - `cstat')
	local Q2 = 2 * `cstat'^2 / (1 + `cstat')
	local hanley_c = `cstat'
	local cstat_se = sqrt((`cstat' * (1 - `cstat') + (`nonevents' - 1) * (`Q1' - `cstat'^2) + (`events' - 1) * (`Q2' - `cstat'^2)) / (`nonevents' * `events'))
	local cstat_lb = `cstat' - (1.96*`cstat_se')
	local cstat_ub = `cstat' + (1.96*`cstat_se')
}


***************************************** OUTPUT

// Creating matrix of results
local res cstat 

	tempname rmat
	matrix `rmat' = J(1,5,.)
	local i=0
	foreach r of local res {
		local ++i
		matrix `rmat'[`i',1] = `obs'
		matrix `rmat'[`i',2] = ``r''
		matrix `rmat'[`i',3] = ``r'_se'
		matrix `rmat'[`i',4] = ``r'_lb'
		matrix `rmat'[`i',5] = ``r'_ub'

		
		}
		mat colnames `rmat' = Obs Estimate SE Lower_CI Upper_CI
		mat rownames `rmat' = "C-Statistic" //`rown'
		
// print matrix 
if "`matrix'"!="" {
			matrix `matrix' = `rmat'
			
			if "`print'"!="noprint" {
				
				matlist `matrix', border(all) 
							
				}
				
			}
			else { 
				if "`print'"!="noprint" {
					
					matlist `rmat', border(all) 
							
					}
				
				}
				
// Return scalars
local res cstat cstat_se cstat_lb cstat_ub cord disc ties  obs
 
		foreach r of local res {
			return scalar `r' = ``r''
			}
			
		if "`matrix'"!="" {
		    matrix `matrix' = `rmat'
			return matrix `matrix' = `rmat'
		}
		else {
		    return matrix rmat = `rmat'
		}
		
restore

end
  