
**************************************************
*												 *
*  PROGRAM TO CALCULATE C-STAT					 *
*  16/06/21 									 *
*  			 									 *
*	Updated: 17/01/2025							 *
*	- moved to using frames 					 *
*	- substantial speed gains					 *
*												 *
*  1.0.1 J. Ensor								 *
**************************************************

*! 1.0.1 J.Ensor 17Jan2025


program define pmcstat, rclass

version 16

/* Syntax
	VARLIST = A list of two variables, the linear predictor for the model,
			and the event indicator (observed outcome)
	NOPRINT = suppress the onscreen output of performance stats
	MATRIX = specify the name of a matrix storing the performance stats 
*/

syntax varlist(min=1 max=2 numeric) [if] [in], [noPRINT  ///
				MATrix(name local) HANley]

*********************************************** SETUP/CHECKS
*SET UP TEMPs
tempvar p rank_disc rank2_disc diff_disc inv_outcome rank_cord rank2_cord diff_cord

// store current frame
local curframe = c(frame)

frame `curframe' {
	
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

// move to new frame
tempname pmcstat_frame_105
frame put `lp' `outcome' if `touse', into(`pmcstat_frame_105')
frame change `pmcstat_frame_105'

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
cap assert `outcome'==0 | `outcome'==1 
        if _rc~=0 {
                noi di as err "Event indicator `outcome' must be coded 0 or 1"
                error 450
        }


*********************************************** C-STAT

// discordant pairs
sort `p' `outcome' 		
qui gen `rank_disc' = _n 

sort `outcome' `p' `rank_disc' 
qui gen `rank2_disc' = _n 

qui gen `diff_disc' = (`rank_disc' - `rank2_disc') if (`outcome'==0) 

// concordant pairs
qui gen `inv_outcome' = (`outcome'==0) 
sort `p' `inv_outcome' 
qui gen `rank_cord' = _n 

sort `inv_outcome' `p' `rank_cord' 
qui gen `rank2_cord' = _n 

qui gen `diff_cord' = (`rank_cord' - `rank2_cord') if (`inv_outcome'==0) 

// total possible pairs
qui su `outcome' if (`outcome'!=.), meanonly 
local obs = r(N)
local prev = r(mean)
local events = r(sum)
local nonevents = r(N) - r(sum)
local pairs = `events'*`nonevents'  

// compute c-stat (allowing for ties)
qui su `diff_disc' 
local disc = r(sum)
qui su `diff_cord' 
local cord = r(sum)

local ties = `pairs'-`disc'-`cord'

local cstat = (`cord'+(0.5*`ties'))/(`pairs')
	
	
***************************************** CI
/*
local logit_c = logit(`cstat')

local var_logit_c = (1+(`obs'/2-1)*(1-`cstat')/(2-`cstat')+(`obs'/2-1)*`cstat'/(1+`cstat'))/(`cstat'*(1-`cstat')*`events'*(`obs'-`events'))

local logit_c_se = `var_logit_c'^.5
local logit_c_lb = `logit_c' - (1.96*`logit_c_se')
local logit_c_ub = `logit_c' + (1.96*`logit_c_se')

local cstat_se = (`var_logit_c'*(`cstat'*(1-`cstat'))^2)^.5 // incorrect - formula should use var(c) but we do not have this - see debray appendix eq.55
local cstat_lb = invlogit(`logit_c_lb')
local cstat_ub = invlogit(`logit_c_ub')

local norm_c = `cstat'
local norm_c_se = ((`cstat'*(1-`cstat'))/`obs')^.5
local norm_c_lb = `cstat' - (1.96*`norm_c_se')
local norm_c_ub = `cstat' + (1.96*`norm_c_se')

local newcombe_c = `cstat'
local newcombe_c_se = ((`cstat'*(1-`cstat'))*(1+(((`obs'/2)-1)*((1-`cstat')/(2-`cstat'))) ///
+((((`obs'/2)-1)*`cstat')/(1+`cstat')))/((`obs'^2)*`prev'*(1-`prev')))^.5
local newcombe_c_lb = `cstat' - (1.96*`newcombe_c_se')
local newcombe_c_ub = `cstat' + (1.96*`newcombe_c_se')

local Q1 = `cstat' / (2 - `cstat')
local Q2 = 2 * `cstat'^2 / (1 + `cstat')
local hanley_c = `cstat'
local hanley_c_se = sqrt((`cstat' * (1 - `cstat') + (`nonevents' - 1) * (`Q1' - `cstat'^2) + (`events' - 1) * (`Q2' - `cstat'^2)) / (`nonevents' * `events'))
local hanley_c_lb = `cstat' - (1.96*`hanley_c_se')
local hanley_c_ub = `cstat' + (1.96*`hanley_c_se')
*/

if "`hanley'"=="" {
	// default use necombe SE formula
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
		mat rownames `rmat' = "C-Statistic" 

		
// print matrix 
if "`matrix'"!="" {
			matrix `matrix' = `rmat'
			
			//return matrix `matrix' = `rmat' 
			if "`print'"!="noprint" {
				//di as res _n "Discrimination statistics ..."
				matlist `matrix', border(all) 
							
				}
				
			}
			else { 
				if "`print'"!="noprint" {
					//di as res _n "Discrimination statistics ..."
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


}

end
