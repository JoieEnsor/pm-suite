
program define pmstabilityss, rclass

version 16

/* Syntax
	VARLIST = A list of variables in the linear predictor for the model
	 
*/

syntax varlist , PREValence(real)   ///
			[ALPHA(real 0) ///
			LP(varname max=1 numeric) ///
			SUBgroup(varname max=1) ///
			NODRAW ///
			NOCOMPare ///
			COLor(string) ///
			STandardised ///
			PMSS ///
			N(numlist) ///
			BETAs(numlist) ///
			CSTATistic(real 0)  ///
			SIMOBS(int 500000)  ///
			SEED(int 1234)  ///
			PCIwidth(numlist) ///
			PCUTpoints(numlist ascending >=0 <=1) /// 	 
			THRESHold(real 0) ///
			TOLerance(real 0.005) ///
			logitpvarincrement(real 0.001) ///
			saving(string)]

*********************************************** SETUP/CHECKS
* store current frame

local curframe = c(frame)

frame `curframe' {

if "`color'"=="" {
	local color "gs10"
}


*******************************
// check number of vars match up with betas
local vars = 0

local b = 1
 foreach a in `varlist' {
 	local vars = `vars'+1
	
	// set varname locals for creating fishers I matrix
	local b = `b'+1
 	local varname`b' = "`a'"
	
 }
 
 
 if "`lp'"=="" {
 	if "`standardised'"!="standardised" {
	 local beta_count = 0
	 foreach b in `betas' {
		local beta_count = `beta_count'+1
	 }
	 
	 if `vars'!=`beta_count' {
		di as err "Number of beta effects does not match number of input variables"
		error 459
	 }
 }
 }
 
local cuts: word count `pcutpoints'
local user_widths: word count `pciwidth'

if `cuts'!=`user_widths' {
	di as err "Number of cutpoints must match number of target interval widths - note pcutpoints() numlist should end with 1" 
		error 459
}
 
// check we only have continuous or binary variables


// check subgroup var format 


*******************************
// create new frame summarising the user dataset for input into pmsim

if "`lp'"=="" { 	// start of loop for if user has not specified LP directly

tempname pmstabilityss_sumstats
cap frame drop `pmstabilityss_sumstats'
frame create `pmstabilityss_sumstats' str25(predictors) double(beta prop mean sd)

if "`standardised'"=="standardised" {


foreach name in `varlist' {
	
	qui levelsof `name'
	if r(r)==2 { 
		// binary vars
		cap assert `name' == 0 | `name' == 1
		if _rc==9 {
			di as err "Binary variables must be coded 0/1"
			error 459
		}
		qui su `name'
		frame post `pmstabilityss_sumstats' ("`name'") (1) (r(mean)) (.) (.)
		
		
// 			qui replace `name' = (`name' - r(mean))/r(sd)
		
	}
	else { 
		// continuous vars
		qui su `name'
		
			qui replace `name' = (`name' - r(mean))/r(sd)
	
		
		qui su `name'
		frame post `pmstabilityss_sumstats' ("`name'") (1) (.) (r(mean)) (r(sd))
		
	}
	
}
}
else {
	tokenize `betas'

foreach name in `varlist' {
	local token = `token' +1
	qui levelsof `name'
	if r(r)==2 { 
		// binary vars
		cap assert `name' == 0 | `name' == 1
		if _rc==9 {
			di as err "Binary variables must be coded 0/1"
			error 459
		}
		qui su `name'
		frame post `pmstabilityss_sumstats' ("`name'") (``token'') (r(mean)) (.) (.)
		
		
	}
	else { 
		// continuous vars
		
		
		qui su `name'
		frame post `pmstabilityss_sumstats' ("`name'") (``token'') (.) (r(mean)) (r(sd))
		
	}
	
}
}


frame change `pmstabilityss_sumstats'
// frame rename `pmstabilityss_sumstats' test
mkmat beta , mat(input_betas) rown(predictors)

*******************************
// run iterative proc to get alpha/scalar required using pmsim/prev

// parse whether user input prev &/or c-stat

if "`standardised'"=="standardised" {
	if `cstatistic'!=0 {
		di "Running pmsim"
		
		pmsim , prev(`prevalence') alpha(`alpha') cstat(`cstatistic')  simobs(`simobs') pred(predictors) beta(beta) prop(prop) mean(mean) sd(sd)  tol(`tolerance') seed(`seed') standardised
		
		local intercept = r(intercept)
		local scalar = r(scalar)
	}
	else {
		di "Running pmsimprev"
		
		pmsimprev , prev(`prevalence') alpha(`alpha') simobs(`simobs') pred(predictors) beta(beta) prop(prop) mean(mean) sd(sd)  tol(`tolerance') seed(`seed')  standardised
		
		local intercept = r(intercept)
		local scalar = 1
	}
}
else {	
	if `cstatistic'!=0 {
		di "Running pmsim"
		
		pmsim , prev(`prevalence') alpha(`alpha') cstat(`cstatistic')  simobs(`simobs') pred(predictors) beta(beta) prop(prop) mean(mean) sd(sd)  tol(`tolerance') seed(`seed')
		
		local intercept = r(intercept)
		local scalar = r(scalar)
	}
	else {
		di "Running pmsimprev"
		
		pmsimprev , prev(`prevalence') alpha(`alpha') simobs(`simobs') pred(predictors) beta(beta) prop(prop) mean(mean) sd(sd)  tol(`tolerance') seed(`seed') 
		
		local intercept = r(intercept)
		local scalar = 1
	}
}

*************************************
// back to default frame to generate LP variable

frame change `curframe'

// create PI 

tempvar pi

 forvalues a=1/`vars' {
 	local b = `a'+1
	
	if `a'==1 {
		qui gen pi = input_betas[`a',1]*`varname`b''
		
	}
	else {
		qui replace pi = pi + input_betas[`a',1]*`varname`b''
		
	}
	
 	}
	
*************************************
// create LP using bisection alpha/scalar

tempvar lp

qui gen lp = `intercept' + `scalar'*(pi)

} // end of loop for if user has not specified LP directly
else {
	qui gen lp = `lp'
}

********************************
// Calculate Information matrix elements
tempvar id
qui gen `id' = _n

local elements = `vars'+1

// create proxy variable for intercept 
local varname1 = "intercept"

tempvar intercept
qui gen intercept = 1


// set up empty vector to fill
local stack_rows = (`elements' + `elements'^2)/2
matrix stack = J(`stack_rows',1,.)

// compute off-diagonals
local row_start = 1
local row_iter = 0
forvalues col=1/`elements' {
	forvalues row=`row_start'/`elements' {
		local row_iter = `row_iter'+1
		qui gen M_`row'`col' = (`varname`row'')*(`varname`col'')*(exp(lp)/((1+exp(lp))^2))
		qui su M_`row'`col'
		matrix stack[`row_iter',1] = r(mean)
		qui drop M_`row'`col'
	}
	local row_start = `row_start'+1
}

mata: st_matrix("I", invvech(st_matrix("stack")))

// inverse
mat invI = inv(I)

// Calculate individuals unit variance

local mata_var_list "intercept `varlist'"
mata: D = st_data(.,("`mata_var_list'"))

mata: invIm = st_matrix("invI")

mata: unit_v = J(rows(D),1,.)

qui su intercept
forvalues i = 1/`r(N)' {
	mata: unit_v[`i',1] = D[`i',.]*invIm*D[`i',.]'
}

mata: st_matrix("unit_v", unit_v)

svmat unit_v


********************************


*********************************
// gen predicted probabilities
qui gen p_true = invlogit(lp)

local cur_data_n = _N
di as txt _n "Fixed SS of input dataset = `cur_data_n'" _n
********************************
// parse user cutpoints for p and desired widths 

local cuts: word count `pcutpoints'

if `cuts'!=0 {
	
mat widths = J(`cuts',4,.)

local q = 1
foreach i in `pcutpoints' {
	matrix widths[`q',1] = `i'
	local ++q
	}
	
	local q = 1
foreach i in `pciwidth' {
	matrix widths[`q',2] = `i'
	matrix widths[`q',3] = `q'
	local ++q
	}
	
mat colnames widths = "categories" "width" "cat_code" "avg_p"

qui gen width = .
qui gen cat_code = .
forvalues i = 1/`cuts' {
	qui replace width = widths[`i',2] if p_true<=widths[`i',1] & width==.
	qui replace cat_code = widths[`i',3] if p_true<=widths[`i',1] & cat_code==.
}

qui gen avg_p = .
forvalues i = 1/`cuts' {
	qui su p_true if cat_code==`i' , det
	qui replace avg_p = `r(p50)' if cat_code==`i'
	mat widths[`i',4] = `r(p50)'
}


***************************************	
// matrix of values for p groupings
mat pvars = J(100,4,.)
mat colnames pvars = "pcode" "pcat" "width" "logit_p_var"

mata: st_matrix("widths_desc",sort(st_matrix("widths"),-1))

qui gen double p_round = round(p_true, 0.01)
qui replace p_round = .01 if p_round==0
qui replace p_round = .99 if p_round==1
qui gen var_logit_p = .

local j = 0
forvalues i = 0.01 (0.01) 1.01 {
	local j = `j'+1
	mat pvars[`j',1] = `j'
	mat pvars[`j',2] = `i'
		
	
	forvalues k = 1/`cuts' {
		if round(pvars[`j',2],0.001) <= round(widths_desc[`k',1],0.001) {
			mat pvars[`j',3] = widths_desc[`k',2] 
			}
	}

	// identify appropriate target var(logit(p)) given pciwidth
	 local width_p = 0 
	 local var_logit_p = 0
	 
	 
	 	while `width_p' < pvars[`j',3] {
			local var_logit_p = `var_logit_p'+`logitpvarincrement'
			local ub_p = invlogit(((logit(pvars[`j',2]) + (1.96*sqrt(`var_logit_p')))))
			local lb_p = invlogit(((logit(pvars[`j',2]) - (1.96*sqrt(`var_logit_p')))))
			local width_p = `ub_p'-`lb_p'
		
	 }
	 mat pvars[`j',4] = `var_logit_p' 
	 qui replace var_logit_p = `var_logit_p' if p_round==round(`i',0.01)
}


********************************
// SS required to meet desired precision
qui gen ss_target_var = ((1/var_logit_p)*(unit_v))
qui su ss_target_var
local target_var_min_N = ceil(r(max))

di as txt "Minimum SS required to meet target UI widths = `target_var_min_N'" _n

} // end of if clause for when user inputs target widths

if "`pmss'"!="" {
	// check for packages 
	local package pmsampsize
	foreach pack of local package {
		capture which `pack'
		if _rc==111 {	
			di as error "Package 'pmsampsize' is required for option pmss" _n "Installation will begin now ..."
			ssc install `pack', replace
			}
		}
	
// work out minimum SS required from original criteria as basis for option 1 calculations
qui pmsampsize, type(b) cstatistic(`cstatistic') parameters(`vars') prevalence(`prevalence')

local pmss_min_N = r(sample_size)
di as txt "Minimum SS required by pmsampsize = `pmss_min_N'" _n
}

**********************************
// build up the list of N's to investigate
local cuts: word count `pcutpoints'
local user_widths: word count `pciwidth'
local user_ns: word count `n'
local cur_data_n = _N

local ss_tests = "`cur_data_n'"
if `cuts'!=0 {
	local ss_tests = "`ss_tests' `target_var_min_N'"
}

if "`pmss'"!="" {
	local ss_tests = "`ss_tests' `pmss_min_N'"
}

if `user_ns'!=0 {
	local ss_tests = "`ss_tests' `n'"
}


**********************************
// summary & plots for chosen N's

local total_ns: word count `ss_tests'
if `total_ns'<=2 {
	local text_size = "medsmall"
}
else {
	local text_size = "small"
}

tempname pmstabilityss_overall_stats
	cap frame drop pmstabilityss_overall_stats
	frame create pmstabilityss_overall_stats double(N Mean Min Median Max)

if `threshold'!=0 {	
	tempname pmstabilityss_threshold_stats
	cap frame drop pmstabilityss_threshold_stats
	frame create pmstabilityss_threshold_stats double(N Threshold Mean Min Median Max)
}

if `cuts'!=0 {	
	tempname pmstabilityss_cuts_stats
	cap frame drop pmstabilityss_cuts_stats
	frame create pmstabilityss_cuts_stats double(N P_category Target_width Mean Min Median Max Prop_target_width_met)

}

local cuts_rspec = "--"
local i = 0
foreach num in `ss_tests' {
		local i = `i'+1
		
qui gen var_ind_`num' = (1/`num')*unit_v

* convert to Ci on logit then p scales
qui gen lower_logitp_`num' = lp - (1.96 * sqrt(var_ind_`num'))
qui gen upper_logitp_`num' = lp + (1.96 * sqrt(var_ind_`num'))

qui gen lower_p_`num' = invlogit(lower_logitp_`num')
qui gen upper_p_`num' = invlogit(upper_logitp_`num')

qui gen width_`num' = upper_p_`num' - lower_p_`num'
* summarise widths overall 
qui summ width_`num', det

frame post pmstabilityss_overall_stats (`num') (r(mean)) (r(min)) (r(p50)) (r(max)) 

if `cuts'!=0 {
	
qui gen width_good_`num' = 1
qui replace width_good_`num' = 0 if width_`num' > width

* summarise widths overall
qui summ width_good_`num'
local met_target_width_`num' : di %9.2gc `r(mean)'

forvalues c = 1/`cuts' { 
		qui summ width_`num' if cat_code==`c', det
		
		frame post pmstabilityss_cuts_stats (`num') (widths[`c',1]) (widths[`c',2]) (r(mean)) (r(min)) (r(p50)) (r(max)) (`met_target_width_`num'')
		
		local cuts_rspec = "- `cuts_rspec' "
	}
	
}

* create z values based on a threshold value being high 
if `threshold'!=0 {
qui gen z_`num' = (logit(`threshold') - lp)/sqrt(var_ind_`num')
qui gen prob_above_`num' = 1- normal(z_`num')
qui gen prob_below_`num' = normal(z_`num')

qui gen prob_different_`num' = .
qui replace prob_different_`num' = prob_above_`num' if p_true < `threshold'
qui replace prob_different_`num' = prob_below_`num' if p_true >= `threshold'
qui su prob_different_`num', det
local textpos = r(max)

frame post pmstabilityss_threshold_stats (`num') (`threshold') (r(mean)) (r(min)) (r(p50)) (r(max)) 

/*
* summarise probs overall 
qui summ prob_above_`num' if p_true < `threshold'
qui summ prob_below_`num' if p_true >= `threshold'


* calculate the proportion of people who have a CI that crosses a threshold 
qui gen swap_`num' = 0
qui replace swap_`num' = 1 if p_true >= `threshold' & lower_p_`num' < `threshold'
qui replace swap_`num' = 1 if p_true < `threshold' & upper_p_`num' >= `threshold'
qui su swap_`num', det
*/

// classification instability plot
twoway (scatter prob_different_`num' p_true, sort `nodraw' jitter(0) msym(Oh) msize(tiny) plotr(lcol(black)) mcol(`color') legend(off) xtitle(True risk, size("`text_size'")) ytitle("Proportion of intervals" "with misclassification", size("`text_size'")) xlab(#5, angle(h) grid nogextend format(%3.1f) labsize("`text_size'")) ylab(#5, angle(h) grid nogextend format(%3.1f) labsize("`text_size'")) graphregion(col(white)) name(class_instability_`num', replace) text(`textpos' 1 "N = `num'", size("`text_size'") place(w) just(right))) //title("Mislassification instability plot for N = `num'", size(medsmall))) 

local comp_class_plot_list = "`comp_class_plot_list' class_instability_`num'"
}

* plot CI for each individual versus their true risk (prediction instability plot)
twoway (rspike lower_p_`num' upper_p_`num' p_true, `nodraw' sort jitter(0) lcol(`color') plotr(lcol(black)) text(1 0 "N = `num'", size("`text_size'") place(se) just(left))) (lowess lower_p_`num' p_true, sort lcol(black) lpattern(dash) bwidth(0.2)) (lowess upper_p_`num' p_true, sort lcol(black) lpattern(dash) bwidth(0.2)) || function y = x, clpat(solid) clcol(black) legend(off) xlab(#5, angle(h) grid nogextend format(%3.1f) labsize("`text_size'")) ylab(#5, angle(h) grid nogextend format(%3.1f) labsize("`text_size'")) xtitle(True risk, size("`text_size'")) ytitle("95% Uncertainty interval" "for true risk", size("`text_size'")) aspect(1) graphr(col(white)) name(instability_plot_`num', replace) //title("Prediction instability plot" "for N = `num'", size(medsmall)) 


local comp_insta_plot_list = "`comp_insta_plot_list' instability_plot_`num'"

}

**********************************
// create matrices of results for displyed output
frame change pmstabilityss_overall_stats 
mkmat Mean Min Median Max, mat(pmstss_overall) rown(N)

		 
matlist pmstss_overall, aligncolnames(r) lines(row) border(row) format(%9.2gc) underscore rowtitle("N") title("Overall summary UI widths")

if `threshold'!=0 {	
	frame change pmstabilityss_threshold_stats 
	mkmat Threshold Mean Min Median Max, mat(pmstss_threshold) rown(N)

	
	matlist pmstss_threshold, aligncolnames(r) lines(row) border(row) format(%9.2gc) underscore rowtitle("N") title("Summary of proportion of misclassified")
}

if `cuts'!=0 {
	frame change pmstabilityss_cuts_stats 
	mkmat P_category Target_width Mean Min Median Max Prop_target_width_met, mat(pmstss_target) rown(N)
	
			 
	matlist pmstss_target, aligncolnames(r) cspec(o2& w12 | w18 | w18 | %9.2gc | %9.2gc | %9.2gc | %9.2gc | w22 o2&) rspec("`cuts_rspec'") underscore rowtitle("N") title("Summary UI widths by probability categories")
}

**********************************
if "`nocompare'"!="nocompare" {
graph combine `comp_insta_plot_list', xcom ycom name(instability_plots_pmstss, replace)

if `threshold'!=0 {
	graph combine `comp_class_plot_list', xcom ycom name(classification_plots_pmstss, replace)
}

}

// return to current frame
frame change `curframe'

****************************
// subgroups 
if "`subgroup'"!="" {

capture confirm numeric variable `subgroup'
if !_rc {
	tempname pmstabilityss_sub_stats
	cap frame drop pmstabilityss_sub_stats
	frame create pmstabilityss_sub_stats double(Subgroup N Mean Min Median Max)
	
	foreach num in `ss_tests' {
		
	qui levelsof `subgroup' , local(sub)    
    foreach s of local sub {
		qui summ width_`num' if `subgroup'==`s', det

		frame post pmstabilityss_sub_stats (`s') (`num') (r(mean)) (r(min)) (r(p50)) (r(max))
		
		* plot CI for each individual versus their true risk (prediction instability plot)
		twoway (rspike lower_p_`num' upper_p_`num' p_true if `subgroup'==`s', `nodraw' sort jitter(0) lcol(`color') plotr(lcol(black)) text(1 0 "N = `num'" "Subgroup = `s'", size("`text_size'") place(se) just(left))) (lowess lower_p_`num' p_true if `subgroup'==`s', sort lcol(black) lpattern(dash) bwidth(0.2)) (lowess upper_p_`num' p_true if `subgroup'==`s', sort lcol(black) lpattern(dash) bwidth(0.2)) || function y = x, clpat(solid) clcol(black) legend(off) xlab(#5, angle(h) grid nogextend format(%3.1f) labsize("`text_size'")) ylab(#5, angle(h) grid nogextend format(%3.1f) labsize("`text_size'")) xtitle(True risk, size("`text_size'")) ytitle("95% Uncertainty interval" "for true risk", size("`text_size'")) aspect(1) graphr(col(white)) name(inst_`num'_`s', replace) //title("Prediction instability plot" "for N = `num'", size(medsmall)) 


	local comp_sub_`s'_list = "`comp_sub_`s'_list' inst_`num'_`s'"
	}
}
}
else {
	tempname pmstabilityss_sub_stats
	cap frame drop pmstabilityss_sub_stats
	frame create pmstabilityss_sub_stats str32(Subgroup) double(N Mean Min Median Max)

foreach num in `ss_tests' {
		
	qui levelsof `subgroup' , local(sub)    
    foreach s of local sub {
		qui summ width_`num' if `subgroup'=="`s'", det

		frame post pmstabilityss_sub_stats ("`s'") (`num') (r(mean)) (r(min)) (r(p50)) (r(max))
		
		* plot CI for each individual versus their true risk (prediction instability plot)
		twoway (rspike lower_p_`num' upper_p_`num' p_true if `subgroup'=="`s'", `nodraw' sort jitter(0) lcol(`color') plotr(lcol(black)) text(1 0 "Subgroup = `s'" "N = `num'", size("`text_size'") place(se) just(left))) (lowess lower_p_`num' p_true if `subgroup'=="`s'", sort lcol(black) lpattern(dash) bwidth(0.2)) (lowess upper_p_`num' p_true if `subgroup'=="`s'", sort lcol(black) lpattern(dash) bwidth(0.2)) || function y = x, clpat(solid) clcol(black) legend(off) xlab(#5, angle(h) grid nogextend format(%3.1f) labsize("`text_size'")) ylab(#5, angle(h) grid nogextend format(%3.1f) labsize("`text_size'")) xtitle(True risk, size("`text_size'")) ytitle("95% Uncertainty interval" "for true risk", size("`text_size'")) aspect(1) graphr(col(white)) name(inst_`num'_`s', replace) //title("Prediction instability plot" "for N = `num'", size(medsmall)) 


	local comp_sub_`s'_list = "`comp_sub_`s'_list' inst_`num'_`s'"
	}
}
}

frame change pmstabilityss_sub_stats 
	mkmat N Mean Min Median Max, mat(pmstss_subgroup) rown(Subgroup)
	
			 
	matlist pmstss_subgroup, tw(32) aligncolnames(r) lines(row) border(row) format(%9.2gc) underscore rowtitle("`subgroup'") title("Summary UI widths by subgroup - `subgroup'")
	
	
if "`nocompare'"!="nocompare" {
	foreach s of local sub {
		graph combine `comp_sub_`s'_list', xcom ycom name(instability_sub_`s', replace)
	}
}

}
***************

// saving frames
if "`saving'"!="" {
	foreach r in overall threshold cuts sub {
		qui frames save pmstss_`r'_`saving', frames(pmstabilityss_`r'_stats) replace
	}
	
}

***************
}

end



program define pmsimprev, rclass

/* Syntax
	PREVALENCE = Target prevalence required 
	ALPHA = Intercept of published/existing model. Or start point for intercept.
	SIMOBS = Number of observations to simulate super population. 
	SEED = Random number seed.
	PREDICTORS = Varname for variable specifying names of predictors in existing model.
	PROPORTION = Varname for variable specifying proportions for simulating binary variables.
	BETA = Varname for variable specifying the beta values for predictors from the existing model.
	MEAN = Varname for variable specifying the means for simulating continuous variables.
	SD = Varname for variable specifying the SDs for simulating continuous variables.
	TOLERANCE = Tolerance acceptable to achieve the target prevalence.
*/

syntax  , PREValence(real) ALPHA(real)  ///
			[STandardised ///
			SIMOBS(int 500000) ///
			SEED(int 123456) ///
			PREDictors(varname max=1 string) ///
			PROPortion(varname max=1 numeric) ///
			BETA(varname max=1 numeric) ///
			Mean(varname max=1 numeric) ///
			SD(varname max=1 numeric) ///
			TOLerance(real 0.004) ///
			DOTS]


*********************************************** SETUP/CHECKS
*SET UP TEMPs
tempvar 
tempname 

// check for user varnames, if none, use defaults
if "`predictors'"=="" {
	local predictors="predictors"
}

if "`proportion'"=="" {
	local proportion="prop"
}

if "`beta'"=="" {
	local beta="beta"
}

if "`mean'"=="" {
	local mean="mean"
}

if "`sd'"=="" {
	local sd="sd"
}

*******************************
// set target prevalence
local target_prev = `prevalence'

*******************************
// PARSE DATASET TO STORE NECESSARY INFO
set seed `seed'

qui count if `predictors'!=""
local vars = r(N)

// store matrix of inputs 
mkmat `beta' `proportion' `mean' `sd' , mat(input_pmsimprev) rown(`predictors')

// now set sim size
qui set obs `simobs'


forvalues a=1/`vars' {
	local b = `a'+1
	local varname`b' = `predictors'[`a']
	if `proportion'[`a']!=. {
		qui gen `varname`b'' = rbinomial(1,`proportion'[`a'])
	}
	if `mean'[`a']!=. {
		qui gen `varname`b'' = rnormal(`mean'[`a'],`sd'[`a'])
	}
	
	if "`standardised'"=="standardised" {
		qui su `varname`b''
		qui replace `varname`b'' = (`varname`b'' - r(mean))/r(sd)
	}
	
	if `a'==1 {
		qui gen lp = `beta'[`a']*`varname`b''
		local beta`a' : di %4.3f `beta'[`a']
		local PI = "`PI' + (`beta`a''*`varname`b'')"
	}
	else {
		qui replace lp = lp + `beta'[`a']*`varname`b''
		local beta`a' : di %4.3f `beta'[`a']
		local PI = "`PI' + (`beta`a''*`varname`b'')"
	}
	
	local predictor_list = "`predictor_list' `varname`b''"
	}


// qui gen p_true = invlogit(lp)

// drop `predictors' `beta' `proportion' `mean' `sd' 
********************************


**********************************
// bisection approach

// set up shells
qui gen intercept = .
qui gen p = .
qui gen outcome = .

// Define endpoints of interval to bisect.
local int_low =  `alpha' - abs(`alpha'*5)
local int_high =  `alpha' + abs(`alpha'*5)

// local iter = 1

// Initial value of the empirical value of outcome prevalence.
// This can be any value as long as it differs from the target value.
local emp_prev = 1 
if `emp_prev'==`target_prev' {
	local emp_prev = `emp_prev'-0.1
}

local iter = 0
while (abs(`emp_prev' - `target_prev') > `tolerance') {
  local iter = `iter'+1 
  local last_diff = abs(`emp_prev' - `target_prev')
  
  set seed `iter'  
  
  local int_mid = (`int_low' + `int_high')/2
  
  qui replace intercept = `int_mid'
  qui gen emp_lp = lp + intercept
  qui replace p = invlogit(emp_lp)
  qui replace outcome = rbinomial(1,p)
  
  qui su outcome 
  local emp_prev = r(mean)
  
  local current_diff = abs(`emp_prev' - `target_prev')
	
	if "`dots'"=="" {
	if `iter'==1 {
		mat output_pmsimprev = (`iter' , `target_prev', `int_low', `int_mid', `int_high', `emp_prev')
		mat colnames output_pmsimprev = Iteration Target_prev Lower_bound Intercept Upper_bound Emp_prev 
		mat rownames output_pmsimprev = `iter'
			
		matlist output_pmsimprev,  aligncolnames(r) name(c) cspec(o2& w12 | w12 | w12 | w12 | w12 | w12 o2&) rspec(&-&)
	} 
	else {
		mat iter_`iter' = (`iter' , `target_prev', `int_low', `int_mid', `int_high', `emp_prev')
		mat colnames iter_`iter' = Iteration Target_prev Lower_bound Intercept Upper_bound Emp_prev 
		mat rownames iter_`iter' = `iter'
		mat rowjoin output_pmsimprev = output_pmsimprev iter_`iter'
	
		matlist iter_`iter',  aligncolnames(r) name(n) cspec(o2& w12 | w12 | w12 | w12 | w12 | w12 o2&) rspec(-&)
	}
	}
	else {
		if `iter'==1 {
		mat output_pmsimprev = (`iter' , `target_prev', `int_low', `int_mid', `int_high', `emp_prev')
		mat colnames output_pmsimprev = Iteration Target_prev Lower_bound Intercept Upper_bound Emp_prev 
		mat rownames output_pmsimprev = `iter'
			
		di _n "Prevalence: " _c
		nois _dots `iter' 0
	} 
	else {
		mat iter_`iter' = (`iter' , `target_prev', `int_low', `int_mid', `int_high', `emp_prev')
		mat colnames iter_`iter' = Iteration Target_prev Lower_bound Intercept Upper_bound Emp_prev 
		mat rownames iter_`iter' = `iter'
		mat rowjoin output_pmsimprev = output_pmsimprev iter_`iter'
	
		nois _dots `iter' 0
	}
	}

if `emp_prev'!=. {
  if (`emp_prev' < `target_prev') {
  	local int_low = `int_mid' 
  }
  else {
  	local int_high = `int_mid'
  }
}
else {
	local int_low = `int_low'/2
	local int_high = `int_high'/2
}
  
  qui drop emp_lp
}	// end of while loop 

if "`dots'"=="dots" {
	di _c "	`iter'"
}

// Intercept for the regression model to produce the desired outcome prevalence
local intercept = `int_mid'
// di as res _n _n "Final model for simulation "
// di as res _n "Outcome = `intercept' `PI'"

**********************************
// Return code
return scalar intercept = `intercept'
return scalar empirical_prevalence = `emp_prev'
return scalar target_prevalence = `target_prev'
return scalar total_iterations = `iter'
ret sca simobs = `simobs'
ret sca tolerance = `tolerance'
return local predictor_list `"`predictor_list'"'

return mat pmsimprev_output = output_pmsimprev 
return mat pmsimprev_input = input_pmsimprev

qui keep `predictors' `beta' `proportion' `mean' `sd' 
qui drop if _n>`vars'
**********
end




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
keep if `touse'

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

		//local rown "`rown' `r'"
		}
		mat colnames `rmat' = Obs Estimate SE Lower_CI Upper_CI
		mat rownames `rmat' = "C-Statistic" //`rown'
/*
	local res cstat logit_c norm_c newcombe_c hanley_c

	tempname rmat
	matrix `rmat' = J(5,5,.)
	local i=0
	foreach r of local res {
		local ++i
		matrix `rmat'[`i',1] = `obs'
		matrix `rmat'[`i',2] = ``r''
		matrix `rmat'[`i',3] = ``r'_se'
		matrix `rmat'[`i',4] = ``r'_lb'
		matrix `rmat'[`i',5] = ``r'_ub'

		//local rown "`rown' `r'"
		}
		mat colnames `rmat' = Obs Estimate SE Lower_CI Upper_CI
		mat rownames `rmat' = "C-Statistic" "Logit(C)" "Normal" "Newcombe" "Hanley" //`rown'
 */
		
// print matrix 
if "`matrix'"!="" {
			matrix `matrix' = `rmat'
			
			//return matrix `matrix' = `rmat' 
			if "`print'"!="noprint" {
				//di as res _n "Discrimination statistics ..."
				matlist `matrix', border(all) //format(%9.3f)
							
				}
				*return matrix `matrix' = `rmat'
			}
			else { 
				if "`print'"!="noprint" {
					//di as res _n "Discrimination statistics ..."
					matlist `rmat', border(all) //format(%9.3f)
							
					}
				*return matrix rmat = `rmat'
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
  
  
program define pmsimcstat, rclass

/* Syntax
	CSTATISTIC = Target c-statistic required 
	ALPHA = Intercept of published/existing model. 
	UPPERBOUND = Upper bound for search for scalar.
	LOWERBOUND = Lower bound for search for scalar.
	SIMOBS = Number of observations to simulate super population. 
	SEED = Random number seed.
	PREDICTORS = Varname for variable specifying names of predictors in existing model.
	PROPORTION = Varname for variable specifying proportions for simulating binary variables.
	BETA = Varname for variable specifying the beta values for predictors from the existing model.
	MEAN = Varname for variable specifying the means for simulating continuous variables.
	SD = Varname for variable specifying the SDs for simulating continuous variables.
	TOLERANCE = Tolerance acceptable to achieve the target prevalence.
*/

syntax  ,  CSTATistic(real) ALPHA(real) ///
			[STandardised ///
			LOWERbound(real 0) ///
			UPPERbound(real 5) ///			
			SIMOBS(int 500000) ///
			SEED(int 123456) ///
			PREDictors(varname max=1 string) ///
			PROPortion(varname max=1 numeric) ///
			BETA(varname max=1 numeric) ///
			Mean(varname max=1 numeric) ///
			SD(varname max=1 numeric) ///
			TOLerance(real 0.004) ///
			DOTS]


*********************************************** SETUP/CHECKS
*SET UP TEMPs
tempvar 
tempname 

// check for user varnames, if none, use defaults
if "`predictors'"=="" {
	local predictors="predictors"
}

if "`proportion'"=="" {
	local proportion="prop"
}

if "`beta'"=="" {
	local beta="beta"
}

if "`mean'"=="" {
	local mean="mean"
}

if "`sd'"=="" {
	local sd="sd"
}

*******************************
// Set target c-statistic 
local target_cstat = `cstatistic'

*******************************
// PARSE DATASET TO STORE NECESSARY INFO
set seed `seed'

qui count if `predictors'!=""
local vars = r(N)

// store matrix of inputs 
mkmat `beta' `proportion' `mean' `sd' , mat(input_pmsimcstat) rown(`predictors')

// now set sim size
qui set obs `simobs'


forvalues a=1/`vars' {
	local b = `a'+1
	local varname`b' = `predictors'[`a']
	if `proportion'[`a']!=. {
		qui gen `varname`b'' = rbinomial(1,`proportion'[`a'])
	}
	if `mean'[`a']!=. {
		qui gen `varname`b'' = rnormal(`mean'[`a'],`sd'[`a'])
	}
	
	if "`standardised'"=="standardised" {
		qui su `varname`b''
		qui replace `varname`b'' = (`varname`b'' - r(mean))/r(sd)
	}
	
	if `a'==1 {
		qui gen lp = `beta'[`a']*`varname`b''
		local beta`a' : di %4.3f `beta'[`a']
		local PI = "`PI' + (`beta`a''*`varname`b'')"
	}
	else {
		qui replace lp = lp + `beta'[`a']*`varname`b''
		local beta`a' : di %4.3f `beta'[`a']
		local PI = "`PI' + (`beta`a''*`varname`b'')"
	}
	
	local predictor_list = "`predictor_list' `varname`b''"
	}


// qui gen p_true = invlogit(lp)

// drop `predictors' `beta' `proportion' `mean' `sd' 
********************************


**********************************
// bisection approach

// set up shells
qui gen cstat = .
qui gen p = .
qui gen outcome = .

// Define endpoints of interval to bisect.
local scalar_low = `lowerbound'
local scalar_high = `upperbound'


// Initial value of the empirical value of outcome prevalence.
// This can be any value as long as it differs from the target value.
local empirical_cstat = 0 // so need a check here on user input
if `empirical_cstat'==`target_cstat' {
	local empirical_cstat = `empirical_cstat'+0.1
}

local iter_cstat = 0

while (abs(`empirical_cstat' - `target_cstat') > `tolerance') {
  
  local iter_cstat = `iter_cstat'+1 
  
  local scalar = (`scalar_low' + `scalar_high')/2
  
  qui gen emp_lp = `alpha'+(lp*`scalar')
  qui replace p = invlogit(emp_lp)
  qui replace outcome = rbinomial(1,p)
  
  qui logit outcome `predictor_list'
  qui predict p_dev2
  qui predict lp2, xb

  qui pmcstat lp2 outcome 
  local empirical_cstat = r(cstat)
  
  qui drop p_dev2 lp2 emp_lp
  
  if "`dots'"=="" {
  if `iter_cstat'==1 {
		mat output_pmsimcstat = (`iter_cstat' , `target_cstat', `scalar_low', `scalar', `scalar_high', `empirical_cstat')
		mat colnames output_pmsimcstat = Iteration Target_cstat Lower_bound Scalar Upper_bound Emp_cstat 
		mat rownames output_pmsimcstat = `iter_cstat'
		
		matlist output_pmsimcstat,  aligncolnames(r) name(c) cspec(o2& w12 | w12 | w12 | w12 | w12 | w12 o2&) rspec(&-&)
	} 
	else {
		mat iter_cstat_`iter_cstat' = (`iter_cstat' , `target_cstat', `scalar_low', `scalar', `scalar_high', `empirical_cstat')
		mat colnames iter_cstat_`iter_cstat' = Iteration Target_cstat Lower_bound Scalar Upper_bound Emp_cstat 
		mat rownames iter_cstat_`iter_cstat' = `iter_cstat'
		mat rowjoin output_pmsimcstat = output_pmsimcstat iter_cstat_`iter_cstat'
	
		matlist iter_cstat_`iter_cstat',  aligncolnames(r) name(n) cspec(o2& w12 | w12 | w12 | w12 | w12 | w12 o2&) rspec(-&)
	}
  }
  else {
  	if `iter_cstat'==1 {
		mat output_pmsimcstat = (`iter_cstat' , `target_cstat', `scalar_low', `scalar', `scalar_high', `empirical_cstat')
		mat colnames output_pmsimcstat = Iteration Target_cstat Lower_bound Scalar Upper_bound Emp_cstat 
		mat rownames output_pmsimcstat = `iter_cstat'
		
		di _n "C-statistic: " _c
		nois _dots `iter_cstat' 0
	} 
	else {
		mat iter_cstat_`iter_cstat' = (`iter_cstat' , `target_cstat', `scalar_low', `scalar', `scalar_high', `empirical_cstat')
		mat colnames iter_cstat_`iter_cstat' = Iteration Target_cstat Lower_bound Scalar Upper_bound Emp_cstat 
		mat rownames iter_cstat_`iter_cstat' = `iter_cstat'
		mat rowjoin output_pmsimcstat = output_pmsimcstat iter_cstat_`iter_cstat'
	
		nois _dots `iter_cstat' 0
	}
  }
  
  
  if (`empirical_cstat' < `target_cstat') {
  	local scalar_low = `scalar' 
  }
  else {
  	local scalar_high = `scalar'
  }
    

}	// end of while loop

if "`dots'"=="dots" {
	di _c "	`iter_cstat'" _n
}


forvalues a=1/`vars' {
	local b = `a'+1
	
	local scaled_beta`a' : di %4.3f `beta`a''*`scalar'
	local scaled_PI = "`scaled_PI' + (`scaled_beta`a''*`varname`b'')"
	
	}

// Scalar for the regression model to produce the desired c-stat.
local cstat = `empirical_cstat'
// di as res _n _n "Final model for simulation "
// di as res _n "Outcome = `alpha' `scaled_PI'"
// di as res _n "Resulting c-statistic = `cstat'"

**********************************
// Return code
return scalar intercept = `alpha'
ret sca scalar = `scalar'
ret sca simobs = `simobs'
ret sca tolerance = `tolerance'
return scalar empirical_cstat = `empirical_cstat'
return scalar target_prevalence = `target_cstat'
return scalar total_iterations = `iter_cstat'
return local predictor_list `"`predictor_list'"'
	
return mat pmsimcstat_output = output_pmsimcstat
return mat pmsimcstat_input = input_pmsimcstat

qui keep `predictors' `beta' `proportion' `mean' `sd' 
qui drop if _n>`vars'
**********
end



program define pmsimcheck, rclass

/* Syntax
	PREVALENCE = Target prevalence required 
	ALPHA = Intercept of published/existing model. Or start point for intercept.
	SIMOBS = Number of observations to simulate super population. 
	SEED = Random number seed.
	PREDICTORS = Varname for variable specifying names of predictors in existing model.
	PROPORTION = Varname for variable specifying proportions for simulating binary variables.
	BETA = Varname for variable specifying the beta values for predictors from the existing model.
	MEAN = Varname for variable specifying the means for simulating continuous variables.
	SD = Varname for variable specifying the SDs for simulating continuous variables.
	TOLERANCE = Tolerance acceptable to achieve the target prevalence.
*/

syntax  , PREValence(real) ALPHA(real) CSTATistic(real) ///
			[STandardised ///
			SIMOBS(int 500000) ///
			SEED(int 123456) ///
			PREDictors(varname max=1 string) ///
			PROPortion(varname max=1 numeric) ///
			BETA(varname max=1 numeric) ///
			Mean(varname max=1 numeric) ///
			SD(varname max=1 numeric) ///
			TOLerance(real 0.004) KEEP]


*********************************************** SETUP/CHECKS
*SET UP TEMPs
tempvar 
tempname 

// check for user varnames, if none, use defaults
if "`predictors'"=="" {
	local predictors="predictors"
}

if "`proportion'"=="" {
	local proportion="prop"
}

if "`beta'"=="" {
	local beta="beta"
}

if "`mean'"=="" {
	local mean="mean"
}

if "`sd'"=="" {
	local sd="sd"
}

*******************************
// set target prevalence
local target_prev = `prevalence'
local target_cstat = `cstatistic'

*******************************
// PARSE DATASET TO STORE NECESSARY INFO
set seed `seed'

qui count if `predictors'!=""
local vars = r(N)

// now set sim size
qui set obs `simobs'


forvalues a=1/`vars' {
	local b = `a'+1
	local varname`b' = `predictors'[`a']
	if `proportion'[`a']!=. {
		qui gen `varname`b'' = rbinomial(1,`proportion'[`a'])
	}
	if `mean'[`a']!=. {
		qui gen `varname`b'' = rnormal(`mean'[`a'],`sd'[`a'])
	}
	
	if "`standardised'"=="standardised" {
		qui su `varname`b''
		qui replace `varname`b'' = (`varname`b'' - r(mean))/r(sd)
	}
	
	if `a'==1 {
		qui gen lp = `beta'[`a']*`varname`b''
		local beta`a' : di %4.3f `beta'[`a']
		local PI = "`PI' + (`beta`a''*`varname`b'')"
	}
	else {
		qui replace lp = lp + `beta'[`a']*`varname`b''
		local beta`a' : di %4.3f `beta'[`a']
		local PI = "`PI' + (`beta`a''*`varname`b'')"
	}
	
	local predictor_list = "`predictor_list' `varname`b''"
	}

qui replace lp = lp + `alpha'


// drop `predictors' `beta' `proportion' `mean' `sd' 
********************************

qui gen p_true = invlogit(lp)
qui gen outcome = rbinomial(1, p_true )
qui su `predictor_list' lp p_true outcome
local emp_prev = r(mean)

qui logit outcome `predictor_list'
qui predict p_dev2
qui predict lp2, xb

qui pmcstat lp2 outcome 
local empirical_cstat = r(cstat)

if "`keep'"==""{
qui drop p_dev2 lp2 outcome p_true

qui keep `predictors' `beta' `proportion' `mean' `sd' 
qui drop if _n>`vars'
}
else {
	qui drop p_dev2 lp2
}
**********************************
// Return code
return scalar intercept = `alpha'
// ret sca scalar = `scalar'
// ret sca simobs = `simobs'
// ret sca tolerance = `tolerance'
return scalar empirical_prevalence = `emp_prev'
return scalar empirical_cstat = `empirical_cstat'
// return scalar target_prevalence = `target_cstat'
// return scalar total_iterations = `iter_cstat'
// return local predictor_list `"`predictor_list'"'

end



program define pmsim, rclass

/* Syntax
	PREVALENCE = Target prevalence required 
	ALPHA = Intercept of published/existing model. Or start point for intercept.
	SIMOBS = Number of observations to simulate super population. 
	SEED = Random number seed.
	PREDICTORS = Varname for variable specifying names of predictors in existing model.
	PROPORTION = Varname for variable specifying proportions for simulating binary variables.
	BETA = Varname for variable specifying the beta values for predictors from the existing model.
	MEAN = Varname for variable specifying the means for simulating continuous variables.
	SD = Varname for variable specifying the SDs for simulating continuous variables.
	TOLERANCE = Tolerance acceptable to achieve the target prevalence.
*/

syntax  , PREValence(real) ALPHA(real) CSTATistic(real) ///
			[STandardised ///
			LOWERbound(real 0) ///
			UPPERbound(real 5) ///	
			SIMOBS(int 500000) ///
			SEED(int 123456) ///
			PREDictors(varname max=1 string) ///
			PROPortion(varname max=1 numeric) ///
			BETA(varname max=1 numeric) ///
			Mean(varname max=1 numeric) ///
			SD(varname max=1 numeric) ///
			TOLerance(real 0.005)]

// local tolerance = 0.001
// local prevalence = 0.052
// local alpha = -3.81
// local cstatistic = 0.83
// local simobs = 500000

*********************************************** SETUP/CHECKS
*SET UP TEMPs
tempvar 
tempname 

// check for user varnames, if none, use defaults
if "`predictors'"=="" {
	local predictors="predictors"
}

if "`proportion'"=="" {
	local proportion="prop"
}

if "`beta'"=="" {
	local beta="beta"
}

if "`mean'"=="" {
	local mean="mean"
}

if "`sd'"=="" {
	local sd="sd"
}

*******************************
// set targets 
local target_prev = `prevalence'
local target_cstat = `cstatistic'

// set achieve 
local achieved = 0
local p_achieved = 0
local c_achieved = 0
*******************************

pmsimprev , prev(`prevalence') alpha(`alpha') simobs(`simobs') pred(`predictors') beta(`beta') prop(`proportion') mean(`mean') sd(`sd')  tol(`tolerance') seed(`seed') dots standardised

local final_alpha = r(intercept)
local iter_prev = r(total_iterations)
mat input_pmsim = r(pmsimprev_input)

pmsimcstat , alpha(`final_alpha') cstat(`cstatistic')  simobs(`simobs') pred(`predictors') beta(`beta') prop(`proportion') mean(`mean') sd(`sd')  tol(`tolerance') seed(`seed') lowerbound(`lowerbound') upperbound(`upperbound') dots standardised

local iter_cstat = r(total_iterations)
local final_scalar = r(scalar)
local current_scalar = r(scalar)
qui replace beta = beta * `final_scalar'

pmsimcheck , prev(`prevalence') alpha(`final_alpha') cstat(`cstatistic')  simobs(`simobs') pred(`predictors') beta(`beta') prop(`proportion') mean(`mean') sd(`sd')  tol(`tolerance') seed(`seed') standardised

local current_prev = r(empirical_prevalence)
local current_cstat = r(empirical_cstat)

if (abs(`current_prev' - `target_prev') < `tolerance') { 
// 	di as res _n "prev  met"
	
	local p_achieved = 1
}
	
if (abs(`current_cstat' - `target_cstat') < `tolerance') {
// 		di as res _n "c-stat  met"
		
		local c_achieved = 1
}

if `p_achieved'==1 & `c_achieved'==1 {
// 	di as res _n "both met"
	local achieved = 1
}

local pmsim_iter = 1
mat output_pmsim = (`pmsim_iter', `final_alpha', `target_prev', `current_prev', `final_scalar', `target_cstat', `current_cstat', `simobs', `iter_prev', `iter_cstat')
mat colnames output_pmsim = Iteration Intercept Target_prev Emp_prev Scalar Target_cstat Emp_cstat Obs Iter_prev Iter_cstat
mat rownames output_pmsim = `pmsim_iter'
			
matlist output_pmsim,  aligncolnames(r) name(c) cspec(o2& w12 | w12 | w12 | w12 | w12 | w12 | w12 | w12 | w12 | w12 o2&) rspec(&-&)

********************************
// While loop 
while `achieved'==0 {
	local pmsim_iter = `pmsim_iter'+1
	
	if `p_achieved'==0 {
// 		di "do prev sim again"
	pmsimprev , prev(`prevalence') alpha(`final_alpha')  simobs(`simobs') pred(`predictors') beta(`beta') prop(`proportion') mean(`mean') sd(`sd')  tol(`tolerance') seed(`seed') dots standardised

	local final_alpha = r(intercept)
	local iter_prev = r(total_iterations)
	
	local p_achieved = 0
	local c_achieved = 0

	pmsimcheck , prev(`prevalence') alpha(`final_alpha') cstat(`cstatistic')  simobs(`simobs') pred(`predictors') beta(`beta') prop(`proportion') mean(`mean') sd(`sd')  tol(`tolerance') seed(`seed') standardised

	local current_prev = r(empirical_prevalence)
	local current_cstat = r(empirical_cstat)
if (abs(`current_prev' - `target_prev') < `tolerance') { 
// 	di as res _n "prev  met"
	
	local p_achieved = 1
}

if (abs(`current_cstat' - `target_cstat') < `tolerance') {
// 	di as res _n "c-stat  met"
	
	local c_achieved = 1
}

if `p_achieved'==1 & `c_achieved'==1 {
// 	di as res _n "both met"
	local achieved = 1
	
	mat iter_pmsim_`pmsim_iter' = (`pmsim_iter', `final_alpha', `target_prev', `current_prev', `current_scalar', `target_cstat', `current_cstat', `simobs', `iter_prev', `iter_cstat')
mat colnames iter_pmsim_`pmsim_iter' = Iteration Intercept Target_prev Emp_prev Scalar Target_cstat Emp_cstat Obs Iter_prev Iter_cstat
mat rownames iter_pmsim_`pmsim_iter' = `pmsim_iter'
mat rowjoin output_pmsim = output_pmsim iter_pmsim_`pmsim_iter'

matlist iter_pmsim_`pmsim_iter',  aligncolnames(r) name(c) cspec(o2& w12 | w12 | w12 | w12 | w12 | w12 | w12 | w12 | w12 | w12 o2&) rspec(&-&)

// exit
continue, break 
}
}

if `c_achieved'==0 {
	
// 	di _n "do cstat sim again"
	pmsimcstat , alpha(`final_alpha') cstat(`cstatistic')  simobs(`simobs') pred(`predictors') beta(`beta') prop(`proportion') mean(`mean') sd(`sd')  tol(`tolerance') seed(`seed') dots standardised

	local iter_cstat = r(total_iterations)
	local final_scalar = `final_scalar'*r(scalar)
	local current_scalar = r(scalar)
	qui replace beta = beta * `current_scalar'

	local p_achieved = 0
	local c_achieved = 0
	
	pmsimcheck , prev(`prevalence') alpha(`final_alpha') cstat(`cstatistic')  simobs(`simobs') pred(`predictors') beta(`beta') prop(`proportion') mean(`mean') sd(`sd')  tol(`tolerance') seed(`seed') standardised

	local current_prev = r(empirical_prevalence)
	local current_cstat = r(empirical_cstat)
if (abs(`current_prev' - `target_prev') < `tolerance') { 
// 	di as res _n "prev  met"
	
	local p_achieved = 1
}

if (abs(`current_cstat' - `target_cstat') < `tolerance') {
// 	di as res _n "c-stat  met"
	
	local c_achieved = 1
}

if `p_achieved'==1 & `c_achieved'==1 {
// 	di as res _n "both met"
	local achieved = 1
	
	mat iter_pmsim_`pmsim_iter' = (`pmsim_iter', `final_alpha', `target_prev', `current_prev', `current_scalar', `target_cstat', `current_cstat', `simobs', `iter_prev', `iter_cstat')
mat colnames iter_pmsim_`pmsim_iter' = Iteration Intercept Target_prev Emp_prev Scalar Target_cstat Emp_cstat Obs Iter_prev Iter_cstat
mat rownames iter_pmsim_`pmsim_iter' = `pmsim_iter'
mat rowjoin output_pmsim = output_pmsim iter_pmsim_`pmsim_iter'

matlist iter_pmsim_`pmsim_iter',  aligncolnames(r) name(c) cspec(o2& w12 | w12 | w12 | w12 | w12 | w12 | w12 | w12 | w12 | w12 o2&) rspec(&-&)

// exit
continue, break 
}
}
	
mat iter_pmsim_`pmsim_iter' = (`pmsim_iter', `final_alpha', `target_prev', `current_prev', `current_scalar', `target_cstat', `current_cstat', `simobs', `iter_prev', `iter_cstat')
mat colnames iter_pmsim_`pmsim_iter' = Iteration Intercept Target_prev Emp_prev Scalar Target_cstat Emp_cstat Obs Iter_prev Iter_cstat
mat rownames iter_pmsim_`pmsim_iter' = `pmsim_iter'
mat rowjoin output_pmsim = output_pmsim iter_pmsim_`pmsim_iter'

matlist iter_pmsim_`pmsim_iter',  aligncolnames(r) name(c) cspec(o2& w12 | w12 | w12 | w12 | w12 | w12 | w12 | w12 | w12 | w12 o2&) rspec(&-&)

// matlist output_pmsim,  aligncolnames(r) name(n) cspec(o2& w12 | w12 | w12 | w12 | w12 | w12 | w12 | w12 | w12 | w12 o2&) rspec(-&)
	
	
}

// Intercept for the regression model to produce the desired outcome prevalence
di as res _n _n "Final model for simulation "
di as res _n "Intercept = `final_alpha'"
di as res _n "Achieved empirical prevalence = `r(empirical_prevalence)'"
di as res _n "Scalar = `final_scalar'"
di as res _n "Achieved empirical c-statistic = `r(empirical_cstat)'"

************************
// Return code
return scalar intercept = `final_alpha'
ret sca scalar = `final_scalar'
ret sca simobs = `simobs'
ret sca tolerance = `tolerance'
return scalar empirical_prevalence = `current_prev'
return scalar empirical_cstat = `current_cstat'
return scalar target_prevalence = `target_prev'
return scalar target_cstat = `target_cstat'
return scalar total_iterations = `pmsim_iter'
// return local predictor_list `"`predictor_list'"'
		
return mat pmsim_output = output_pmsim
return mat pmsim_input = input_pmsim

************************
// Leave behind the simulated data
// pmsimcheck , prev(`prevalence') alpha(`final_alpha') cstat(`cstatistic')  simobs(`simobs') pred(`predictors') beta(`beta') prop(`proportion') mean(`mean') sd(`sd')  tol(`tolerance') seed(`seed') keep

end

********************

