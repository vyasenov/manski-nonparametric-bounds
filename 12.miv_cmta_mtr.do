/*=========================================*/
/* 12. PROGRAM FOR MIV + cMTA + MTR BOUNDS */
/*=========================================*/
capture program drop miv_mtr_mta_bounds
program miv_mtr_mta_bounds, rclass 

	version 14
	args y d z mtrsign MTAsign
	
	sum `d', meanonly
	local max_d = r(max)
	local max_d_1 = r(max) -1
	sum `z', meanonly
	local max_z = r(max)
	sum `y', meanonly
	local N = r(N)

	forvalues k = 1 / `max_z' {
		tempname n_z`k'
		qui count if `z' == `k' 
		gen `n_z`k'' = r(N)
	}
	
	forvalues i = 1 / `max_d' {		// auxiliary
		forvalues k = 1 / `max_z' {
			tempname lb_d`i'z`k' ub_d`i'z`k'
			
			sum `y' if `d' > `i' & `z' == `k', meanonly		// <---------
			local meang = r(mean) 
			local cond_probg = r(N) / `n_z`k''		// <---------

			sum `y' if `d' == `i' & `z' == `k', meanonly	// <---------
			local mean = r(mean)
			
			sum `y' if `d' < `i' & `z' == `k', meanonly		// <---------
			local meanl = r(mean)
			local cond_probl = r(N) / `n_z`k''		// <---------

			if `mtrsign' == 1 & `MTAsign' == 1 {
				if `i' == 1 { 		// correction for the upper bound of the lowest treatment level
					gen `lb_d`i'z`k'' = `meang' * `cond_probg' + `mean' * (1 - `cond_probg')
					gen `ub_d`i'z`k'' =		          0 	   + `mean'
				}
				
				else if `i' == `max_d' { 		// correction for the lower bound of the highest treatment level
					gen `lb_d`i'z`k'' = 			  0        + `mean'
					gen `ub_d`i'z`k'' = `meanl' * `cond_probl' + `mean' * (1 - `cond_probl')
				}
				
				else {
					gen `lb_d`i'z`k'' = `meang' * `cond_probg' + `mean' * (1 - `cond_probg')
					gen `ub_d`i'z`k'' = `meanl' * `cond_probl' + `mean' * (1 - `cond_probl')
				}
			}
			
			else if `mtrsign' == 2 & `MTAsign' == 2 {
				if `i' == 1 { 		// correction for the lower bound of the lowest treatment level
					gen `lb_d`i'z`k'' = 	          0 	   + `mean' 
					gen `ub_d`i'z`k'' =	`meang' * `cond_probg' + `mean' * (1 - `cond_probg')
				}
				
				else if `i' == `max_d' { 		// correction for the upper bound of the highest treatment level
					gen `lb_d`i'z`k'' = `meanl' * `cond_probl' + `mean' * (1 - `cond_probl')
					gen `ub_d`i'z`k'' = 			  0        + `mean'
				}
				
				else {
					gen `lb_d`i'z`k'' = `meanl' * `cond_probl' + `mean' * (1 - `cond_probl')
					gen `ub_d`i'z`k'' = `meang' * `cond_probg' + `mean' * (1 - `cond_probg')
				}	
			}
			
				else {
					gen `lb_d`i'z`k'' = .
					gen `ub_d`i'z`k'' = .			
			}
		}
	}	
	
* this part is not fully automated
	forvalues i = 1 / `max_d' {	 // auxiliary
		tempname lb_max_d`i'z1 lb_max_d`i'z2 lb_max_d`i'z3 lb_max_d`i'z4 lb_max_d`i'z5 lb_max_d`i'z6 lb_max_d`i'z7 lb_max_d`i'z8 lb_max_d`i'z9 lb_max_d`i'z10
		tempname ub_min_d`i'z1 ub_min_d`i'z2 ub_min_d`i'z3 ub_min_d`i'z4 ub_min_d`i'z5 ub_min_d`i'z6 ub_min_d`i'z7 ub_min_d`i'z8 ub_min_d`i'z9 ub_min_d`i'z10
			
		gen `lb_max_d`i'z1' = `lb_d`i'z1'
		gen `lb_max_d`i'z2' = max(`lb_d`i'z2', `lb_max_d`i'z1')
		gen `lb_max_d`i'z3' = max(`lb_d`i'z3', `lb_max_d`i'z2')
		gen `lb_max_d`i'z4' = max(`lb_d`i'z4', `lb_max_d`i'z3')
		gen `lb_max_d`i'z5' = max(`lb_d`i'z5', `lb_max_d`i'z4')
		gen `lb_max_d`i'z6' = max(`lb_d`i'z6', `lb_max_d`i'z5')
		gen `lb_max_d`i'z7' = max(`lb_d`i'z7', `lb_max_d`i'z6')
		gen `lb_max_d`i'z8' = max(`lb_d`i'z8', `lb_max_d`i'z7')
		gen `lb_max_d`i'z9' = max(`lb_d`i'z9', `lb_max_d`i'z8')
		gen `lb_max_d`i'z10' = max(`lb_d`i'z10', `lb_max_d`i'z9')
			
		gen `ub_min_d`i'z10' = `ub_d`i'z10'
		gen `ub_min_d`i'z9' = min(`ub_d`i'z9', `ub_min_d`i'z10')
		gen `ub_min_d`i'z8' = min(`ub_d`i'z8', `ub_min_d`i'z9')
		gen `ub_min_d`i'z7' = min(`ub_d`i'z7', `ub_min_d`i'z8')
		gen `ub_min_d`i'z6' = min(`ub_d`i'z6', `ub_min_d`i'z7')
		gen `ub_min_d`i'z5' = min(`ub_d`i'z5', `ub_min_d`i'z6')
		gen `ub_min_d`i'z4' = min(`ub_d`i'z4', `ub_min_d`i'z5')
		gen `ub_min_d`i'z3' = min(`ub_d`i'z3', `ub_min_d`i'z4')
		gen `ub_min_d`i'z2' = min(`ub_d`i'z2', `ub_min_d`i'z3')
		gen `ub_min_d`i'z1' = min(`ub_d`i'z1', `ub_min_d`i'z2')		
	}
			
* this part is not fully automated
	forvalues i = 1 / `max_d' {  // mean potential outcomes
		return scalar lb_Y`i' =  round((`n_z1' / `N') * `lb_max_d`i'z1'  + ///
			(`n_z2' / `N') * `lb_max_d`i'z2' + ///
			(`n_z3' / `N') * `lb_max_d`i'z3' + ///
			(`n_z4' / `N') * `lb_max_d`i'z4' + ///
			(`n_z5' / `N') * `lb_max_d`i'z5' + ///
			(`n_z6' / `N') * `lb_max_d`i'z6' + ///
			(`n_z7' / `N') * `lb_max_d`i'z7' + ///
			(`n_z8' / `N') * `lb_max_d`i'z8' + ///
			(`n_z9' / `N') * `lb_max_d`i'z9' + ///
			(`n_z10' / `N') * `lb_max_d`i'z10', 0.001)
			
		return scalar ub_Y`i' =  round((`n_z1' / `N') * `ub_min_d`i'z1'  + ///
			(`n_z2' / `N') * `ub_min_d`i'z2' + ///
			(`n_z3' / `N') * `ub_min_d`i'z3' + ///
			(`n_z4' / `N') * `ub_min_d`i'z4' + ///
			(`n_z5' / `N') * `ub_min_d`i'z5' + ///
			(`n_z6' / `N') * `ub_min_d`i'z6' + ///
			(`n_z7' / `N') * `ub_min_d`i'z7' + ///
			(`n_z8' / `N') * `ub_min_d`i'z8' + ///
			(`n_z9' / `N') * `ub_min_d`i'z9' + ///
			(`n_z10' / `N') * `ub_min_d`i'z10', 0.001)
			
		di "MIV+MTR+MTA Y`i' = (`return(lb_Y`i')' , `return(ub_Y`i')')"
	}

	forvalues j = 1 / `max_d_1' {		// treatment effects
		if `mtrsign' == 1 & `MTAsign' == 1 { // make sure the upper bound of MTR- sign is non-positive
			return scalar lb_ATE_`max_d'`j' =  round(return(lb_Y`max_d') - return(ub_Y`j'), 0.001)
			return scalar ub_ATE_`max_d'`j' =  round(min(return(ub_Y`max_d') - return(lb_Y`j'), 0), 0.001)	
		}

		else if `mtrsign' == 2 & `MTAsign' == 2 { // make sure the upper bound of MTR- sign is non-positive
			return scalar lb_ATE_`max_d'`j' =  round(max(return(lb_Y`max_d') - return(ub_Y`j'), 0), 0.001)
			return scalar ub_ATE_`max_d'`j' =  round(return(ub_Y`max_d') - return(lb_Y`j'), 0.001)	
		}
		
		else {
			return scalar lb_ATE_`max_d'`j' =  .
			return scalar ub_ATE_`max_d'`j' =  .		
		}
		
		di "MIV+MTR+MTA ATE `max_d'-`j' = (`return(lb_ATE_`max_d'`j')' , `return(ub_ATE_`max_d'`j')')"
	}
end