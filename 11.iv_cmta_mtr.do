/*========================================*/
/* 11. PROGRAM FOR IV + cMTA + MTR BOUNDS */
/*========================================*/
capture program drop iv_mtr_mta_bounds
program iv_mtr_mta_bounds, rclass 

	version 14
	args y d z mtrsign MTAsign
	
	sum `d', meanonly
	local max_d = r(max)
	local max_d_1 = r(max) - 1
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
					gen lb_d`i'z`k' = `meang' * `cond_probg' + `mean' * (1 - `cond_probg')
					gen ub_d`i'z`k' =		          0 	 + `mean'
				}
				
				else if `i' == `max_d' {   // correction for the lower bound of the highest treatment level
					gen lb_d`i'z`k' = 			  0          + `mean'
					gen ub_d`i'z`k' = `meanl' * `cond_probl' + `mean' * (1 - `cond_probl')
				}
				
				else {
					gen lb_d`i'z`k' = `meang' * `cond_probg' + `mean' * (1 - `cond_probg')
					gen ub_d`i'z`k' = `meanl' * `cond_probl' + `mean' * (1 - `cond_probl')
				}
			}
			
			else if `mtrsign' == 2 & `MTAsign' == 2 {
				if `i' == 1 { 	// correction for the lower bound of the lowest treatment level
					gen lb_d`i'z`k' = 	          0 	       + `mean' 
					gen ub_d`i'z`k' =	`meang' * `cond_probg' + `mean' * (1 - `cond_probg')
				}
				
				else if `i' == `max_d' { 	// correction for the upper bound of the highest treatment level
					gen lb_d`i'z`k' = `meanl' * `cond_probl' + `mean' * (1 - `cond_probl')
					gen ub_d`i'z`k' = 			  0          + `mean'
				}
				
				else {
					gen lb_d`i'z`k' = `meanl' * `cond_probl' + `mean' * (1 - `cond_probl')
					gen ub_d`i'z`k' = `meang' * `cond_probg' + `mean' * (1 - `cond_probg')
				}	
			}
			
			else {
				gen lb_d`i'z`k' = .
				gen ub_d`i'z`k' = .			
			}
		}
	}	

	forvalues i = 1 / `max_d' {	
		tempvar aux1 aux2

		egen `aux1' = rowmax(lb_d`i'z*)
		egen `aux2' = rowmin(ub_d`i'z*)
		
		return scalar lb_Y`i' = round(`aux1', 0.001)
		return scalar ub_Y`i' = round(`aux2', 0.001)		
		
		di "IV+MTR+MTA Y`i' = (`return(lb_Y`i')' , `return(ub_Y`i')')"
	}
	
	drop lb* ub*

	forvalues j = 1 / `max_d_1' {		// treatment effects
		if `mtrsign' == 1 & `MTAsign' == 1 { 	// make sure the upper bound of MTR- sign is non-positive
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
		
		di "IV+MTR+MTA ATE `max_d'-`j' = (`return(lb_ATE_`max_d'`j')' , `return(ub_ATE_`max_d'`j')')"
	}
end
