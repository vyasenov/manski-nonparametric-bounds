/*==================================*/
/* 10. PROGRAM FOR MTR + MTA BOUNDS */
/*==================================*/
capture program drop mtr_mta_bounds
program mtr_mta_bounds, rclass 

	version 14
	args y d mtrsign MTAsign
	
	sum `d', meanonly
	local max_d = r(max)
	local max_d_1 = r(max) - 1
	sum `y', meanonly
	local N = r(N)
	
	di "Potential Outcomes:"
	forvalues i = 1 / `max_d' {	// potential outcomes
		sum `y' if `d' > `i', meanonly		// <---------
		local meang = r(mean)
		local probg = r(N) / `N'
		
		sum `y' if `d' == `i', meanonly		// <---------
		local mean = r(mean)
		
		sum `y' if `d' < `i', meanonly
		local meanl = r(mean)
		local probl = r(N) / `N'
		
		if (`mtrsign' == 2 & `MTAsign' == 2) {     // potential outcomes MTR+, MTA+ 
			if `i' == 1 {   // correction for the lower bound of the lowest treatment level
				return scalar lb_Y`i' = round(`mean', 0.001)
				return scalar ub_Y`i' = round(`meang' * `probg' + `mean' * (1 - `probg'), 0.001)
			}
			
			else if `i' == `max_d' {   // correction for the upper bound of the highest treatment level
				return scalar lb_Y`i' = round(`meanl' * `probl' + `mean' * (1 - `probl'), 0.001)
				return scalar ub_Y`i' = round(`mean' , 0.001)
			}
			
			else {
				return scalar lb_Y`i' = round(`meanl' * `probl' + `mean' * (1 - `probl'), 0.001)
				return scalar ub_Y`i' = round(`meang' * `probg' + `mean' * (1 - `probg'), 0.001)
			}
		}
		
		else if (`mtrsign' == 1 & `MTAsign' == 1) {     // potential outcomes MTR-, MTA-
				if `i' == 1 {   // correction for the upper bound of the lowest treatment level
				return scalar lb_Y`i' = round(`meang' * `probg' + `mean' * (1 - `probg'), 0.001)
				return scalar ub_Y`i' = round(`mean' , 0.001)
			}
			
			else if `i' == `max_d' {   // correction for the lower bound of the highest treatment level
				return scalar lb_Y`i' = round(`mean' , 0.001)
				return scalar ub_Y`i' = round(`meanl' * `probl' + `mean' * (1 - `probl'), 0.001)
			}
			
			else {
				return scalar lb_Y`i' = round(`meang' * `probg' + `mean' * (1 - `probg'), 0.001)
				return scalar ub_Y`i' = round(`meanl' * `probl' + `mean' * (1 - `probl'), 0.001)
			}
		}
		
		else {
				return scalar lb_Y`i' = .
				return scalar ub_Y`i' = .
		}
		
		di "MTR+MTA Y`i' = (`return(lb_Y`i')' , `return(ub_Y`i')')"
	}

	forvalues j = 1 / `max_d_1' {  // treatment effects
		if `mtrsign' == 1 & `MTAsign' == 1 { // make sure the upper bound of MTR- sign is non-positive
				return scalar lb_ATE_`max_d'`j' = round(return(lb_Y`max_d') - return(ub_Y`j'), 0.001)
				return scalar ub_ATE_`max_d'`j' = round(min(return(ub_Y`max_d') - return(lb_Y`j'), 0), 0.001)	
				
				if `j' <= 3 {
				return scalar lb_ATE_4`j' = round(return(lb_Y4) - return(ub_Y`j'), 0.001)
				return scalar ub_ATE_4`j' = round(min(return(ub_Y4) - return(lb_Y`j'),0), 0.001)
			}	

				if `j' <= 2 {
				return scalar lb_ATE_3`j' = round(return(lb_Y3) - return(ub_Y`j'), 0.001)
				return scalar ub_ATE_3`j' = round(min(return(ub_Y3) - return(lb_Y`j'),0), 0.001)
			}
				
				return scalar lb_ATE_21 = round(return(lb_Y2) - return(ub_Y1), 0.001)
				return scalar ub_ATE_21 = round(min(return(ub_Y2) - return(lb_Y1),0), 0.001)	
			}
		
		else if `mtrsign' == 2 & `MTAsign' == 2 {  // make sure the lower bound of MTR+ sign is non-negative
			return scalar lb_ATE_`max_d'`j' = round(max(return(lb_Y`max_d') - return(ub_Y`j'), 0), 0.001)
			return scalar ub_ATE_`max_d'`j' = round(return(ub_Y`max_d') - return(lb_Y`j'), 0.001)	
			
			if `j' <= 3 {
				return scalar lb_ATE_4`j' = round(max(return(lb_Y4) - return(ub_Y`j'),0), 0.001)
				return scalar ub_ATE_4`j' = round(return(ub_Y4) - return(lb_Y`j'), 0.001)
			}	

			if `j' <= 2 {
				return scalar lb_ATE_3`j' = round(max(return(lb_Y3) - return(ub_Y`j'),0), 0.001)
				return scalar ub_ATE_3`j' = round(return(ub_Y3) - return(lb_Y`j'), 0.001)
			}	
			
			return scalar lb_ATE_21 = round(max(return(lb_Y2) - return(ub_Y1),0), 0.001)
			return scalar ub_ATE_21 = round(return(ub_Y2) - return(lb_Y1), 0.001)	
		}

		else {
			return scalar lb_ATE_`max_d'`j' = .
			return scalar ub_ATE_`max_d'`j' = .
		}
	}
	
		di " "
		di "Average Treatment Effects:"
		di "MTR+MTA ATE 5-4 = (`return(lb_ATE_54)' , `return(ub_ATE_54)')"	
		di "MTR+MTA ATE 5-3 = (`return(lb_ATE_53)' , `return(ub_ATE_53)')"
		di "MTR+MTA ATE 5-2 = (`return(lb_ATE_52)' , `return(ub_ATE_52)')"
		di "MTR+MTA ATE 5-1 = (`return(lb_ATE_51)' , `return(ub_ATE_51)')"
		di "------------"
		di "MTR+MTA ATE 4-3 = (`return(lb_ATE_43)' , `return(ub_ATE_43)')"
		di "MTR+MTA ATE 4-2 = (`return(lb_ATE_42)' , `return(ub_ATE_42)')"
		di "MTR+MTA ATE 4-1 = (`return(lb_ATE_41)' , `return(ub_ATE_41)')"
		di "------------"
		di "MTR+MTA ATE 3-2 = (`return(lb_ATE_32)' , `return(ub_ATE_32)')"
		di "MTR+MTA ATE 3-1 = (`return(lb_ATE_31)' , `return(ub_ATE_21)')"
		di "------------"
		di "MTR+MTA ATE 2-1 = (`return(lb_ATE_21)' , `return(ub_ATE_21)')"
	
	* now convert to beta's
		forvalues j = 1/5 {    // this is to get the mean of d_cont for each d
		qui sum d_continuous if d == `j', meanonly
		scalar d`j' = r(mean)
		}

	forvalues j = 1/5 {  // this is to get the multiplicaiton factors to convert to beta
		local counter = `j' - 1
		while `counter' > 0 {
			scalar d`j'`counter' = d`j' - d`counter'
			scalar factor_d`j'`counter' = 1 / d`j'`counter'
			local counter = `counter' - 1
		}
	}
		di " "
		di "Betas:"
		di "MTR+MTA Beta 5-4 = (" round(factor_d54 * `return(lb_ATE_54)', .001) "," round(factor_d54 * `return(ub_ATE_54)', .001) ")"	
		di "MTR+MTA Beta 5-3 = (" round(factor_d53 * `return(lb_ATE_53)', .001) "," round(factor_d53 * `return(ub_ATE_53)', .001) ")"
		di "MTR+MTA Beta 5-2 = (" round(factor_d52 * `return(lb_ATE_52)', .001) "," round(factor_d52 * `return(ub_ATE_52)', .001) ")"
		di "MTR+MTA Beta 5-1 = (" round(factor_d51 * `return(lb_ATE_51)', .001) "," round(factor_d51 * `return(ub_ATE_51)', .001) ")"
		di "------------"
		di "MTR+MTA Beta 4-3 = (" round(factor_d43 * `return(lb_ATE_43)', .001) "," round(factor_d43 * `return(ub_ATE_43)', .001) ")"
		di "MTR+MTA Beta 4-2 = (" round(factor_d42 * `return(lb_ATE_42)', .001) "," round(factor_d42 * `return(ub_ATE_42)', .001) ")"
		di "MTR+MTA Beta 4-1 = (" round(factor_d41 * `return(lb_ATE_41)', .001) "," round(factor_d41 * `return(ub_ATE_41)', .001) ")"
		di "------------"
		di "MTR+MTA Beta 3-2 = (" round(factor_d32 * `return(lb_ATE_32)', .001) "," round(factor_d32 * `return(ub_ATE_32)', .001) ")"
		di "MTR+MTA Beta 3-1 = (" round(factor_d31 * `return(lb_ATE_31)', .001) "," round(factor_d31 * `return(ub_ATE_21)', .001) ")"
		di "------------"
		di "MTR+MTA Beta 2-1 = (" round(factor_d21 * `return(lb_ATE_21)', .001) "," round(factor_d21 * `return(ub_ATE_21)', .001) ")"

	scalar drop _all
end
