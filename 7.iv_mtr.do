/*=================================*/
/* 7. PROGRAM FOR IV + MTR BOUNDS */
/*=================================*/
capture program drop iv_mtr_bounds
program iv_mtr_bounds, rclass 

	version 14
	args y d z sign ymin ymax
	
	sum `d', meanonly
	local max_d = r(max)
	local max_d_1 = r(max) - 1
	sum `y', meanonly
	local N = r(N)
	sum `z', meanonly
	local max_z = r(max)
	
	forvalues k = 1 / `max_z' {
		tempname n_z`k'	
		qui count if `z' == `k' 
		gen `n_z`k'' = r(N)
	}
	
	forvalues i = 1 / `max_d' {	
		forvalues k = 1 / `max_z' {
			sum `y' if `d' >= `i' & `z' == `k', meanonly 		// <---------
			local meangeq = r(mean) 
			local ngeq = r(N)
			
			sum `y' if `d' <= `i' & `z' == `k', meanonly		// <---------
			local meanleq = r(mean)
			local nleq = r(N)
			
			local cond_probgeq = `ngeq' / `n_z`k''		// <---------
			local cond_probleq = `nleq' / `n_z`k''		// <---------
			
			if `sign' == 1 {  // for MTR-
				gen lb_d`i'z`k' = `meangeq' * (`cond_probgeq') + `ymin' * (1 - `cond_probgeq')
				gen ub_d`i'z`k' = `meanleq' * (`cond_probleq') + `ymax' * (1 - `cond_probleq')
			}
			
			else if `sign' == 2 {  // for MTR+
				gen lb_d`i'z`k' = `meanleq' * (`cond_probleq') + `ymin' * (1 - `cond_probleq')
				gen ub_d`i'z`k' = `meangeq' * (`cond_probgeq') + `ymax' * (1 - `cond_probgeq')
			}
			
			else {
				gen lb_d`i'z`k' = .
				gen ub_d`i'z`k' = .
			}
		}
	}	
	
	di "Potential Outcomes:"
	forvalues i = 1 / `max_d' {	
		tempvar aux1 aux2

		egen `aux1' = rowmax(lb_d`i'z*)
		egen `aux2' = rowmin(ub_d`i'z*)
		
		return scalar lb_Y`i' = round(`aux1', 0.001)
		return scalar ub_Y`i' = round(`aux2', 0.001)		
		
		di "IV+MTR Y`i' = (`return(lb_Y`i')' , `return(ub_Y`i')')"
	}
	
	drop lb* ub*
	
forvalues j = 1 / `max_d_1' {  // treatment effects
		
		if `sign' == 1 {   // make sure the upper bound of MTR- sign is non-positive
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
		
		
		else if `sign' == 2 {   // make sure the lower bound of MTR+ sign is non-negative
		
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
		di "IV+MTR ATE 5-4 = (`return(lb_ATE_54)' , `return(ub_ATE_54)')"	
		di "IV+MTR ATE 5-3 = (`return(lb_ATE_53)' , `return(ub_ATE_53)')"
		di "IV+MTR ATE 5-2 = (`return(lb_ATE_52)' , `return(ub_ATE_52)')"
		di "IV+MTR ATE 5-1 = (`return(lb_ATE_51)' , `return(ub_ATE_51)')"
		di "------------"
		di "IV+MTR ATE 4-3 = (`return(lb_ATE_43)' , `return(ub_ATE_43)')"
		di "IV+MTR ATE 4-2 = (`return(lb_ATE_42)' , `return(ub_ATE_42)')"
		di "IV+MTR ATE 4-1 = (`return(lb_ATE_41)' , `return(ub_ATE_41)')"
		di "------------"
		di "IV+MTR ATE 3-2 = (`return(lb_ATE_32)' , `return(ub_ATE_32)')"
		di "IV+MTR ATE 3-1 = (`return(lb_ATE_31)' , `return(ub_ATE_21)')"
		di "------------"
		di "IV+MTR ATE 2-1 = (`return(lb_ATE_21)' , `return(ub_ATE_21)')"
		
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
		di "IV+MTR Beta 5-4 = (" round(factor_d54 * `return(lb_ATE_54)', .001) "," round(factor_d54 * `return(ub_ATE_54)', .001) ")"	
		di "IV+MTR Beta 5-3 = (" round(factor_d53 * `return(lb_ATE_53)', .001) "," round(factor_d53 * `return(ub_ATE_53)', .001) ")"
		di "IV+MTR Beta 5-2 = (" round(factor_d52 * `return(lb_ATE_52)', .001) "," round(factor_d52 * `return(ub_ATE_52)', .001) ")"
		di "IV+MTR Beta 5-1 = (" round(factor_d51 * `return(lb_ATE_51)', .001) "," round(factor_d51 * `return(ub_ATE_51)', .001) ")"
		di "------------"
		di "IV+MTR Beta 4-3 = (" round(factor_d43 * `return(lb_ATE_43)', .001) "," round(factor_d43 * `return(ub_ATE_43)', .001) ")"
		di "IV+MTR Beta 4-2 = (" round(factor_d42 * `return(lb_ATE_42)', .001) "," round(factor_d42 * `return(ub_ATE_42)', .001) ")"
		di "IV+MTR Beta 4-1 = (" round(factor_d41 * `return(lb_ATE_41)', .001) "," round(factor_d41 * `return(ub_ATE_41)', .001) ")"
		di "------------"
		di "IV+MTR Beta 3-2 = (" round(factor_d32 * `return(lb_ATE_32)', .001) "," round(factor_d32 * `return(ub_ATE_32)', .001) ")"
		di "IV+MTR Beta 3-1 = (" round(factor_d31 * `return(lb_ATE_31)', .001) "," round(factor_d31 * `return(ub_ATE_21)', .001) ")"
		di "------------"
		di "IV+MTR Beta 2-1 = (" round(factor_d21 * `return(lb_ATE_21)', .001) "," round(factor_d21 * `return(ub_ATE_21)', .001) ")"

	scalar drop _all
end
