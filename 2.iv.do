/*==========================*/
/* 2. PROGRAM FOR IV BOUNDS */
/*==========================*/
capture program drop iv_bounds
program iv_bounds, rclass 

	version 14
	args y d z ymin ymax
	
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
			sum `y' if `d' == `i' & `z' == `k', meanonly		// <---------
			local cond_prob = r(N) / `n_z`k''		// <---------
			
			gen lb_d`i'z`k' = r(mean) * `cond_prob' + `ymin' * (1 - `cond_prob')
			gen ub_d`i'z`k' = r(mean) * `cond_prob' + `ymax' * (1 - `cond_prob')
		}
	}	
	
	di "Potential Outcomes:"
	forvalues i = 1 / `max_d' {	
		tempvar aux1 aux2

		egen `aux1' = rowmax(lb_d`i'z*)
		egen `aux2' = rowmin(ub_d`i'z*)
		
		return scalar lb_Y`i' = round(`aux1', 0.001)
		return scalar ub_Y`i' = round(`aux2', 0.001)		
		
		di "IV Y`i' = (`return(lb_Y`i')' , `return(ub_Y`i')')"
	}
	
	drop lb* ub*
	
	forvalues j = 1 / `max_d_1' {	
	
		return scalar lb_ATE_`max_d'`j' = round(return(lb_Y`max_d') - return(ub_Y`j'), 0.001)
		return scalar ub_ATE_`max_d'`j' = round(return(ub_Y`max_d') - return(lb_Y`j'), 0.001)	
		
		if `j' <= 3 {
			return scalar lb_ATE_4`j' = round(return(lb_Y4) - return(ub_Y`j'), 0.001)
			return scalar ub_ATE_4`j' = round(return(ub_Y4) - return(lb_Y`j'), 0.001)
		}	

		if `j' <= 2 {
			return scalar lb_ATE_3`j' = round(return(lb_Y3) - return(ub_Y`j'), 0.001)
			return scalar ub_ATE_3`j' = round(return(ub_Y3) - return(lb_Y`j'), 0.001)
		}
	}	
		return scalar lb_ATE_21 = round(return(lb_Y2) - return(ub_Y1), 0.001)
		return scalar ub_ATE_21 = round(return(ub_Y2) - return(lb_Y1), 0.001)	

		di " "
		di "Average Treatment Effects:"	
		di "IV ATE 5-4 = (`return(lb_ATE_54)' , `return(ub_ATE_54)')"	
		di "IV ATE 5-3 = (`return(lb_ATE_53)' , `return(ub_ATE_53)')"
		di "IV ATE 5-2 = (`return(lb_ATE_52)' , `return(ub_ATE_52)')"
		di "IV ATE 5-1 = (`return(lb_ATE_51)' , `return(ub_ATE_51)')"
		di "------------"
		di "IV ATE 4-3 = (`return(lb_ATE_43)' , `return(ub_ATE_43)')"
		di "IV ATE 4-2 = (`return(lb_ATE_42)' , `return(ub_ATE_42)')"
		di "IV ATE 4-1 = (`return(lb_ATE_41)' , `return(ub_ATE_41)')"
		di "------------"
		di "IV ATE 3-2 = (`return(lb_ATE_32)' , `return(ub_ATE_32)')"
		di "IV ATE 3-1 = (`return(lb_ATE_31)' , `return(ub_ATE_21)')"
		di "------------"
		di "IV ATE 2-1 = (`return(lb_ATE_21)' , `return(ub_ATE_21)')"
		
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
		di "IV Beta 5-4 = (" round(factor_d54 * `return(lb_ATE_54)', .001) "," round(factor_d54 * `return(ub_ATE_54)', .001) ")"	
		di "IV Beta 5-3 = (" round(factor_d53 * `return(lb_ATE_53)', .001) "," round(factor_d53 * `return(ub_ATE_53)', .001) ")"
		di "IV Beta 5-2 = (" round(factor_d52 * `return(lb_ATE_52)', .001) "," round(factor_d52 * `return(ub_ATE_52)', .001) ")"
		di "IV Beta 5-1 = (" round(factor_d51 * `return(lb_ATE_51)', .001) "," round(factor_d51 * `return(ub_ATE_51)', .001) ")"
		di "------------"
		di "IV Beta 4-3 = (" round(factor_d43 * `return(lb_ATE_43)', .001) "," round(factor_d43 * `return(ub_ATE_43)', .001) ")"
		di "IV Beta 4-2 = (" round(factor_d42 * `return(lb_ATE_42)', .001) "," round(factor_d42 * `return(ub_ATE_42)', .001) ")"
		di "IV Beta 4-1 = (" round(factor_d41 * `return(lb_ATE_41)', .001) "," round(factor_d41 * `return(ub_ATE_41)', .001) ")"
		di "------------"
		di "IV Beta 3-2 = (" round(factor_d32 * `return(lb_ATE_32)', .001) "," round(factor_d32 * `return(ub_ATE_32)', .001) ")"
		di "IV Beta 3-1 = (" round(factor_d31 * `return(lb_ATE_31)', .001) "," round(factor_d31 * `return(ub_ATE_21)', .001) ")"
		di "------------"
		di "IV Beta 2-1 = (" round(factor_d21 * `return(lb_ATE_21)', .001) "," round(factor_d21 * `return(ub_ATE_21)', .001) ")"

	scalar drop _all
end
