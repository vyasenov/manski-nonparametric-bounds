/*==================================*/
/* 1. PROGRAM FOR WORST CASE BOUNDS */
/*==================================*/
capture program drop worst_case_bounds
program worst_case_bounds, rclass
 
	version 14
	args y d ymin ymax
	
*	syntax varlist   // this is if we want to use the syntax command instead of args
*	local y `1'
*	local d `2'
*	local ymin `3'
*	local ymax `4'
	
	sum `d', meanonly
	local max_d = r(max)
	local max_d_1 = r(max) - 1
	sum `y', meanonly
	local N = r(N)
	
	di "Potential Outcomes:"
	forvalues i = 1 / `max_d' {		
		sum `y' if `d' == `i', meanonly			// <---------
		local mean = r(mean)
		local n = r(N)
		local prob = `n' / `N'			// <---------
		
		return scalar lb_Y`i' = round(`mean' * `prob' + `ymin' * (1 -`prob'), 0.001)
		return scalar ub_Y`i' = round(`mean' * `prob' + `ymax' * (1 -`prob'), 0.001)
		
	*	di "Worst Case Y`i' = (`return(lb_Y`i')' , `return(ub_Y`i')')"
	}
	
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
		di "Worst Case ATE 5-4 = (`return(lb_ATE_54)' , `return(ub_ATE_54)')"	
		di "Worst Case ATE 5-3 = (`return(lb_ATE_53)' , `return(ub_ATE_53)')"
		di "Worst Case ATE 5-2 = (`return(lb_ATE_52)' , `return(ub_ATE_52)')"
		di "Worst Case ATE 5-1 = (`return(lb_ATE_51)' , `return(ub_ATE_51)')"
		di "------------"
		di "Worst Case ATE 4-3 = (`return(lb_ATE_43)' , `return(ub_ATE_43)')"
		di "Worst Case ATE 4-2 = (`return(lb_ATE_42)' , `return(ub_ATE_42)')"
		di "Worst Case ATE 4-1 = (`return(lb_ATE_41)' , `return(ub_ATE_41)')"
		di "------------"
		di "Worst Case ATE 3-2 = (`return(lb_ATE_32)' , `return(ub_ATE_32)')"
		di "Worst Case ATE 3-1 = (`return(lb_ATE_31)' , `return(ub_ATE_21)')"
		di "------------"
		di "Worst Case ATE 2-1 = (`return(lb_ATE_21)' , `return(ub_ATE_21)')"
	
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
		return scalar ave_beta_lb = round( factor_d51 * `return(lb_ATE_51)'+ factor_d52 * `return(lb_ATE_52)'+ factor_d53 * `return(lb_ATE_53)'+ factor_d54 * `return(lb_ATE_54)', 0.001)/4
		return scalar ave_beta_ub = round( factor_d51 * `return(ub_ATE_51)'+ factor_d52 * `return(ub_ATE_52)'+ factor_d53 * `return(ub_ATE_53)'+ factor_d54 * `return(ub_ATE_54)', 0.001)/4

		di " "
		di "Betas:"
		di "Worst Case Beta 5-4 = (" round(factor_d54 * `return(lb_ATE_54)', .001) "," round(factor_d54 * `return(ub_ATE_54)', .001) ")"	
		di "Worst Case Beta 5-3 = (" round(factor_d53 * `return(lb_ATE_53)', .001) "," round(factor_d53 * `return(ub_ATE_53)', .001) ")"
		di "Worst Case Beta 5-2 = (" round(factor_d52 * `return(lb_ATE_52)', .001) "," round(factor_d52 * `return(ub_ATE_52)', .001) ")"
		di "Worst Case Beta 5-1 = (" round(factor_d51 * `return(lb_ATE_51)', .001) "," round(factor_d51 * `return(ub_ATE_51)', .001) ")"
		di "Worst Case Average Beta = (" `return(ave_beta_lb)' " , "`return(ave_beta_ub)' " )"
		di "------------"
		di "Worst Case Beta 4-3 = (" round(factor_d43 * `return(lb_ATE_43)', .001) "," round(factor_d43 * `return(ub_ATE_43)', .001) ")"
		di "Worst Case Beta 4-2 = (" round(factor_d42 * `return(lb_ATE_42)', .001) "," round(factor_d42 * `return(ub_ATE_42)', .001) ")"
		di "Worst Case Beta 4-1 = (" round(factor_d41 * `return(lb_ATE_41)', .001) "," round(factor_d41 * `return(ub_ATE_41)', .001) ")"
		di "------------"
		di "Worst Case Beta 3-2 = (" round(factor_d32 * `return(lb_ATE_32)', .001) "," round(factor_d32 * `return(ub_ATE_32)', .001) ")"
		di "Worst Case Beta 3-1 = (" round(factor_d31 * `return(lb_ATE_31)', .001) "," round(factor_d31 * `return(ub_ATE_21)', .001) ")"
		di "------------"
		di "Worst Case Beta 2-1 = (" round(factor_d21 * `return(lb_ATE_21)', .001) "," round(factor_d21 * `return(ub_ATE_21)', .001) ")"

	scalar drop _all
end
