/*=================================*/
/* 6. PROGRAM FOR IV + MTA BOUNDS */
/*=================================*/
capture program drop iv_mta_bounds
program iv_mta_bounds, rclass 

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
			qui count if `d' >= `i' & `z' == `k'
			local ngeq = r(N)
			
			qui count if `d' <= `i' & `z' == `k'
			local nleq = r(N)
			
			local cond_probgeq = `ngeq' / `n_z`k''		// <---------
			local cond_probleq = `nleq' / `n_z`k''		// <---------
			
			sum `y' if `d' == `i' & `z' == `k', meanonly		// <---------
			
			if `sign' == 2 {   // for MTA+
				gen lb_d`i'z`k' = r(mean) * `cond_probgeq' + `ymin' * (1 - `cond_probgeq')
				gen ub_d`i'z`k' = r(mean) * `cond_probleq' + `ymax' * (1 - `cond_probleq')
			}
			
			else if `sign' == 1 {  // for MTA-
				gen lb_d`i'z`k' = r(mean) * `cond_probleq' + `ymin' * (1 - `cond_probleq')
				gen ub_d`i'z`k' = r(mean) * `cond_probgeq' + `ymax' * (1 - `cond_probgeq')		
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
		
		di "IV+MTA Y`i' = (`return(lb_Y`i')' , `return(ub_Y`i')')"
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
		di "IV+MTA ATE 5-4 = (`return(lb_ATE_54)' , `return(ub_ATE_54)')"	
		di "IV+MTA ATE 5-3 = (`return(lb_ATE_53)' , `return(ub_ATE_53)')"
		di "IV+MTA ATE 5-2 = (`return(lb_ATE_52)' , `return(ub_ATE_52)')"
		di "IV+MTA ATE 5-1 = (`return(lb_ATE_51)' , `return(ub_ATE_51)')"
		di "------------"
		di "IV+MTA ATE 4-3 = (`return(lb_ATE_43)' , `return(ub_ATE_43)')"
		di "IV+MTA ATE 4-2 = (`return(lb_ATE_42)' , `return(ub_ATE_42)')"
		di "IV+MTA ATE 4-1 = (`return(lb_ATE_41)' , `return(ub_ATE_41)')"
		di "------------"
		di "IV+MTA ATE 3-2 = (`return(lb_ATE_32)' , `return(ub_ATE_32)')"
		di "IV+MTA ATE 3-1 = (`return(lb_ATE_31)' , `return(ub_ATE_21)')"
		di "------------"
		di "IV+MTA ATE 2-1 = (`return(lb_ATE_21)' , `return(ub_ATE_21)')"
	
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
		di "IV+MTA Beta 5-4 = (" round(factor_d54 * `return(lb_ATE_54)', .001) "," round(factor_d54 * `return(ub_ATE_54)', .001) ")"	
		di "IV+MTA Beta 5-3 = (" round(factor_d53 * `return(lb_ATE_53)', .001) "," round(factor_d53 * `return(ub_ATE_53)', .001) ")"
		di "IV+MTA Beta 5-2 = (" round(factor_d52 * `return(lb_ATE_52)', .001) "," round(factor_d52 * `return(ub_ATE_52)', .001) ")"
		di "IV+MTA Beta 5-1 = (" round(factor_d51 * `return(lb_ATE_51)', .001) "," round(factor_d51 * `return(ub_ATE_51)', .001) ")"
		di "------------"
		di "IV+MTA Beta 4-3 = (" round(factor_d43 * `return(lb_ATE_43)', .001) "," round(factor_d43 * `return(ub_ATE_43)', .001) ")"
		di "IV+MTA Beta 4-2 = (" round(factor_d42 * `return(lb_ATE_42)', .001) "," round(factor_d42 * `return(ub_ATE_42)', .001) ")"
		di "IV+MTA Beta 4-1 = (" round(factor_d41 * `return(lb_ATE_41)', .001) "," round(factor_d41 * `return(ub_ATE_41)', .001) ")"
		di "------------"
		di "IV+MTA Beta 3-2 = (" round(factor_d32 * `return(lb_ATE_32)', .001) "," round(factor_d32 * `return(ub_ATE_32)', .001) ")"
		di "IV+MTA Beta 3-1 = (" round(factor_d31 * `return(lb_ATE_31)', .001) "," round(factor_d31 * `return(ub_ATE_21)', .001) ")"
		di "------------"
		di "IV+MTA Beta 2-1 = (" round(factor_d21 * `return(lb_ATE_21)', .001) "," round(factor_d21 * `return(ub_ATE_21)', .001) ")"

	scalar drop _all
	end
	
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
	
/*=================================*/
/* 8. PROGRAM FOR MIV + MTR BOUNDS */
/*=================================*/
capture program drop miv_mtr_bounds
program miv_mtr_bounds, rclass 

	version 14
	args y d z sign ymin ymax
	
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
	
	forvalues i = 1 / `max_d' {	
		forvalues k = 1 / `max_z' {
			tempname lb_d`i'z`k' ub_d`i'z`k' n_z`k'
			
			qui count if `z' == `k' 
			gen `n_z`k'' = r(N)
		
			sum `y' if `d' >= `i' & `z' == `k', meanonly 		// <---------
			local meangeq = r(mean) 
			local ngeq = r(N)
			
			sum `y' if `d' <= `i' & `z' == `k', meanonly		// <---------
			local meanleq = r(mean)
			local nleq = r(N)
			
			local cond_probgeq = `ngeq' / `n_z`k''		// <---------
			local cond_probleq = `nleq' / `n_z`k''		// <---------
			
			if `sign' == 1 {  // for MTR-
				gen `lb_d`i'z`k'' = `meangeq' * (`cond_probgeq') + `ymin' * (1 - `cond_probgeq')
				gen `ub_d`i'z`k'' = `meanleq' * (`cond_probleq') + `ymax' * (1 - `cond_probleq')
			}
			
			else if `sign' == 2 {   // for MTR+
				gen `lb_d`i'z`k'' = `meanleq' * (`cond_probleq') + `ymin' * (1 - `cond_probleq')
				gen `ub_d`i'z`k'' = `meangeq' * (`cond_probgeq') + `ymax' * (1 - `cond_probgeq')
			}
			
			else {
				gen `lb_d`i'z`k'' = .
				gen `ub_d`i'z`k'' = .
			}
		}
	}	
	
* this part is not fully automated
	forvalues i = 1 / `max_d' {	
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
		
		di "Potential Outcomes:"
* this part is not fully automated
	forvalues i = 1 / `max_d' {   // potential outcomes
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
			
		di "MIV+MTR Y`i' = (`return(lb_Y`i')' , `return(ub_Y`i')')"
	}

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
		di "MIV+MTR ATE 5-4 = (`return(lb_ATE_54)' , `return(ub_ATE_54)')"	
		di "MIV+MTR ATE 5-3 = (`return(lb_ATE_53)' , `return(ub_ATE_53)')"
		di "MIV+MTR ATE 5-2 = (`return(lb_ATE_52)' , `return(ub_ATE_52)')"
		di "MIV+MTR ATE 5-1 = (`return(lb_ATE_51)' , `return(ub_ATE_51)')"
		di "------------"
		di "MIV+MTR ATE 4-3 = (`return(lb_ATE_43)' , `return(ub_ATE_43)')"
		di "MIV+MTR ATE 4-2 = (`return(lb_ATE_42)' , `return(ub_ATE_42)')"
		di "MIV+MTR ATE 4-1 = (`return(lb_ATE_41)' , `return(ub_ATE_41)')"
		di "------------"
		di "MIV+MTR ATE 3-2 = (`return(lb_ATE_32)' , `return(ub_ATE_32)')"
		di "MIV+MTR ATE 3-1 = (`return(lb_ATE_31)' , `return(ub_ATE_21)')"
		di "------------"
		di "MIV+MTR ATE 2-1 = (`return(lb_ATE_21)' , `return(ub_ATE_21)')"
		
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
		di "MIV+MTR Beta 5-4 = (" round(factor_d54 * `return(lb_ATE_54)', .001) "," round(factor_d54 * `return(ub_ATE_54)', .001) ")"	
		di "MIV+MTR Beta 5-3 = (" round(factor_d53 * `return(lb_ATE_53)', .001) "," round(factor_d53 * `return(ub_ATE_53)', .001) ")"
		di "MIV+MTR Beta 5-2 = (" round(factor_d52 * `return(lb_ATE_52)', .001) "," round(factor_d52 * `return(ub_ATE_52)', .001) ")"
		di "MIV+MTR Beta 5-1 = (" round(factor_d51 * `return(lb_ATE_51)', .001) "," round(factor_d51 * `return(ub_ATE_51)', .001) ")"
		di "------------"
		di "MIV+MTR Beta 4-3 = (" round(factor_d43 * `return(lb_ATE_43)', .001) "," round(factor_d43 * `return(ub_ATE_43)', .001) ")"
		di "MIV+MTR Beta 4-2 = (" round(factor_d42 * `return(lb_ATE_42)', .001) "," round(factor_d42 * `return(ub_ATE_42)', .001) ")"
		di "MIV+MTR Beta 4-1 = (" round(factor_d41 * `return(lb_ATE_41)', .001) "," round(factor_d41 * `return(ub_ATE_41)', .001) ")"
		di "------------"
		di "MIV+MTR Beta 3-2 = (" round(factor_d32 * `return(lb_ATE_32)', .001) "," round(factor_d32 * `return(ub_ATE_32)', .001) ")"
		di "MIV+MTR Beta 3-1 = (" round(factor_d31 * `return(lb_ATE_31)', .001) "," round(factor_d31 * `return(ub_ATE_21)', .001) ")"
		di "------------"
		di "MIV+MTR Beta 2-1 = (" round(factor_d21 * `return(lb_ATE_21)', .001) "," round(factor_d21 * `return(ub_ATE_21)', .001) ")"

	scalar drop _all
	end
		
/*==================================*/
/* 9. PROGRAM FOR MIV + cMTA BOUNDS */
/*==================================*/
capture program drop miv_mta_bounds
program miv_mta_bounds, rclass 

	version 14
	args y d z sign ymin ymax
	
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
	
	forvalues i = 1 / `max_d' {
		forvalues k = 1 / `max_z' {
			tempname lb_d`i'z`k' ub_d`i'z`k'
			
			qui count if `d' >= `i' & `z' == `k'
			local ngeq = r(N)
			
			qui count if `d' <= `i' & `z' == `k'
			local nleq = r(N)
			
			local cond_probgeq = `ngeq' / `n_z`k''		// <---------
			local cond_probleq = `nleq' / `n_z`k''		// <---------
			
			sum `y' if `d' == `i' & `z' == `k', meanonly		// <---------
			
			if `sign' == 2 {  // for MTA+
				gen `lb_d`i'z`k'' = r(mean) * `cond_probgeq' + `ymin' * (1 - `cond_probgeq')
				gen `ub_d`i'z`k'' = r(mean) * `cond_probleq' + `ymax' * (1 - `cond_probleq')
			}
			
			else if `sign' == 1 {  // for MTA-
				gen `lb_d`i'z`k'' = r(mean) * `cond_probleq' + `ymax' * (1 - `cond_probleq')
				gen `ub_d`i'z`k'' = r(mean) * `cond_probgeq' + `ymin' * (1 - `cond_probgeq')
			}
			
			else {
				gen `lb_d`i'z`k'' = .
				gen `ub_d`i'z`k'' = .
			}
		}
	}	
	
* this part is not fully automated
	forvalues i = 1 / `max_d' {	
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
		
		di "Potential Outcomes:"
* this part is not fully automated
	forvalues i = 1 / `max_d' {
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
				
		di "MIV+MTA Y`i' = (`return(lb_Y`i')' , `return(ub_Y`i')')"
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
		di "MIV+MTA ATE 5-4 = (`return(lb_ATE_54)' , `return(ub_ATE_54)')"	
		di "MIV+MTA ATE 5-3 = (`return(lb_ATE_53)' , `return(ub_ATE_53)')"
		di "MIV+MTA ATE 5-2 = (`return(lb_ATE_52)' , `return(ub_ATE_52)')"
		di "MIV+MTA ATE 5-1 = (`return(lb_ATE_51)' , `return(ub_ATE_51)')"
		di "------------"
		di "MIV+MTA ATE 4-3 = (`return(lb_ATE_43)' , `return(ub_ATE_43)')"
		di "MIV+MTA ATE 4-2 = (`return(lb_ATE_42)' , `return(ub_ATE_42)')"
		di "MIV+MTA ATE 4-1 = (`return(lb_ATE_41)' , `return(ub_ATE_41)')"
		di "------------"
		di "MIV+MTA ATE 3-2 = (`return(lb_ATE_32)' , `return(ub_ATE_32)')"
		di "MIV+MTA ATE 3-1 = (`return(lb_ATE_31)' , `return(ub_ATE_21)')"
		di "------------"
		di "MIV+MTA ATE 2-1 = (`return(lb_ATE_21)' , `return(ub_ATE_21)')"
		
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
		di "MIV+MTA Beta 5-4 = (" round(factor_d54 * `return(lb_ATE_54)', .001) "," round(factor_d54 * `return(ub_ATE_54)', .001) ")"	
		di "MIV+MTA Beta 5-3 = (" round(factor_d53 * `return(lb_ATE_53)', .001) "," round(factor_d53 * `return(ub_ATE_53)', .001) ")"
		di "MIV+MTA Beta 5-2 = (" round(factor_d52 * `return(lb_ATE_52)', .001) "," round(factor_d52 * `return(ub_ATE_52)', .001) ")"
		di "MIV+MTA Beta 5-1 = (" round(factor_d51 * `return(lb_ATE_51)', .001) "," round(factor_d51 * `return(ub_ATE_51)', .001) ")"
		di "------------"
		di "MIV+MTA Beta 4-3 = (" round(factor_d43 * `return(lb_ATE_43)', .001) "," round(factor_d43 * `return(ub_ATE_43)', .001) ")"
		di "MIV+MTA Beta 4-2 = (" round(factor_d42 * `return(lb_ATE_42)', .001) "," round(factor_d42 * `return(ub_ATE_42)', .001) ")"
		di "MIV+MTA Beta 4-1 = (" round(factor_d41 * `return(lb_ATE_41)', .001) "," round(factor_d41 * `return(ub_ATE_41)', .001) ")"
		di "------------"
		di "MIV+MTA Beta 3-2 = (" round(factor_d32 * `return(lb_ATE_32)', .001) "," round(factor_d32 * `return(ub_ATE_32)', .001) ")"
		di "MIV+MTA Beta 3-1 = (" round(factor_d31 * `return(lb_ATE_31)', .001) "," round(factor_d31 * `return(ub_ATE_21)', .001) ")"
		di "------------"
		di "MIV+MTA Beta 2-1 = (" round(factor_d21 * `return(lb_ATE_21)', .001) "," round(factor_d21 * `return(ub_ATE_21)', .001) ")"

	scalar drop _all
	end
	
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

/*==============================================*/
/* PROGRAM FOR IMBENS & MANSKI CONFIDENCE SETS  */
/*==============================================*/
clear mata
mata
void myeqb(todo, p, gamma, lnf, g, H)
{
	x   = p[1]
	lnf = (normal(x + gamma)-normal(-x)-0.95 )^2
}

void imconfb(string scalar ngamma)
{
	S = optimize_init()
	optimize_init_evaluator(S, &myeqb())
	optimize_init_evaluatortype(S, "d0")
	optimize_init_params(S, 0)
	gamma = st_numscalar(ngamma)
	optimize_init_argument(S, 1, gamma)
	optimize_init_which(S,  "min" )
	optimize_init_tracelevel(S,"none")	
	optimize_init_conv_ptol(S, 1e-6)   // originally it was 1e-8
	optimize_init_conv_vtol(S, 1e-6)   // originally it was 1e-8
	p = optimize(S)
	st_numscalar("cnew", p)	
}
end
