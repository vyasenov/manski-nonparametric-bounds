/*==============================*/
/* 0. PROGRAM FOR ETS ESTIMATES */
/*==============================*/
capture program drop ets_estimates
program ets_estimates, rclass

	version 14
	args y d
	
	sum `d', meanonly
	local max_d = r(max)
	local max_d_1 = r(max) - 1

	di "Potential Outcomes:"
	forvalues i = 1 / `max_d' {		
		sum `y' if `d' == `i', meanonly     // <---------
		return scalar ets_Y`i' = round(r(mean), 0.001)
		di "ETS Y`i' = " return(ets_Y`i')
	}

	forvalues j = 1 / `max_d_1' {	
		return scalar ets_ATE_`max_d'`j' = round(return(ets_Y`max_d') - return(ets_Y`j'), 0.001)
		
		if `j' <= 3 {
		return scalar ets_ATE_4`j' = round(return(ets_Y4) - return(ets_Y`j'), 0.001)
		
		}
		
		if `j' <= 2 {
		return scalar ets_ATE_3`j' = round(return(ets_Y3) - return(ets_Y`j'), 0.001)
		}
		
		return scalar ets_ATE_21 = round(return(ets_Y2) - return(ets_Y1), 0.001)
	}
	
	di " "
	di "------------"
	di "Average Treatment Effects:"
	di "ETS ATE 54 = " return(ets_ATE_54)
	di "ETS ATE 53 = " return(ets_ATE_53)
	di "ETS ATE 52 = " return(ets_ATE_52)
	di "ETS ATE 51 = " return(ets_ATE_51)
	di "------------"
	di "ETS ATE 43 = " return(ets_ATE_43)
	di "ETS ATE 42 = " return(ets_ATE_42)
	di "ETS ATE 41 = " return(ets_ATE_41)
	di "------------"
	di "ETS ATE 32 = " return(ets_ATE_32)
	di "ETS ATE 31 = " return(ets_ATE_31)
	di "------------"
	di "ETS ATE 21 = " return(ets_ATE_21)
	
	end
