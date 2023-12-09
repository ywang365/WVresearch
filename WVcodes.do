* Guaranteed Cash Incentives Boosted COVID-19 Vaccinations of Young Adults: Natural Experiment Evidence from West Virginia

* Content
* 1. Data preparation
* 2. Programs for drawing figures of vax trends & unadjusted vax proportions
* 3. Analysis of the overall effect
*    3.1 Programs for main results with/without kernel function weights
*    3.2 Programs for test for placebo cutoffs 
*    3.3 Programs for test for placebo treated states
* 4. Programs for estimating heterogeneous effects




* Please replace the directories with your own
global dir "..."   //for storing data
global fruit "..."   //for storing output tables/figures


cd $dir
set more off

******************************************************************************
* 1. Data preparation
use $dir\WVdata_raw, clear   //Imported from Household Pulse Survey weeks 22-32 https://www.census.gov/programs-surveys/household-pulse-survey/datasets.html

* Generate two outcome vars
gen vaxEver = 1 if recvdvacc == 1
replace vaxEver = 0 if recvdvacc == 2
lab var vaxEver "Ever vaccinated against COVID-19"

gen vaxFull = 1 if doses == 1
replace vaxFull = 0 if vaxFull == . & vaxEver != . 
lab var vaxFull "Completed or intended to complete the primary series"

* Generate vars for difference-in-discontinuities design
gen age = 2021 - tbirth_year
gen ctrAge = age -35

gen D = ctrAge <= 0   //Identify the treated age group

gen T = week >= 29    //Identify the treated peiod

gen ctrAgeD = ctrAge * D  
gen DT = D * T
gen ctrAgeT = ctrAge * T
gen ctrAgeDT = ctrAgeD * T

* Modify possible covariates
gen race = 1 if rrace == 1
replace race = 2 if rrace == 2
replace race = 3 if inrange(rrace, 3, 4)

gen educ = 1 if inrange(eeduc, 6, 7)
replace educ = 2 if inrange(eeduc, 1, 5)

gen married = 1 if ms == 1
replace married = 2 if inrange(ms, 2, 5)
drop if married == .

gen workCat = 1 if kindwork == 1 | kindwork == 3
replace workCat = 2 if kindwork == 2
replace workCat = 3 if inrange(kindwork, 4, 5)
replace workCat = 4 if inrange(kindwork, -99, -88)  
replace workCat = 5 if anywork == 2  //jobless

gen incCat = 1 if inrange(income, 1, 2)
replace incCat = 2 if inrange(income, 3, 5)
replace incCat = 3 if inrange(income, 6, 8)
replace incCat = 4 if inrange(income, -99, -88)  

gen haskid = thhld_numkid > 0 

replace hadcovid = . if inrange(hadcovid, -99, -88)
drop if hadcovid == .

save $dir\WVdata_refine, replace


******************************************************************************
* 2. Programs for drawing figures of vax trends & unadjusted vax proportions
cap program drop figTrend
program figTrend
	version 15
	syntax, y(string)

	set scheme s1mono
	
	if "`y'" == "vaxEver" {
		local yt1 "Proportion ever vaccinated against COVID-19"
	}
	else if "`y'" == "vaxFull" {
		local yt1 "Proportion that completed"
		local yt2 "or intended to complete primary series"	
	}

	collapse (mean) `y' if inrange(age, 23, 48) [pw = pweight], by(age T)

	# delimit;
	twoway scatter `y' age if T == 0 & age>35, msymbol(O) mcolor("0 168 180") msize(*1.2)  ||
		   scatter `y' age if T == 0 & age<=35, msymbol(O) mcolor("0 168 180") msize(*1.2)  ||
		   scatter `y' age if T == 1 & age>35, msymbol(O) mcolor(dkorange*0.9) msize(*1.2) ||
		   scatter `y' age if T == 1 & age<=35, msymbol(O) mcolor(dkorange*0.9) msize(*1.2) ||
			
		   lfit `y' age if T == 0 & age>35, lcolor("0 168 180") lpattern(solid) lwidth(*1.4) ||
		   lfit `y' age if T == 0 & age<=35, lcolor("0 168 180") lpattern(solid) lwidth(*1.4)  ||
		   lfit `y' age if T == 1 & age>35, lcolor(dkorange*0.9) lpattern(solid) lwidth(*1.4) ||
		   lfit `y' age if T == 1 & age<=35, lcolor(dkorange*0.9) lpattern(solid) lwidth(*1.4) ||
		   
		   , 
		   xline(35, style(extended) lpattern(dash_dot) lcolor(black*0.4))
		   yscale(range(0, 1))
		   ylabel(0(0.1)1, labsize(*0.9) angle(0) format(%02.1f))
		   ytitle("`yt1'" "`yt2'", size(*1.1) linegap(1.5)) xtitle("Age", size(*1.1))
		   xscale(range(23, 48) titlegap(1.5))
		   xlabel(23(1)48, labsize(*0.8) angle(0))
						   
		   legend(order(1 3) label(1 "Pre-intervention") label(3 "Intervention")
				rowgap(*0.5) colgap(*0.7) keygap(*0.5) symxsize(*0.5) symysize(*0.1) size(*0.9) pos(1) ring(0) row(1))
		   ysize(2.7)
		   xsize(4);
		   
	graph export "$fruit\rawTrend-`y'.png", as(png) name("Graph") replace;
   
	#delimit cr
end

* Run the programs
use $dir\WVdata_refine, clear
figTrend, y(vaxEver)

use $dir\WVdata_refine, clear
figTrend, y(vaxFull) 


******************************************************************************
* 3. Analysis of the overall effect

* 3.1 Programs for main results with/without kernel function weights
cap program drop regMain  
program define regMain, eclass
	syntax [, v(string) hL(real 12) hU(real 12) w(string) fig(string))]  //v: varname, hL: lower bound of ctrAge, hU: upper bound of ctrAge, w: kernel function weights, fig(T) for drawing figures that check robustness for different bandwidths
	
	local bL = `hL' + 1
	local bU = `hU' + 1
	
	cap erase "$fruit\R-`v'-kernel`w'-`bL'-`bU'.csv"

	matrix A = J(`=`hU'-`hL'+1', 3, .)

	forvalues h = `hL'/`hU' {
	
		cap drop sup
		gen sup = ctrAge/`h'
		
		if "`w'" == "No" {         //No kernel function weights incorporated.
			cap drop finalW
			gen finalW = pweight
		}	
		else if "`w'" == "gaus" {     //Gaussian kernel function weights
			cap drop gausW
			gen gausW = 1/sqrt(2*_pi)*exp(-0.5*sup^2)
			cap drop finalW
			gen finalW = pweight*gausW
		}
		else if "`w'" == "logi" {      //Logistic kernel function weights
			cap drop logiW
			gen logiW = 1/(exp(sup) + 2 + exp(-sup))
			cap drop finalW
			gen finalW = pweight*logiW
		}
        else if "`w'" == "tria" {     //Triangular kernel function weights
			cap drop triaW
			gen triaW = 1- abs(sup)
			cap drop finalW
			gen finalW = pweight*triaW
		}
		else if "`w'" == "epan" {     //Epanechnikov kernel function weights
			cap drop epanW
			gen epanW = 0.75*(1-sup^2)
			cap drop finalW
			gen finalW = pweight*epanW
		}		
		
		eststo m`h'No:   reg `v' ctrAge D ctrAgeD ctrAgeT DT ctrAgeDT i.week         [pw = finalW] if inrange(ctrAge, -`h', `h'+1), cluster(age)
			cap drop esample 
			gen esample = e(sample) == 1
			
			widget, h(`h') w(`w') co(No) 
			
		eststo m`h'Cov:  reg `v' ctrAge D ctrAgeD ctrAgeT DT ctrAgeDT i.week $covars [pw = finalW] if inrange(ctrAge, -`h', `h'+1), cluster(age)
						
				ereturn local beta = _b[DT]
				ereturn local beta_L = _b[DT] - invttail(e(df_r),0.025)*_se[DT]
				ereturn local beta_U = _b[DT] + invttail(e(df_r),0.025)*_se[DT]
				
				if "`fig'" == "T" {
					mat A[`=`h'-`hL'+1', 1] = _b[DT]
					mat A[`=`h'-`hL'+1', 2] = _b[DT] - invttail(e(df_r),0.025)*_se[DT]
					mat A[`=`h'-`hL'+1', 3] = _b[DT] + invttail(e(df_r),0.025)*_se[DT]
				}
				else {
					dis "Use fig(T) to draw figures."
				}
		
			widget, h(`h') w(`w') co(Cov)			
				
		esttab m* using "$fruit\R-`v'-kernel`w'-`bL'-`bU'.csv", keep(D DT) cells(b(fmt(%10.4f) star) ci(fmt(%10.4f) par) se(fmt(%9.4f) par) p(fmt(%9.4f) par) ) starlevels(* .1 ** .05 *** .01) stats(BandWidth KernelWeight Covars Icpt_PrTr Icpt_PrCt Icpt_PstTr Icpt_PstCt Impact_Icpt_PstTr N, fmt(%4.0f %10s %6s %11.4f %11.4f %11.4f %11.4f %11.4f %9.0f)) collabels(, none) title("Outcome: `v'") legend addnotes(" " " ") append plain    //We added 2 space in notes to make tables separated by 2 rows in Excel.

		eststo clear
	}
	
	if "`fig'" == "T" {
		mat colnames A = lb lb_l95 lb_u95
		clear
		svmat A, names(col)
	
	    outsheet lb lb_l95 lb_u95 using "$fruit\coefci-`v'-`bL'-`bU'.csv", comma replace 
		
		gen bw = _n + `hL'
		
		set scheme s1mono
		
		# delimit;
		twoway line lb bw, lp(solid) lcolor("0 168 180") ||
			   line lb_l95 bw, lp(dash) lcolor("0 168 180") ||
			   line lb_u95 bw, lp(dash) lcolor("0 168 180") ||
			   ,
			   xscale(range(`bL', `bU') titlegap(1))
			   yscale(range(-0.5, 1.0))
			   xlabel(`bL' (1) `bU', labsize(*0.8) angle(0))
			   ylabel(-0.5 (0.1) 1.0, format(%02.1f) labsize(*0.8) angle(0))
			   yline(0, style(extended) lpattern(dash_dot) lcolor(grey))

			   xtitle("Bandwidth", size(*1))
			   ytitle("Coefficient", size(*1))
			   
			   legend(order(1 2) label(1 "Policy effect estimate") label(2 "95% CI") 
					rowgap(*0.4) colgap(*0.2) keygap(*0.4) symxsize(*0.4) symysize(*0.1) size(*1) pos(2) ring(0) row(2))
			   title("", size(*1));
		graph export "$fruit\coef-`v'-`bL'-`bU'.png", as(png) name("Graph") replace;
		# delimit cr
	}
	else {
		dis "Use fig(T) to draw figures."
	}
end

cap program drop widget
program define widget
	syntax [, h(real 12) w(string) co(string)] 
	
	estadd local BandWidth `=`h'+1' : m`h'`co'
	estadd local KernelWeight `w' : m`h'`co'
	estadd local Covars `co' : m`h'`co'
	
	margins if inrange(week, 22, 28) & esample == 1, at(ctrAge = 0 D = 1 ctrAgeD = 0 ctrAgeT = 0 DT = 0 ctrAgeDT = 0)   //T=0 age=35 D=1 
	estadd scalar Icpt_PrTr = r(b)[1, 1] : m`h'`co'    //Icpt_PrTr: intercept that indicates pre-period vax % at age 35 (as it is treated)
			
	margins if inrange(week, 22, 28) & esample == 1, at(ctrAge = 0 D = 0 ctrAgeD = 0 ctrAgeT = 0 DT = 0 ctrAgeDT = 0)   //Counterfactual condition: T=0 age=35 D=0 
	estadd scalar Icpt_PrCt = r(b)[1, 1] : m`h'`co'    //Icpt_PrCt: intercept that indicates pre-period vax % at age 35 (assume it were in control group)
	
	margins if inrange(week, 29, 32) & esample == 1, at(ctrAge = 0 D = 1 ctrAgeD = 0 ctrAgeT = 0 DT = 1 ctrAgeDT = 0)   //T=1 age=35 D=1 
	estadd scalar Icpt_PstTr = r(b)[1, 1] : m`h'`co'    //Icpt_PstTr: intercept that indicates post-period vax % at age 35 (as it is treated)
	scalar Icpt_PstTr = r(b)[1, 1] 
	
	margins if inrange(week, 29, 32) & esample == 1, at(ctrAge = 0 D = 0 ctrAgeD = 0 ctrAgeT = 0 DT = 0 ctrAgeDT = 0)   //Counterfactual condition: T=1 age=35 D=0 
	estadd scalar Icpt_PstCt = r(b)[1, 1] : m`h'`co'    //Icpt_PstCt: intercept that indicates post-period vax % at age 35 (assume it were in control group)
	
	estadd scalar Impact_Icpt_PstTr = _b[DT] / Icpt_PstTr * 100 : m`h'`co'

end


* Run the programs 
global covars i.egender i.rhispanic i.race i.educ i.married i.workCat i.haskid i.hadcovid

* Obtain the overall effect on 2 outcomes (Note: The bandwidth labels in the paper corresponds to h+1 here, so h 1-15 for bandwidth 2-16.)
use $dir\WVdata_refine, clear
regMain, v(vaxEver) hL(1) hU(15) w(No) fig(T)   
use $dir\WVdata_refine, clear
regMain, v(vaxFull) hL(1) hU(15) w(No) fig(T)

* Test for different kernel function weights (bandwidth 13 used by default)
use $dir\WVdata_refine, clear
regMain, v(vaxEver) w(gaus)   
regMain, v(vaxEver) w(logi) 
regMain, v(vaxEver) w(tria)  
regMain, v(vaxEver) w(epan) 

regMain, v(vaxFull) w(gaus)    
regMain, v(vaxFull) w(logi) 
regMain, v(vaxFull) w(tria) 
regMain, v(vaxFull) w(epan) 


*******************************************
* 3.2 Programs for test for placebo cutoffs 
cap program drop placeboCutoffTest
program def placeboCutoffTest
	syntax, v(string) ageL(real) ageU(real) h(real)

	matrix B = J(`=`ageU'-`ageL'+1', 4, .)
	
	forvalues i = `ageL'/`ageU' {

		placeboAge, a(`i')
		reg `v' ctrAge D ctrAgeD ctrAgeT DT ctrAgeDT i.week $covars [pw = pweight] if inrange(ctrAge, -12, 13), cluster(age)
		
		mat B[`=`i'-`ageL'+1', 1] = _b[DT]
		mat B[`=`i'-`ageL'+1', 2] = _b[DT] - invttail(e(df_r),0.025)*_se[DT]
		mat B[`=`i'-`ageL'+1', 3] = _b[DT] + invttail(e(df_r),0.025)*_se[DT]
		mat B[`=`i'-`ageL'+1', 4] = `i'
	}

	mat colnames B = lb lb_l95 lb_u95 age
	clear
	svmat B, names(col)

	outsheet lb lb_l95 lb_u95 age using "$fruit\coefci-`v'-placeboAge-`ageL'-`ageU'.csv", comma replace 
	
	set scheme s1mono
	
	# delimit;
	twoway line lb age, lp(solid) lcolor("0 158 170") ||
		   line lb_l95 age, lp(dash) lcolor("0 158 170") ||
		   line lb_u95 age, lp(dash) lcolor("0 158 170") ||
		   ,
		   xscale(range(`ageL', `ageU') titlegap(2))
		   yscale(range(-0.4, 0.5))
		   xlabel(`ageL' (1) `ageU', labsize(*1.1) angle(0))
		   ylabel(-0.4 (0.1) 0.5, format(%02.1f) labsize(*1) angle(0))
		   yline(0, style(extended) lpattern(dash_dot) lcolor(black) lwidth(vthin))

		   xtitle("Cutoff", size(*1.2))
		   ytitle("", size(*1))
		   
		   legend(order(1 2) label(1 "Coefficient of Interest") label(2 "95% CI")
				rowgap(*0.4) colgap(*0.2) keygap(*0.4) symxsize(*0.4) symysize(*0.1) size(*1) pos(2) ring(0) row(2))
		   title("", size(*1))
		   ysize(2.7)
		   xsize(4);        //BandWidth = `h'
	graph export "$fruit\placeboAge-`v'-`ageL'-`ageU'.png", as(png) name("Graph") replace;
	# delimit cr
end 

cap program drop placeboAge
program def placeboAge
	syntax, a(real)
	
	use $dir\WVdata_refine, clear
	cap drop ctrAge
	cap drop D
	cap drop ctrAgeT
	cap drop ctrAgeD
	cap drop DT
	cap drop ctrAgeDT

	gen ctrAge = age - `a'
	gen D = (age <= `a')
	
	gen ctrAgeD = ctrAge * D 
	gen DT = D * T
	gen ctrAgeT = ctrAge * T
	gen ctrAgeDT = ctrAgeD * T
end

* Run the programs 
global covars i.egender i.rhispanic i.race i.educ i.married i.workCat i.haskid i.hadcovid

use $dir\WVdata_refine, clear
placeboCutoffTest, v(vaxEver) ageL(30) ageU(50) h(12)  //h=12 for bandwidth 13.

use $dir\WVdata_refine, clear
placeboCutoffTest, v(vaxFull) ageL(30) ageU(50) h(12)


*************************************
* 3.3 Programs for test for placebo treated states
cap program drop placeboStateTest
program define placeboStateTest
	version 15
	syntax, v(string)
	
	cap erase "$fruit\R-`v'-placeboTreatedStates.csv"
		
	foreach s of numlist 1/2 4/6 8/13 15/42 44/51 53/56 { 
					
		eststo m`s': reg `v' ctrAge D ctrAgeD ctrAgeT DT ctrAgeDT i.week $covars [pw = pweight] if inrange(ctrAge, -12, 13) & est_st == `s', cluster(age)
			
		esttab m* using "$fruit\R-`v'-placeboTreatedStates.csv", keep(DT) cells("b(fmt(%9.4f) star) se(fmt(%9.4f)) t(fmt(%9.4f)) p(fmt(%9.4f))") starlevels(* .1 ** .05 *** .01) collabels("Coef" "SE" "t" "p-value") mtitles(`s') stats(N, fmt(%9.0f) labels("Obs")) legend append style(fixed)

		eststo clear
	}

end

* Run the programs
global covars i.egender i.rhispanic i.race i.educ i.married i.workCat i.haskid i.hadcovid 

use $dir\WVdata_placeboTreatedStates, clear   //The dataset was generated using similar codes as for generating WVdata_refine.dta.
placeboStateTest, v(vaxEver)
placeboStateTest, v(vaxFull)


*******************************************************************************
* 4. Programs for estimating heterogeneous effects
cap program drop regHetero
program define regHetero
	syntax, v(string) xs(string) xk(string) catv(string) catname(string) categ(string)
	
	cap erase "$fruit\R-`v'-`categ'.csv"
	
	# delimit ;	
	
	eststo mH: reg `v' `xs' i.week $covars [pw = pweight] if inrange(ctrAge, -12, 13), cluster(age);
		widgetHetero, y(`v') catv(`catv') catname(`catname');	
														
	esttab mH using "$fruit\R-`v'-`categ'.csv", keep(`xk') cells(b(fmt(%10.4f) star) ci(fmt(%10.4f) par) se(fmt(%9.4f) par) p(fmt(%9.4f) par) ) starlevels(* .1 ** .05 *** .01) stats(Icpt_PrTr_C0 Icpt_PrCt_C0 Icpt_PstTr_C0 Icpt_PstCt_C0 Icpt_PrTr_C1 Icpt_PrCt_C1 Icpt_PstTr_C1 Icpt_PstCt_C1 N, fmt(%11.4f %11.4f %11.4f %11.4f %11.4f %11.4f %11.4f %11.4f %9.0f)) collabels(, none) title("Outcome: `v'") legend addnotes(" " " ") append plain;

	eststo clear;

	# delimit cr
end

cap program drop widgetHetero 
program define widgetHetero 
	syntax, y(string) catv(string) catname(string)

	* (1)
	margins if inrange(week, 22, 28) & `catname' == 0, at(ctrAge = 0 D = 1 ctrAgeD = 0 ctrAgeT = 0 DT = 0 ctrAgeDT = 0 ctrAge`catv' = 0 D`catv' = 0 ctrAgeD`catv' = 0 ctrAgeT`catv' = 0 DT`catv' = 0 ctrAgeDT`catv' = 0)
	estadd scalar Icpt_PrTr_C0 = r(b)[1, 1] : mH        //C0 for reference category
			
	margins if inrange(week, 22, 28) & `catname' == 0, at(ctrAge = 0 D = 0 ctrAgeD = 0 ctrAgeT = 0 DT = 0 ctrAgeDT = 0 ctrAge`catv' = 0 D`catv' = 0 ctrAgeD`catv' = 0 ctrAgeT`catv' = 0 DT`catv' = 0 ctrAgeDT`catv' = 0)
	estadd scalar Icpt_PrCt_C0 = r(b)[1, 1] : mH
	
	margins if inrange(week, 29, 32) & `catname' == 0, at(ctrAge = 0 D = 1 ctrAgeD = 0 ctrAgeT = 0 DT = 1 ctrAgeDT = 0 ctrAge`catv' = 0 D`catv' = 0 ctrAgeD`catv' = 0 ctrAgeT`catv' = 0 DT`catv' = 0 ctrAgeDT`catv' = 0)
	estadd scalar Icpt_PstTr_C0 = r(b)[1, 1] : mH
	
	margins if inrange(week, 29, 32) & `catname' == 0, at(ctrAge = 0 D = 0 ctrAgeD = 0 ctrAgeT = 0 DT = 0 ctrAgeDT = 0 ctrAge`catv' = 0 D`catv' = 0 ctrAgeD`catv' = 0 ctrAgeT`catv' = 0 DT`catv' = 0 ctrAgeDT`catv' = 0)
	estadd scalar Icpt_PstCt_C0 = r(b)[1, 1] : mH
	
	* (2)
	margins if inrange(week, 22, 28) & `catname' == 1, at(ctrAge = 0 D = 1 ctrAgeD = 0 ctrAgeT = 0 DT = 0 ctrAgeDT = 0 ctrAge`catv' = 0 D`catv' = 1 ctrAgeD`catv' = 0 ctrAgeT`catv' = 0 DT`catv' = 0 ctrAgeDT`catv' = 0)
	estadd scalar Icpt_PrTr_C1 = r(b)[1, 1] : mH        //C1 for the other category
			
	margins if inrange(week, 22, 28) & `catname' == 1, at(ctrAge = 0 D = 0 ctrAgeD = 0 ctrAgeT = 0 DT = 0 ctrAgeDT = 0 ctrAge`catv' = 0 D`catv' = 0 ctrAgeD`catv' = 0 ctrAgeT`catv' = 0 DT`catv' = 0 ctrAgeDT`catv' = 0)
	estadd scalar Icpt_PrCt_C1 = r(b)[1, 1] : mH
	
	margins if inrange(week, 29, 32) & `catname' == 1, at(ctrAge = 0 D = 1 ctrAgeD = 0 ctrAgeT = 0 DT = 1 ctrAgeDT = 0 ctrAge`catv' = 0 D`catv' = 1 ctrAgeD`catv' = 0 ctrAgeT`catv' = 0 DT`catv' = 1 ctrAgeDT`catv' = 0)
	estadd scalar Icpt_PstTr_C1 = r(b)[1, 1] : mH
	
	margins if inrange(week, 29, 32) & `catname' == 1, at(ctrAge = 0 D = 0 ctrAgeD = 0 ctrAgeT = 0 DT = 0 ctrAgeDT = 0 ctrAge`catv' = 0 D`catv' = 0 ctrAgeD`catv' = 0 ctrAgeT`catv' = 0 DT`catv' = 0 ctrAgeDT`catv' = 0)
	estadd scalar Icpt_PstCt_C1 = r(b)[1, 1] : mH
	
end

* Run the programs
* a) Educational attainment
global covars i.egender i.rhispanic i.race i.educ i.married i.workCat i.haskid i.hadcovid 

use $dir\WVdata_refine, clear
gen colgrd = educ == 1  
gen ctrAgeE = ctrAge * colgrd
gen DE = D * colgrd
gen ctrAgeDE = ctrAgeD * colgrd

gen ctrAgeTE = ctrAgeT * colgrd
gen DTE = DT * colgrd
gen ctrAgeDTE = ctrAgeDT * colgrd

regHetero, v(vaxEver) xs(ctrAge ctrAgeE D DE ctrAgeD ctrAgeDE ctrAgeT ctrAgeTE DT DTE ctrAgeDT ctrAgeDTE) xk(DT DTE) catv(E) catname(colgrd) categ(educ)
regHetero, v(vaxFull) xs(ctrAge ctrAgeE D DE ctrAgeD ctrAgeDE ctrAgeT ctrAgeTE DT DTE ctrAgeDT ctrAgeDTE) xk(DT DTE) catv(E) catname(colgrd) categ(educ)

* b) Employment status
global covars i.egender i.rhispanic i.race i.educ i.married i.workCat i.haskid i.hadcovid 

use $dir\WVdata_refine, clear
gen nojob = workCat == 5 
gen ctrAgeNj = ctrAge * nojob
gen DNj = D * nojob
gen ctrAgeDNj = ctrAgeD * nojob

gen ctrAgeTNj = ctrAgeT * nojob
gen DTNj = DT * nojob
gen ctrAgeDTNj = ctrAgeDT * nojob

regHetero, v(vaxEver) xs(ctrAge ctrAgeNj D DNj ctrAgeD ctrAgeDNj ctrAgeT ctrAgeTNj DT DTNj ctrAgeDT ctrAgeDTNj) xk(DT DTNj) catv(Nj) catname(nojob) categ(jobStatus)
regHetero, v(vaxFull) xs(ctrAge ctrAgeNj D DNj ctrAgeD ctrAgeDNj ctrAgeT ctrAgeTNj DT DTNj ctrAgeDT ctrAgeDTNj) xk(DT DTNj) catv(Nj) catname(nojob) categ(jobStatus)

* c) Income level
global covars i.egender i.rhispanic i.race i.educ i.married i.workCat i.haskid i.hadcovid i.incCat 

use $dir\WVdata_refine, clear
drop if incCat == 4    //income unknown
gen lowInc = inrange(income, 1, 2)

gen ctrAgeLo = ctrAge * lowInc
gen DLo = D * lowInc
gen ctrAgeDLo = ctrAgeD * lowInc

gen ctrAgeTLo = ctrAgeT * lowInc
gen DTLo = DT * lowInc
gen ctrAgeDTLo = ctrAgeDT * lowInc

regHetero, v(vaxEver) xs(ctrAge ctrAgeLo D DLo ctrAgeD ctrAgeDLo ctrAgeT ctrAgeTLo DT DTLo ctrAgeDT ctrAgeDTLo) xk(DT DTLo) catv(Lo) catname(lowInc) categ(incomeLevel)
regHetero, v(vaxFull) xs(ctrAge ctrAgeLo D DLo ctrAgeD ctrAgeDLo ctrAgeT ctrAgeTLo DT DTLo ctrAgeDT ctrAgeDTLo) xk(DT DTLo) catv(Lo) catname(lowInc) categ(incomeLevel)

* d) COVID-19 infection history
global covars i.egender i.rhispanic i.race i.educ i.married i.workCat i.haskid i.hadcovid 

use $dir\WVdata_refine, clear
drop if hadcovid == 3   //unsure of whether had covid

gen covidyes = hadcovid == 1 

gen ctrAgeC = ctrAge * covidyes
gen DC = D * covidyes
gen ctrAgeDC = ctrAgeD * covidyes

gen ctrAgeTC = ctrAgeT * covidyes
gen DTC = DT * covidyes
gen ctrAgeDTC = ctrAgeDT * covidyes

regHetero, v(vaxEver) xs(ctrAge ctrAgeC D DC ctrAgeD ctrAgeDC ctrAgeT ctrAgeTC DT DTC ctrAgeDT ctrAgeDTC) xk(DT DTC) catv(C) catname(covidyes) categ(covidHistory)
regHetero, v(vaxFull) xs(ctrAge ctrAgeC D DC ctrAgeD ctrAgeDC ctrAgeT ctrAgeTC DT DTC ctrAgeDT ctrAgeDTC) xk(DT DTC) catv(C) catname(covidyes) categ(covidHistory)





