*********************************************************************
* Project: 	Coverage and associated factors for non-compliance 		*
*			to flu vaccination in elders from peri-urban settings	*
*			of northern Peru.										*
* Authors: 	Roberto Niño-Garcia, Ludwing A. Zeta, 					*
*			Luis M. Helguero-Santin, Carlos Culquichicon			*
*			Andres G. Lescano										*
* Data analysists: 	Carlos Culquichicon, Roberto Niño-Garcia		* 
*					and Andres G. Lescano							*
* Last update:		02/05/19										*
*********************************************************************

version 15.1
local dir `autocd'
cd "`autocd'"
use "db_vaccflu_v4.dta", clear
set more off
drop num_enf 
*log using "log_vaccflu",replace

** Cleaning-up data
recode estadocivil 		(0 3 4 = 0 "Sin pareja") (1 2 = 1 "Con pareja"), ///
						gen(civil_rec)
recode nivel_educativo 	(0 1 2 = 0 "Sin carrera profesional") ///
						(3 4 = 1 "Con carrera profesional"), gen (edu_rec)
recode con_quien_vive   (0 = 0 "Sin cuidador") ///
						(1 2 3 = 1 "Con cuidador"), gen (cuidador)
recode enfer_agrup		(0 = 0 "Ninguna") (1 2 =1 "1 o mas"), gen(cronicas)
recode motiv_no_vacuna  (0 = 0 "Under-awareness utility of vaccination") ///
                        (1 3= 1 "Misinformation about vaccine availability") ///
						(2 = 2 "Afraid of vaccination") ///
						(4 = 3 "Vaccine shortage") ///
						(5 6 7 8 9 10 = 4 "Iteraction with health conditions"), ///
						gen (motiv_agrup)
global ssgen 		sexo civil_rec edu_rec
global sssocial 	cuidador trabaja fuma alcohol ejercicios
global sssalud 		alergia_pollo buena_salud cronicas tipo_servicio hosp_2015 

** Descriptive analysis
tab1 $ssgen $sssocial $sssalud
sum edad, d
graph box edad, over (vacunacion_2015)
swilk vacunacion_2015

** Bivariate analysis
ttest edad, by(vacunacion_2015)

foreach var in $ssgen {
tab `var' vacunacion_2015, chi2 row exp
}

foreach var in $sssocial {
tab `var' vacunacion_2015, chi2 row
}

foreach var in $sssalud {
tab `var' vacunacion_2015, chi2 row
}

** Regressions
* Null model
poisson vacunacion_2015 $ssgen $sssocial $sssalud, vce(robust)
gen nomiss = e(sample)
eststo m_0:  poisson 	vacunacion_2015 	if nomiss==1,nolog irr

* Level 1
eststo m_0:  poisson 	vacunacion_2015 	if nomiss==1,nolog irr
foreach i in $ssgen $sssocial $sssalud {
eststo m_`i': quietly poisson 	vacunacion_2015 `i' if nomiss==1,nolog irr
lrtest m_0 m_`i'
}

esttab m_*, varwidth(25) star (* 0.05 ** 0.01 *** 0.001) ///
	compress nogaps p  stats(p df chi2 N) ///
	mtitles("MA0" "MA1" "MA2" "MA3" "MA4" "MA5" "MA6" "MA7" "MA8" "MA9" ///
			"MA10" "MA11" "MA12" "MA13") ///
	addnotes("Exponentiated coeficients are Prevalence Ratio")
eststo clear
/* trabaja selected */

* Level 2
eststo m_0: quietly poisson 	vacunacion_2015 	trabaja ///
								if nomiss==1,nolog irr
foreach i in $ssgen cuidador fuma alcohol ejercicios $sssalud {
eststo m_`i': quietly poisson 	vacunacion_2015 	trabaja `i' ///
								if nomiss==1,nolog irr
quietly lrtest m_0 m_`i'
}

esttab m_*, varwidth(25) star (+ 0.10 * 0.05 ** 0.01 *** 0.001) ///
	compress nogaps p  stats(p df chi2 N) ///
	mtitles("MA0" "MA1" "MA2" "MA3" "MA4" "MA5" "MA6" "MA7" "MA8" "MA9" ///
			"MA10" "MA11" "MA12" "MA13") ///
	addnotes("Exponentiated coeficients are Prevalence Ratio") 
eststo clear
/* civil_rec selected */

* Level 3
eststo m_0: quietly poisson 	vacunacion_2015 	trabaja civil_rec ///
								if nomiss==1,nolog irr
foreach i in sexo edu_rec cuidador fuma alcohol ejercicios $sssalud {
eststo m_`i': quietly poisson 	vacunacion_2015 	trabaja civil_rec `i' ///
								if nomiss==1,nolog irr
quietly lrtest m_0 m_`i'
}

esttab m_*, varwidth(25) star (+ 0.10 * 0.05 ** 0.01 *** 0.001) ///
	compress nogaps p stats(p df chi2 N) ///
	mtitles("MA0" "MA1" "MA2" "MA3" "MA4" "MA5" "MA6" "MA7" "MA8" "MA9" ///
			"MA10" "MA11" "MA12" "MA13") ///
	addnotes("Exponentiated coeficients are Prevalence Ratio") 
eststo clear
/* no-variable selecteded */


*> Crude analysis
foreach var in $ssgen {
glm vacunacion_2015 i.`var',fam(poisson) link(log) eform nolog
}
glm vacunacion_2015 ib1.edu_rec,fam(poisson) link(log) eform nolog

foreach var in $sssocial {
glm vacunacion_2015 i.`var',fam(poisson) link(log) eform nolog
}
glm vacunacion_2015 ib1.trabaja,fam(poisson) link(log) eform nolog 

foreach var in $sssalud {
glm vacunacion_2015 i.`var',fam(poisson) link(log) eform nolog
}

*> Multivariable model
poisson 	vacunacion_2015 ib1.trabaja i.civil_rec, nolog irr
foreach i in 	sexo edu_rec cuidador fuma alcohol ejercicios alergia_pollo ///
				buena_salud cronicas tipo_servicio hosp_2015 {
poisson 	vacunacion_2015 ib1.trabaja i.civil_rec i.`i', nolog irr
}

poisson 	vacunacion_2015 ib1.trabaja i.civil_rec ib1.edu_rec, nolog irr
poisson 	vacunacion_2015 ib1.trabaja i.civil_rec, nolog irr


** Exploring differences between non-adherence groups
$sssalud $sssocial $sssalud
foreach i in vacunacion_2015{
	
	foreach v in 0 1 2 3 4 {
	sktest `i' if motiv_agrup==`v' & motiv_agrup!=.
	}
	
tabstat `i', by(motiv_agrup) stat (n mean sd p50 sk k)
kwallis `i', by(motiv_agrup)
di ("-------------------------------------------------------------------------")
}

exit
