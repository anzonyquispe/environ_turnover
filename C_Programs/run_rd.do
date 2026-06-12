/*=========================================================================
  run_rd.do
  -------------------------------------------------------------------------
  Close-elections RD analysis for the electoral-turnover / land-cover
  project. Uses the rdrobust Stata package (Calonico, Cattaneo, Titiunik).

  Three RD designs are estimated on the same analysis panel:

    RD-1a Incumbent-ideology change (score_inc).
          Sample : 2015, 2019, 2023 (inc_ideology_win observed).
          D = 1 if winner's ideology differs from previous incumbent.

    RD-1b Incumbent re-election (score_inc_won).
          Sample : 2015, 2019, 2023 (inc_won observed; 2011 masked in
          build_panels.py because the raw column is mechanically 0 for
          every 2011 row).
          score_inc_won = margin x inc_won_signed where
          inc_won_signed = +1 if the incumbent did NOT win again
          (turnover) and -1 if the incumbent re-won.
          D = 1{score_inc_won > 0} -> turnover at the cutoff.

    RD-2  Three ideology dichotomies of the winner.
          Sample : 2011, 2015, 2019, 2023 (winner ideology observed,
          dropping the "No info" category).
          Three transformations:
            ideo_l_vs_rc :  Left  vs  Right or Center
            ideo_lc_vs_r :  Right vs  Left or Center
            ideo_lr_vs_c :  Center vs Left or Right
          For each, score = margin * ideo_signed, D = 1{score > 0}.

  Each design also runs a McCrary (rddensity) test per (outcome, sample)
  to check for manipulation of the running variable at the cutoff.

  Specification: rdrobust with p=1, triangular kernel, MSE-optimal
  bandwidth (bwselect("mserd")), SEs clustered at the municipality.
  Reported: bias-corrected point estimate paired with the robust
  standard error and p-value ("Robust" line of rdrobust output).
  Heterogeneity samples: Full / Above / Below 2010 median of the
  outcome's own baseline value.

  Outputs:
    A_MicroData/results/rd_results.csv      long table of all RD cells
                                            (design x outcome x horizon
                                            x sample)
    A_MicroData/results/rd_density.csv      McCrary (rddensity) test
                                            per (design x outcome x sample)
    A_MicroData/figs/rd/rdplot_<design>_<outcome>_<sample>.png
                                            rdplot truncated to the
                                            MSE-optimal bandwidth, where
                                            <design> is RD1a / RD1b /
                                            RD2_<transform>

  Run from the project root:  do C_Programs/run_rd.do
=========================================================================*/

clear all
set more off
version 17
cd "/Users/anzony.quisperojas/Documents/GitHub/environ_turnover"
cap mkdir "A_MicroData/results"
cap mkdir "A_MicroData/figs"
cap mkdir "A_MicroData/figs/rd"
cap log close
log using "A_MicroData/results/run_rd.log", replace text

* one-time install
cap which rdrobust
if _rc ssc install rdrobust, replace
cap which rdplot
if _rc ssc install rdrobust, replace
cap which rddensity
if _rc ssc install rddensity, replace

use "A_MicroData/analysis_panel.dta", clear
egen mpio_id = group(mpio)

local OUTCOMES forest mixed_forest shrublands savannas grassland agriculture crop_nature tree_cover mean_night
local HORIZONS main k1 k2 k3 pre3
local SAMPLES all above below

* postfile for results
postfile rd_post str20 design str20 outcome str10 horizon str10 sample ///
    double coef double se double pvalue double h_l double h_r ///
    long n_h long n_total int p_degree str10 kernel ///
    using "A_MicroData/results/rd_results.dta", replace

* postfile for McCrary density tests (one per design x outcome x sample)
postfile rd_den str20 design str20 outcome str10 sample ///
    double T_q double pv_q long N_l long N_r ///
    using "A_MicroData/results/rd_density.dta", replace

* ----------------------------------------------------------------------
*   helper to apply the sample mask consistently
* ----------------------------------------------------------------------
program drop _all
program define _apply_sample
    args outcome sample
    if "`sample'" == "above"  keep if above_median_`outcome'_2010 == 1
    if "`sample'" == "below"  keep if above_median_`outcome'_2010 == 0
end

* ======================================================================
*   RD-1a: incumbent-IDEOLOGY change (score_inc)
*          sample restricted to elections with inc_ideology_win observed
*          (2015 / 2019 / 2023)
* ======================================================================
foreach v of local OUTCOMES {
    foreach h of local HORIZONS {
        foreach s of local SAMPLES {
            preserve
            keep if !missing(score_inc) & !missing(std_dY_`v'_`h')
            _apply_sample `v' `s'
            cap rdrobust std_dY_`v'_`h' score_inc, ///
                c(0) p(1) kernel(triangular) bwselect(mserd) ///
                vce(cluster mpio_id)
            if _rc {
                post rd_post ("RD-1a") ("`v'") ("`h'") ("`s'") ///
                    (.) (.) (.) (.) (.) (0) (0) (1) ("triangular")
            }
            else {
                local coef = e(tau_bc)
                local se   = e(se_tau_rb)
                local pval = e(pv_rb)
                local hl   = e(h_l)
                local hr   = e(h_r)
                local nh   = e(N_h_l) + e(N_h_r)
                local nt   = e(N)
                post rd_post ("RD-1a") ("`v'") ("`h'") ("`s'") ///
                    (`coef') (`se') (`pval') (`hl') (`hr') ///
                    (`nh') (`nt') (1) ("triangular")
            }
            restore
        }
    }

    * rdplot for the main horizon, restricted to selected MSE bandwidth
    foreach s of local SAMPLES {
        preserve
        keep if !missing(score_inc) & !missing(std_dY_`v'_main)
        _apply_sample `v' `s'
        cap rdrobust std_dY_`v'_main score_inc, c(0) p(1) ///
            kernel(triangular) bwselect(mserd) vce(cluster mpio_id)
        if !_rc {
            local h = e(h_l)
            keep if abs(score_inc) <= `h'
            cap rdplot std_dY_`v'_main score_inc, c(0) p(1) ///
                kernel(triangular) binselect(esmv) ///
                graph_options(title("RD-1a `v', sample=`s', |score|<=`=string(`h', "%6.2f")'") ///
                              xtitle("score_inc = margin x inc\_ideology\_win\_signed") ///
                              ytitle("std_dY_`v'_main"))
            if !_rc graph export "A_MicroData/figs/rd/rdplot_RD1a_`v'_`s'.png", replace width(1200)
        }
        restore
    }

    * McCrary density test (one per outcome x sample at the running variable)
    foreach s of local SAMPLES {
        preserve
        keep if !missing(score_inc)
        _apply_sample `v' `s'
        cap rddensity score_inc, c(0)
        if _rc {
            post rd_den ("RD-1a") ("`v'") ("`s'") (.) (.) (0) (0)
        }
        else {
            local T  = e(T_q)
            local pv = e(pv_q)
            local Nl = e(N_l)
            local Nr = e(N_r)
            post rd_den ("RD-1a") ("`v'") ("`s'") (`T') (`pv') (`Nl') (`Nr')
        }
        restore
    }
}

* ======================================================================
*   RD-1b: incumbent RE-ELECTION (score_inc_won)
*          score = margin * inc_won_signed; positive score = turnover
*          (incumbent did NOT win), negative = incumbent re-won.
*          Sample 2015 / 2019 / 2023 (2011 masked in build_panels.py).
* ======================================================================
foreach v of local OUTCOMES {
    foreach h of local HORIZONS {
        foreach s of local SAMPLES {
            preserve
            keep if !missing(score_inc_won) & !missing(std_dY_`v'_`h')
            _apply_sample `v' `s'
            cap rdrobust std_dY_`v'_`h' score_inc_won, ///
                c(0) p(1) kernel(triangular) bwselect(mserd) ///
                vce(cluster mpio_id)
            if _rc {
                post rd_post ("RD-1b") ("`v'") ("`h'") ("`s'") ///
                    (.) (.) (.) (.) (.) (0) (0) (1) ("triangular")
            }
            else {
                local coef = e(tau_bc)
                local se   = e(se_tau_rb)
                local pval = e(pv_rb)
                local hl   = e(h_l)
                local hr   = e(h_r)
                local nh   = e(N_h_l) + e(N_h_r)
                local nt   = e(N)
                post rd_post ("RD-1b") ("`v'") ("`h'") ("`s'") ///
                    (`coef') (`se') (`pval') (`hl') (`hr') ///
                    (`nh') (`nt') (1) ("triangular")
            }
            restore
        }
    }

    foreach s of local SAMPLES {
        preserve
        keep if !missing(score_inc_won) & !missing(std_dY_`v'_main)
        _apply_sample `v' `s'
        cap rdrobust std_dY_`v'_main score_inc_won, c(0) p(1) ///
            kernel(triangular) bwselect(mserd) vce(cluster mpio_id)
        if !_rc {
            local h = e(h_l)
            keep if abs(score_inc_won) <= `h'
            cap rdplot std_dY_`v'_main score_inc_won, c(0) p(1) ///
                kernel(triangular) binselect(esmv) ///
                graph_options(title("RD-1b `v', sample=`s', |score|<=`=string(`h', "%6.2f")'") ///
                              xtitle("score_inc_won = margin x inc\_won\_signed (positive = turnover)") ///
                              ytitle("std_dY_`v'_main"))
            if !_rc graph export "A_MicroData/figs/rd/rdplot_RD1b_`v'_`s'.png", replace width(1200)
        }
        restore
    }

    foreach s of local SAMPLES {
        preserve
        keep if !missing(score_inc_won)
        _apply_sample `v' `s'
        cap rddensity score_inc_won, c(0)
        if _rc {
            post rd_den ("RD-1b") ("`v'") ("`s'") (.) (.) (0) (0)
        }
        else {
            local T  = e(T_q)
            local pv = e(pv_q)
            local Nl = e(N_l)
            local Nr = e(N_r)
            post rd_den ("RD-1b") ("`v'") ("`s'") (`T') (`pv') (`Nl') (`Nr')
        }
        restore
    }
}

* ======================================================================
*   RD-2: three ideology transformations
*         sample includes 2011 (winner's own ideology observed)
* ======================================================================
local IDEO_TRANSFORMS l_vs_rc lc_vs_r lr_vs_c

foreach tr of local IDEO_TRANSFORMS {
    foreach v of local OUTCOMES {
        foreach h of local HORIZONS {
            foreach s of local SAMPLES {
                preserve
                keep if !missing(score_`tr') & !missing(std_dY_`v'_`h')
                _apply_sample `v' `s'
                cap rdrobust std_dY_`v'_`h' score_`tr', ///
                    c(0) p(1) kernel(triangular) bwselect(mserd) ///
                    vce(cluster mpio_id)
                if _rc {
                    post rd_post ("RD-2:`tr'") ("`v'") ("`h'") ("`s'") ///
                        (.) (.) (.) (.) (.) (0) (0) (1) ("triangular")
                }
                else {
                    local coef = e(tau_bc)
                    local se   = e(se_tau_rb)
                    local pval = e(pv_rb)
                    local hl   = e(h_l)
                    local hr   = e(h_r)
                    local nh   = e(N_h_l) + e(N_h_r)
                    local nt   = e(N)
                    post rd_post ("RD-2:`tr'") ("`v'") ("`h'") ("`s'") ///
                        (`coef') (`se') (`pval') (`hl') (`hr') ///
                        (`nh') (`nt') (1) ("triangular")
                }
                restore
            }
        }

        * rdplot for the main horizon
        foreach s of local SAMPLES {
            preserve
            keep if !missing(score_`tr') & !missing(std_dY_`v'_main)
            _apply_sample `v' `s'
            cap rdrobust std_dY_`v'_main score_`tr', c(0) p(1) ///
                kernel(triangular) bwselect(mserd) vce(cluster mpio_id)
            if !_rc {
                local h = e(h_l)
                keep if abs(score_`tr') <= `h'
                cap rdplot std_dY_`v'_main score_`tr', c(0) p(1) ///
                    kernel(triangular) binselect(esmv) ///
                    graph_options(title("RD-2:`tr' `v', sample=`s', |score|<=`=string(`h', "%6.2f")'") ///
                                  xtitle("score_`tr' = margin x ideo_`tr'") ///
                                  ytitle("std_dY_`v'_main"))
                if !_rc graph export "A_MicroData/figs/rd/rdplot_RD2_`tr'_`v'_`s'.png", replace width(1200)
            }
            restore
        }

        * McCrary density test
        foreach s of local SAMPLES {
            preserve
            keep if !missing(score_`tr')
            _apply_sample `v' `s'
            cap rddensity score_`tr', c(0)
            if _rc {
                post rd_den ("RD-2:`tr'") ("`v'") ("`s'") (.) (.) (0) (0)
            }
            else {
                local T  = e(T_q)
                local pv = e(pv_q)
                local Nl = e(N_l)
                local Nr = e(N_r)
                post rd_den ("RD-2:`tr'") ("`v'") ("`s'") (`T') (`pv') (`Nl') (`Nr')
            }
            restore
        }
    }
}

postclose rd_post
postclose rd_den

* Export both result tables as CSV for the report skill.
use "A_MicroData/results/rd_results.dta", clear
export delimited using "A_MicroData/results/rd_results.csv", replace
use "A_MicroData/results/rd_density.dta", clear
export delimited using "A_MicroData/results/rd_density.csv", replace

di _newline as result "RD analysis complete. Results in A_MicroData/results/, "  ///
                       "figures in A_MicroData/figs/rd/."
log close
