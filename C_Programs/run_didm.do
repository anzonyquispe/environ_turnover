/*=========================================================================
  run_didm.do
  -------------------------------------------------------------------------
  De Chaisemartin & D'Haultfoeuille (DCDH) DID_M switchers-in event
  study via the `did_multiplegt_dyn` Stata package.

  Panel    : balance_panel.dta (mpio x year, 2001-2024).
  Treatment: D = party-level mayoral turnover broadcast from the most
             recent election (zero before 2015; switches at 2015/2019/
             2023).
  Outcomes : 8 in total (7 IGBP land-cover layers + VCF tree_cover).
  Spec     : effects(10) placebo(4) cluster(mpio).
  Samples  : Full / Above / Below own-outcome 2010 baseline median.

  For each (outcome, sample) the do-file:
    1. estimates the dynamic effects and placebos with graph_off
    2. retrieves the effect / placebo matrices and the average total
       effect (and its SE) plus the placebo joint-test p-value
    3. draws an event-study plot annotated with the average total
       effect (top left) and the placebo joint-test p-value
    4. exports the raw coefficient table as CSV.

  Outputs:
    A_MicroData/results/didm_results.csv  long table (lag, coef, se, ...)
    A_MicroData/figs/didm/didm_<outcome>_<sample>.png

  Run from project root:  do C_Programs/run_didm.do
=========================================================================*/

clear all
set more off
version 17
cd "/Users/anzony.quisperojas/Documents/GitHub/environ_turnover"
cap mkdir "A_MicroData/results"
cap mkdir "A_MicroData/figs"
cap mkdir "A_MicroData/figs/didm"
cap log close
log using "A_MicroData/results/run_didm.log", replace text

cap which did_multiplegt_dyn
if _rc ssc install did_multiplegt_dyn, replace

use "A_MicroData/balance_panel.dta", clear
egen mpio_id = group(mpio)
xtset mpio_id year

local OUTCOMES forest mixed_forest shrublands savannas grassland agriculture crop_nature tree_cover mean_night
local SAMPLES all above below
local EFFECTS  10
local PLACEBOS 4

postfile didm_post str20 outcome str10 sample int lag str10 kind ///
    double coef double se long n_swi using ///
    "A_MicroData/results/didm_results.dta", replace

postfile didm_sum str20 outcome str10 sample double avg_te double avg_te_se ///
    double p_joint_placebo long n_total using ///
    "A_MicroData/results/didm_summary.dta", replace

foreach v of local OUTCOMES {
    foreach s of local SAMPLES {
        preserve
        if "`s'" == "above" keep if above_median_`v'_2010 == 1
        if "`s'" == "below" keep if above_median_`v'_2010 == 0

        di _newline as result ///
            "================ outcome=`v'  sample=`s' ================"

        cap did_multiplegt_dyn `v' mpio_id year D, ///
            effects(`EFFECTS') placebo(`PLACEBOS') ///
            cluster(mpio_id) switchers(in) graph_off
        if _rc {
            di as error "did_multiplegt_dyn failed for `v', `s' (rc=`_rc')"
            restore
            continue
        }

        local n_total = e(N)

        * Average total effect and its SE; placebo joint-test p-value.
        * Scalar names verified from the package's ereturn list.
        local avg_te    = e(Av_tot_effect)
        local avg_te_se = e(se_avg_total_effect)
        local p_pl      = e(p_jointplacebo)
        local p_effs    = e(p_jointeffects)

        * Post the dynamic effects and placebos, pulling each from its
        * dedicated scalar (e.g. e(Effect_3), e(se_effect_3), ...).
        forvalues l = 1/`EFFECTS' {
            local coef  = e(Effect_`l')
            local se    = e(se_effect_`l')
            local n_swi = e(N_switchers_effect_`l')
            if missing(`coef') local coef = .
            if missing(`se')   local se   = .
            if missing(`n_swi') local n_swi = 0
            post didm_post ("`v'") ("`s'") (`l') ("effect") ///
                (`coef') (`se') (`n_swi')
        }
        forvalues k = 1/`PLACEBOS' {
            local coef  = e(Placebo_`k')
            local se    = e(se_placebo_`k')
            local n_swi = e(N_switchers_placebo_`k')
            if missing(`coef') local coef = .
            if missing(`se')   local se   = .
            if missing(`n_swi') local n_swi = 0
            post didm_post ("`v'") ("`s'") (-`k') ("placebo") ///
                (`coef') (`se') (`n_swi')
        }
        post didm_sum ("`v'") ("`s'") (`avg_te') (`avg_te_se') ///
                      (`p_pl') (`n_total')

        * ---- event-study plot ----
        * Pull all coefficients into locals BEFORE entering the frame,
        * so the plot block doesn't depend on e() return values being
        * preserved across frame switches.
        forvalues l = 1/`EFFECTS' {
            local C_eff_`l' = cond(missing(e(Effect_`l')), ., e(Effect_`l'))
            local S_eff_`l' = cond(missing(e(se_effect_`l')), ., e(se_effect_`l'))
        }
        forvalues k = 1/`PLACEBOS' {
            local C_pl_`k' = cond(missing(e(Placebo_`k')), ., e(Placebo_`k'))
            local S_pl_`k' = cond(missing(e(se_placebo_`k')), ., e(se_placebo_`k'))
        }
        local avg_te_s   = string(`avg_te', "%6.3f")
        local avg_te_ses = string(`avg_te_se', "%6.3f")
        local p_pl_s     = string(`p_pl', "%6.3f")

        cap frame drop _plot
        frame create _plot
        frame _plot {
            local NROWS = `EFFECTS' + `PLACEBOS' + 1
            set obs `NROWS'
            gen lag  = .
            gen coef = .
            gen se   = .

            * Placebos: Placebo_k -> lag = -k. Walk from -PLACEBOS to -1
            * so the lag axis is sorted ascending.
            local i 1
            forvalues k = `PLACEBOS'(-1)1 {
                replace lag  = -`k' in `i'
                replace coef = `C_pl_`k'' in `i'
                replace se   = `S_pl_`k'' in `i'
                local ++i
            }
            * Omitted reference at lag = 0 (coef = 0, no CI).
            replace lag  = 0 in `i'
            replace coef = 0 in `i'
            replace se   = 0 in `i'
            local ++i
            * Effects: Effect_l -> lag = +l, l = 1..EFFECTS.
            forvalues l = 1/`EFFECTS' {
                replace lag  = `l' in `i'
                replace coef = `C_eff_`l'' in `i'
                replace se   = `S_eff_`l'' in `i'
                local ++i
            }

            gen ci_lo = coef - 1.96 * se
            gen ci_hi = coef + 1.96 * se
            * No CI band for the omitted reference (visualises as a single dot).
            replace ci_lo = . if lag == 0
            replace ci_hi = . if lag == 0

            sort lag

            twoway (rcap ci_lo ci_hi lag, lcolor(gs8)) ///
                   (connected coef lag, msymbol(O) ///
                        mcolor(navy) lcolor(navy) lwidth(medthin)), ///
                yline(0, lcolor(gs10) lwidth(thin)) ///
                xline(0.5, lpattern(dot) lcolor(gs10)) ///
                xlabel(`=-`PLACEBOS''(1)`EFFECTS', labsize(small)) ///
                xtitle("Event time {it:l} (years from switch)") ///
                ytitle("DID{sub:l}{sup:+} on `v' (pp)") ///
                title("DID{sub:M} `v', sample=`s'") ///
                note("Avg total effect = `avg_te_s' (SE `avg_te_ses')" ///
                     "Joint placebo test p = `p_pl_s'", ///
                     size(small) position(11) ring(0)) ///
                legend(off) ///
                graphregion(color(white))
            graph export "A_MicroData/figs/didm/didm_`v'_`s'.png", replace width(1300)
        }
        frame drop _plot

        restore
    }
}

postclose didm_post
postclose didm_sum

use "A_MicroData/results/didm_results.dta", clear
export delimited using "A_MicroData/results/didm_results.csv", replace
use "A_MicroData/results/didm_summary.dta", clear
export delimited using "A_MicroData/results/didm_summary.csv", replace

di _newline as result "DID_M analysis complete. Tables in A_MicroData/results/, " ///
                       "figures in A_MicroData/figs/didm/."
log close
