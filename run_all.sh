#!/usr/bin/env bash
# =============================================================================
# run_all.sh
# -----------------------------------------------------------------------------
# Runs the full empirical pipeline for the electoral-turnover project:
#
#   1. build_panels.py  (Python)
#   2. run_rd.do        (Stata)
#   3. run_didm.do      (Stata)
#
# Per-step output is captured under A_MicroData/logs/ with timestamps so any
# failure can be diagnosed afterwards by reading the latest log. The wrapper
# log A_MicroData/logs/run_all.log records the high-level timing and exit
# codes for every step.
#
# After this script exits 0, ask Claude to "compile the landcover report" to
# regenerate D_Reports/landcover_results.tex via the report skill, then run
# pdflatex twice on it.
#
# Usage (from the repo root):
#   bash run_all.sh             # full pipeline
#   bash run_all.sh --only rd   # only the RD step
#   bash run_all.sh --skip rd   # skip the RD step
# =============================================================================

set -uo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$REPO_ROOT"

LOG_DIR="A_MicroData/logs"
mkdir -p "$LOG_DIR" "A_MicroData/results" "A_MicroData/figs/rd" "A_MicroData/figs/didm"

TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
WRAPPER_LOG="$LOG_DIR/run_all_${TIMESTAMP}.log"
ln -sfn "$(basename "$WRAPPER_LOG")" "$LOG_DIR/run_all.log"

# ----------------------------------------------------------------------------
# argument parsing
# ----------------------------------------------------------------------------
ONLY=""
SKIP=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --only) ONLY="$2"; shift 2 ;;
    --skip) SKIP="$2"; shift 2 ;;
    -h|--help)
      sed -n '2,28p' "$0"
      exit 0 ;;
    *) echo "unknown flag: $1" >&2; exit 2 ;;
  esac
done

should_run() {
  local step="$1"
  if [[ -n "$ONLY" && "$ONLY" != "$step" ]]; then return 1; fi
  if [[ -n "$SKIP" && "$SKIP" == "$step" ]]; then return 1; fi
  return 0
}

# ----------------------------------------------------------------------------
# locate Stata binary
# ----------------------------------------------------------------------------
STATA_BIN=""
for cand in stata stata-mp stata-se stata-be StataMP StataSE Stata; do
  if command -v "$cand" >/dev/null 2>&1; then STATA_BIN="$cand"; break; fi
done

# ----------------------------------------------------------------------------
# logging helpers
# ----------------------------------------------------------------------------
COLOR_RESET=$'\033[0m'
COLOR_BLUE=$'\033[34m'
COLOR_GREEN=$'\033[32m'
COLOR_RED=$'\033[31m'
COLOR_YELLOW=$'\033[33m'

log_w() { printf '%s [run_all] %s\n' "$(date +%H:%M:%S)" "$*" | tee -a "$WRAPPER_LOG"; }

banner() {
  local title="$1"
  printf '\n%s================ %s ================%s\n' \
    "$COLOR_BLUE" "$title" "$COLOR_RESET" | tee -a "$WRAPPER_LOG"
}

step() {
  # step <name> <cmd...>
  local name="$1"; shift
  if ! should_run "$name"; then
    log_w "$COLOR_YELLOW[skip]$COLOR_RESET $name"
    return 0
  fi
  banner "STEP: $name"
  local log="$LOG_DIR/${name}_${TIMESTAMP}.log"
  ln -sfn "$(basename "$log")" "$LOG_DIR/${name}.log"
  log_w "command : $*"
  log_w "log     : $log"
  local t0; t0=$(date +%s)
  ( "$@" ) >"$log" 2>&1
  local rc=$?
  local dt=$(( $(date +%s) - t0 ))
  if [[ $rc -ne 0 ]]; then
    log_w "${COLOR_RED}[FAIL]${COLOR_RESET} $name (rc=$rc, ${dt}s) -- see tail of $log:"
    tail -n 30 "$log" | sed 's/^/  | /' | tee -a "$WRAPPER_LOG"
    exit "$rc"
  fi
  log_w "${COLOR_GREEN}[ok]${COLOR_RESET} $name (${dt}s)"
}

# ----------------------------------------------------------------------------
# preflight
# ----------------------------------------------------------------------------
banner "preflight"
log_w "repo root      : $REPO_ROOT"
log_w "wrapper log    : $WRAPPER_LOG"
log_w "python         : $(command -v python || echo 'missing')"
log_w "stata binary   : ${STATA_BIN:-not found in PATH}"
log_w "raw elections  : $(test -f B_RawData/elections_long_clean_with_details.dta && echo present || echo MISSING)"
log_w "IGBP CSV       : $(test -f A_MicroData/IGBP_LandCover_Municipios_Colombia_2001_2024.csv && echo present || echo MISSING)"
log_w "VCF CSV        : $(test -f A_MicroData/VCF_TreeCover_Municipios_Colombia_2000_2025.csv && echo present || echo MISSING)"

if ! command -v python >/dev/null 2>&1; then
  log_w "${COLOR_RED}python missing${COLOR_RESET}"; exit 3
fi
if [[ -z "$STATA_BIN" ]]; then
  log_w "${COLOR_YELLOW}stata not on PATH${COLOR_RESET} -- Stata steps will be skipped."
fi

# ----------------------------------------------------------------------------
# pipeline
# ----------------------------------------------------------------------------
step build_panels  python C_Programs/build_panels.py

if [[ -n "$STATA_BIN" ]]; then
  # Stata batch mode (-b) writes its own .log next to the .do; we additionally
  # capture stdout/stderr to the per-step log for completeness.
  step rd   "$STATA_BIN" -b do C_Programs/run_rd.do
  step didm "$STATA_BIN" -b do C_Programs/run_didm.do
  # Move Stata's auto-generated batch logs into A_MicroData/logs/.
  for f in run_rd.log run_didm.log; do
    if [[ -f "$f" ]]; then
      mv -f "$f" "$LOG_DIR/${f%.log}_stata_${TIMESTAMP}.log"
      ln -sfn "${f%.log}_stata_${TIMESTAMP}.log" "$LOG_DIR/${f%.log}_stata.log"
    fi
  done
else
  log_w "${COLOR_YELLOW}skipping rd${COLOR_RESET} -- no Stata binary"
  log_w "${COLOR_YELLOW}skipping didm${COLOR_RESET} -- no Stata binary"
fi

# ----------------------------------------------------------------------------
# summary
# ----------------------------------------------------------------------------
banner "summary"
log_w "artefacts:"
for f in \
  A_MicroData/analysis_panel.parquet \
  A_MicroData/analysis_panel.dta \
  A_MicroData/balance_panel.parquet \
  A_MicroData/balance_panel.dta \
  A_MicroData/build_panels_info.csv \
  A_MicroData/results/rd_results.csv \
  A_MicroData/results/didm_results.csv \
  A_MicroData/results/didm_summary.csv ; do
  if [[ -f "$f" ]]; then
    size="$(wc -c < "$f" | tr -d ' ')"
    log_w "$(printf '  ok %10s B  %s' "$size" "$f")"
  else
    log_w "  -- missing       $f"
  fi
done

n_rd_png=$(find A_MicroData/figs/rd   -name '*.png' 2>/dev/null | wc -l | tr -d ' ')
n_dm_png=$(find A_MicroData/figs/didm -name '*.png' 2>/dev/null | wc -l | tr -d ' ')
log_w "RD plots   : $n_rd_png PNGs in A_MicroData/figs/rd/"
log_w "DID-M plots: $n_dm_png PNGs in A_MicroData/figs/didm/"

banner "next step"
log_w "ask Claude: \"compile the landcover report\""
log_w "then: pdflatex D_Reports/landcover_results.tex && pdflatex D_Reports/landcover_results.tex"
log_w "wrapper log: $WRAPPER_LOG"
