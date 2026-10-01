#!/usr/bin/env bash
## Run every SCANFI arm, one after another (in parallel they collide on downloads into the
## shared inputs/ directory). Each arm's log and exit status are kept; a failed arm does not
## stop the others.
##
##   baseline    dev module (+ @optinfo guard)  x baseline LandR     -- the status quo
##   forestland  dev module (+ @optinfo guard)  x new LandR          -- #221's rule alone
##   landcover   new module, stratumType=landcover       x new LandR + CWIM
##   sitecomp    new module, stratumType=siteComposition x new LandR + CWIM
##
## Usage: run-all-scanfi.sh [area ...]   (default: big_s2020 small_s2020)

set -u
ROOT=${CLASS240_ROOT:-/mnt/fast/class240-regression}
HERE=$(cd "$(dirname "$0")" && pwd)
RUN="$HERE/run-one-scanfi.R"
DEV_MOD=$HOME/GitHub/worktrees/class240-dev/Biomass_borealDataPrep
NEW_MOD=$HOME/GitHub/worktrees/class240-scanfi/Biomass_borealDataPrep
BASE_LANDR=$HOME/GitHub/worktrees/landr-baseline-dev
NEW_LANDR=$HOME/GitHub/worktrees/landr-240-scanfi
LOGDIR=$ROOT/scratch/runs-s2020
mkdir -p "$LOGDIR"
AREAS=("$@")
[ ${#AREAS[@]} -eq 0 ] && AREAS=(big_s2020 small_s2020)

run_arm () {  # label area landr module stratum wetland
  local label=$1 area=$2 landr=$3 module=$4 stratum=$5 wetland=$6
  local log="$LOGDIR/${area}__${label}.log"
  echo "$(date '+%F %T') START $area $label" | tee -a "$LOGDIR/summary.txt"
  ( cd "$module" && CLASS240_STRATUM="$stratum" CLASS240_WETLAND="$wetland" \
      Rscript "$RUN" "$label" "$area" "$landr" "$module" ) > "$log" 2>&1
  local rc=$?
  echo "$(date '+%F %T') END   $area $label exit=$rc" | tee -a "$LOGDIR/summary.txt"
}

for area in "${AREAS[@]}"; do
  run_arm baseline   "$area" "$BASE_LANDR" "$DEV_MOD" ""                0
  run_arm forestland "$area" "$NEW_LANDR"  "$DEV_MOD" ""                0
  run_arm landcover  "$area" "$NEW_LANDR"  "$NEW_MOD" landcover         1
  run_arm sitecomp   "$area" "$NEW_LANDR"  "$NEW_MOD" siteComposition   1
done
echo "$(date '+%F %T') ALL DONE" | tee -a "$LOGDIR/summary.txt"
