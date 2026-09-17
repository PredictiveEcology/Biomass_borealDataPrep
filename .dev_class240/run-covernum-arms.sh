#!/usr/bin/env bash
## Extra arms: the class-240 rework PLUS the coverNum pixel-count fix (feat/class240-scanfi-covernum),
## to isolate that fix's effect against the landcover / sitecomp arms. Sequential, like the main driver.
set -u
ROOT=/mnt/fast/class240-regression
RUN=$HOME/GitHub/worktrees/class240-scanfi/Biomass_borealDataPrep/.dev_class240/run-one-scanfi.R
MOD=$HOME/GitHub/worktrees/class240-covernum/Biomass_borealDataPrep
LANDR=$HOME/GitHub/worktrees/landr-240-scanfi
LOGDIR=$ROOT/scratch/runs-s2020
for area in big_s2020 small_s2020; do
  for spec in "landcover-cov landcover" "sitecomp-cov siteComposition"; do
    set -- $spec; label=$1; stratum=$2
    echo "$(date '+%F %T') START $area $label" | tee -a "$LOGDIR/summary.txt"
    ( cd "$MOD" && CLASS240_STRATUM="$stratum" CLASS240_WETLAND=1 Rscript "$RUN" "$label" "$area" "$LANDR" "$MOD" ) > "$LOGDIR/${area}__${label}.log" 2>&1
    echo "$(date '+%F %T') END   $area $label exit=$?" | tee -a "$LOGDIR/summary.txt"
  done
done
echo "$(date '+%F %T') COVERNUM DONE" | tee -a "$LOGDIR/summary.txt"
