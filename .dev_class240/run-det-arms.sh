#!/usr/bin/env bash
## Every configuration again, deterministic (all rows in both models, nearestWeighted, seed 1).
set -u
ROOT=/mnt/fast/class240-regression
RUN=$HOME/GitHub/worktrees/class240-scanfi/Biomass_borealDataPrep/.dev_class240/run-one-scanfi.R
DEV=$HOME/GitHub/worktrees/class240-dev/Biomass_borealDataPrep
NEW=$HOME/GitHub/worktrees/class240-scanfi/Biomass_borealDataPrep
COV=$HOME/GitHub/worktrees/class240-covernum/Biomass_borealDataPrep
BL=$HOME/GitHub/worktrees/landr-baseline-dev
NL=$HOME/GitHub/worktrees/landr-240-scanfi
LOG=$ROOT/scratch/runs-s2020
arm () { # label area landr module stratum wetland
  echo "$(date '+%F %T') START $2 $1" | tee -a $LOG/summary-det.txt
  ( cd "$4" && CLASS240_DETERMINISTIC=1 CLASS240_STRATUM="$5" CLASS240_WETLAND="$6" Rscript "$RUN" "$1" "$2" "$3" "$4" ) > "$LOG/$2__$1.log" 2>&1
  echo "$(date '+%F %T') END   $2 $1 exit=$?" | tee -a $LOG/summary-det.txt
}
for a in big_s2020 small_s2020; do
  arm det-baseline     $a $BL $DEV ""              0
  arm det-forestland   $a $NL $DEV ""              0
  arm det-landcover    $a $NL $NEW landcover       1
  arm det-sitecomp     $a $NL $NEW siteComposition 1
  arm det-landcover-cov $a $NL $COV landcover      1
  arm det-sitecomp-cov $a $NL $COV siteComposition 1
done
echo "$(date '+%F %T') DET DONE" | tee -a $LOG/summary-det.txt
