#!/usr/bin/env bash
## Full rerun on the integration builds (all 9 PRs), against current development.
set -u
ROOT=/mnt/fast/class240-regression
H=$HOME/GitHub/worktrees/bbdp-integ/Biomass_borealDataPrep/.dev_class240
INT_MOD=$HOME/GitHub/worktrees/bbdp-integ/Biomass_borealDataPrep
INT_L=$HOME/GitHub/worktrees/landr-integ
DEV_MOD=$HOME/GitHub/worktrees/class240-dev/Biomass_borealDataPrep
DEV_L=$HOME/GitHub/worktrees/landr-dev
MAN=$HOME/arbutus_manifest_SCANFI_v2_clean.csv
LOG=$ROOT/scratch/runs-integ; mkdir -p $LOG; S=$LOG/summary.txt
step () { echo "$(date '+%F %T') $*" | tee -a $S; }
step "START fixtures (integration module + integration LandR, LandR's own mirror)"
( cd $INT_MOD && CLASS240_FIXTURE_SUFFIX=_int CLASS240_REMAP_MANIFEST="" CLASS240_MODULE_PATH=$(dirname $INT_MOD) CLASS240_LANDR_BASE=$INT_L \
    Rscript $H/prepare-inputs-scanfi.R all ) > $LOG/fixtures.log 2>&1
step "END   fixtures exit=$?"
arm () { # label area landr module stratum wetland det manifest
  step "START $2 $1"
  ( cd "$4" && CLASS240_STRATUM="$5" CLASS240_WETLAND="$6" CLASS240_DETERMINISTIC="$7" CLASS240_REMAP_MANIFEST="$8" \
      Rscript $H/run-one-scanfi.R "$1" "$2" "$3" "$4" ) > "$LOG/$2__$1.log" 2>&1
  step "END   $2 $1 exit=$?"
}
for a in big_s2020_int small_s2020_int; do
  arm base-det      $a $DEV_L $DEV_MOD ""              0 1 "$MAN"
  arm int-lc-det    $a $INT_L $INT_MOD landcover       1 1 ""
  arm int-sc-det    $a $INT_L $INT_MOD siteComposition 1 1 ""
done
## default mode (no harness overrides): subsetSeed alone must make these identical
arm int-lc-default-1 big_s2020_int $INT_L $INT_MOD landcover 1 0 ""
arm int-lc-default-2 big_s2020_int $INT_L $INT_MOD landcover 1 0 ""
step "ALL DONE"
