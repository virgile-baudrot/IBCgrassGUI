#!/usr/bin/env bash
# Start the control designs of §4.3 on Scaleway Serverless Jobs: one job run per (biome, design),
# 5 control runs on 6 vCPU each. Outputs go to s3://$BUCKET/$PREFIX/4.3_control_design/.
#   bash Input-EFSA/scaleway/submit_4.3.sh <prefix> [design ...]
#   e.g. bash Input-EFSA/scaleway/submit_4.3.sh control-2026-10-28 full s50 s25 s10
#        BIOMES="BOR_I1" bash Input-EFSA/scaleway/submit_4.3.sh test full
# Designs are those of run_4.3_control_design.R ("nodom" needs the "full" outputs in the same
# job and cannot run alone here). A failed job run is started again with the same arguments.
set -euo pipefail
JOB_ID="${JOB_ID:-b7b866d3-b307-4011-816c-9acfdeacd3af}" # job definition "ibcgrass-efsa"
BUCKET="${BUCKET:-ibcgrass-efsa-runs}"
prefix="${1:?prefix (e.g. control-2026-10-28)}"
if [ $# -gt 1 ]; then designs=("${@:2}"); else designs=(full s50 s25 s10); fi
read -r -a biomes <<< "${BIOMES:-STE_I1 BLS_I1 ALP_I1 ATL_I1 BOR_I1 CON_I1 PAN_I1 MED_I1 MAC_I1}"
log="submitted_${prefix}.tsv"
for b in "${biomes[@]}"; do
  for d in "${designs[@]}"; do
    id=$(scw jobs definition start "$JOB_ID" \
      environment-variables.SCRIPT=run_4.3_control_design.R \
      environment-variables.ARGS="$b $d" environment-variables.NCORES=6 \
      environment-variables.BUCKET="$BUCKET" environment-variables.PREFIX="$prefix/4.3_control_design" \
      -o json | python3 -I -c 'import json,sys;d=json.load(sys.stdin);print(d["job_runs"][0]["id"])')
    printf '%s\t%s\t%s\n' "$b" "$d" "$id" | tee -a "$log"
  done
done
echo "job run ids in $log; follow with: scw jobs run list -o human"
