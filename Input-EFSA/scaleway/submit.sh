#!/usr/bin/env bash
# Start the §4.4 design on Scaleway Serverless Jobs: one job run per community (biome, MC),
# 32 model runs on 6 vCPU each (40 to 100 min). Outputs go to s3://$BUCKET/$PREFIX/.
#   bash Input-EFSA/scaleway/submit.sh <prefix> [nMC] [biome ...]
#   e.g. bash Input-EFSA/scaleway/submit.sh main-2026-10-30 10
#        bash Input-EFSA/scaleway/submit.sh dryrun 2 BOR_I1
#        LEVELS=z06,z07,z08,z09 bash Input-EFSA/scaleway/submit.sh z06-z09-2026-10-27 10
#        (treatments of these levels only, 24 runs per job run, see run_4.4_stressor_effects.R)
# A job run that fails is simply started again with the same (biome, MC): same seed, same result.
set -euo pipefail
JOB_ID="${JOB_ID:?set JOB_ID to the id of your Scaleway job definition}"
BUCKET="${BUCKET:?set BUCKET to your Object Storage bucket}"
prefix="${1:?prefix (e.g. main-2026-10-30)}"
nmc="${2:-10}"
if [ $# -gt 2 ]; then biomes=("${@:3}"); else biomes=(STE_I1 BLS_I1 ALP_I1 ATL_I1 BOR_I1 CON_I1 PAN_I1 MED_I1 MAC_I1); fi
log="submitted_${prefix}.tsv"
for b in "${biomes[@]}"; do
  for mc in $(seq 1 "$nmc"); do
    id=$(scw jobs definition start "$JOB_ID" \
      environment-variables.SCRIPT=run_4.4_stressor_effects.R \
      environment-variables.ARGS="$b $mc" environment-variables.NCORES=6 \
      ${LEVELS:+environment-variables.LEVELS="$LEVELS"} \
      environment-variables.BUCKET="$BUCKET" environment-variables.PREFIX="$prefix/4.4_stressor_effects" \
      -o json | python3 -I -c 'import json,sys;d=json.load(sys.stdin);print(d["job_runs"][0]["id"])')
    printf '%s\t%s\t%s\n' "$b" "$mc" "$id" | tee -a "$log"
  done
done
echo "job run ids in $log; follow with: scw jobs run list -o human"
