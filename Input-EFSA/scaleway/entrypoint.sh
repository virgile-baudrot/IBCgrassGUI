#!/bin/sh
# Run one simulation script and copy its outputs to Object Storage.
# Environment:
#   SCRIPT  run_4.4_stressor_effects.R | run_4.3_control_design.R | run_4.1_single_species.R
#   ARGS    script arguments before out_dir, e.g. "BOR_I1 3" (4.4), "BOR_I1 all" (4.3), "all" (4.1)
#   NCORES  cores used (default: all)
#   BUCKET, PREFIX  destination s3://BUCKET/PREFIX/ (if BUCKET is empty, nothing is copied)
#   AWS_ACCESS_KEY_ID, AWS_SECRET_ACCESS_KEY  Object Storage keys (job secrets)
set -u
cd "$IBC_ROOT"
OUT=/tmp/out
mkdir -p "$OUT"
# shellcheck disable=SC2086
LOG="$OUT/job_$(hostname).log"
Rscript "Input-EFSA/${SCRIPT:-run_4.4_stressor_effects.R}" ${ARGS:-} "$OUT" "${NCORES:-$(nproc)}" > "$LOG" 2>&1
status=$?
cat "$LOG"
if [ -n "${BUCKET:-}" ]; then
  export RCLONE_CONFIG_SCW_TYPE=s3 RCLONE_CONFIG_SCW_PROVIDER=Scaleway RCLONE_CONFIG_SCW_ENV_AUTH=true \
         RCLONE_CONFIG_SCW_REGION=fr-par RCLONE_CONFIG_SCW_ENDPOINT=s3.fr-par.scw.cloud
  rclone copy "$OUT" "scw:$BUCKET/${PREFIX:-test}/" || status=1
fi
exit "$status"
