#!/usr/bin/env bash
# Copy the outputs of a submission from Object Storage to outputs/<prefix>/ (keys from the scw config).
#   bash Input-EFSA/scaleway/fetch.sh <prefix>
set -euo pipefail
prefix="${1:?prefix}"
export RCLONE_CONFIG_SCW_TYPE=s3 RCLONE_CONFIG_SCW_PROVIDER=Scaleway RCLONE_CONFIG_SCW_REGION=fr-par \
       RCLONE_CONFIG_SCW_ENDPOINT=s3.fr-par.scw.cloud
RCLONE_CONFIG_SCW_ACCESS_KEY_ID="$(scw config get access-key)"
RCLONE_CONFIG_SCW_SECRET_ACCESS_KEY="$(scw config get secret-key)"
export RCLONE_CONFIG_SCW_ACCESS_KEY_ID RCLONE_CONFIG_SCW_SECRET_ACCESS_KEY
rclone copy --progress "scw:${BUCKET:?set BUCKET to your Object Storage bucket}/$prefix" "outputs/$prefix"
