#!/usr/bin/env bash
set -Eeuo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/common.sh"
require_command python
if [[ -f "${BUNDLE_DIR}/dataset_manifest.json" ]]; then
  echo "El paquete ya existe; se verificará sin sobrescribirlo."
  python "${SCRIPT_DIR}/verify_bundle.py" --bundle "${BUNDLE_DIR}"
  exit 0
fi
python "${SCRIPT_DIR}/build_dataset_bundle.py" --project "${PROJECT_ID}" --quant-dir "${QUANT_DIR}" --matrix-dir "${MATRIX_DIR}" --multiqc-dir "${MULTIQC_DIR}" --runs "${RUN_MANIFEST}" --samples "${SAMPLE_METADATA}" --reference-checksums "${REFERENCE_CHECKSUMS}" --output "${BUNDLE_DIR}"
