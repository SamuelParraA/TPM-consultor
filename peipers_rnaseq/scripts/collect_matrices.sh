#!/usr/bin/env bash
set -Eeuo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/common.sh"
require_command python
ensure_work_dirs
python "${SCRIPT_DIR}/collect_expression.py" --quant-dir "${QUANT_DIR}" --gff "${GFF_FILE}" --samples "${SAMPLE_METADATA}" --output "${MATRIX_DIR}"
