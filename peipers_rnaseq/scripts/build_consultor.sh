#!/usr/bin/env bash
set -Eeuo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/common.sh"
require_command python
if [[ -f "${CONSULTOR_DIR}/rnaseq_index.sqlite" ]]; then
  python "${SCRIPT_DIR}/verify_consultor.py" --directory "${CONSULTOR_DIR}" --experiment "${EXPERIMENT_ID}" --samples "${SAMPLE_COUNT}"
  echo "El consultor ya existe y no fue sobrescrito."
  exit 0
fi
python "${SCRIPT_DIR}/build_standalone_consultor.py" --repo-root "${REPO_ROOT}" --quant-dir "${QUANT_DIR}" --metadata "${SAMPLE_METADATA}" --gff "${GFF_FILE}" --experiment-id "${EXPERIMENT_ID}" --experiment-name "${EXPERIMENT_NAME}" --sample-count "${SAMPLE_COUNT}" --output "${CONSULTOR_DIR}"
python "${SCRIPT_DIR}/verify_consultor.py" --directory "${CONSULTOR_DIR}" --experiment "${EXPERIMENT_ID}" --samples "${SAMPLE_COUNT}"
echo "Para abrirlo: cd '${CONSULTOR_DIR}' && python app.py"
