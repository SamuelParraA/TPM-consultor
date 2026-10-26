#!/usr/bin/env bash
set -Eeuo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/common.sh"
require_command Rscript
[[ -f "${MATRIX_DIR}/tx2gene.tsv" ]] || die "falta tx2gene.tsv; ejecute primero matrices"
if [[ -f "${DESEQ2_DIR}/results_names.txt" ]]; then
  echo "DESeq2 ya fue ejecutado en ${DESEQ2_DIR}; no se sobrescribirá."
  exit 0
fi
mkdir -p "${DESEQ2_DIR}"
Rscript "${PIPELINE_ROOT}/optional_deseq2/run_deseq2.R" "${QUANT_DIR}" "${SAMPLE_METADATA}" "${MATRIX_DIR}/tx2gene.tsv" "${DESEQ2_DIR}"
