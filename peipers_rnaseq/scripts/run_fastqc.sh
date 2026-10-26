#!/usr/bin/env bash
set -Eeuo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/common.sh"
require_command fastqc
require_command multiqc
ensure_work_dirs

actual_fastq="$(find "${FASTQ_DIR}" -maxdepth 1 -type f -name '*.fastq.gz' | wc -l)"
[[ "${actual_fastq}" -eq "${EXPECTED_FASTQ_COUNT}" ]] || die "se esperaban ${EXPECTED_FASTQ_COUNT} FASTQ y hay ${actual_fastq}; complete la descarga"
if [[ "${RUN_FASTQC}" != "1" ]]; then
  echo "RUN_FASTQC=${RUN_FASTQC}; se omite QC por configuración."
  exit 0
fi
if [[ -f "${MULTIQC_DIR}/multiqc_report.html" ]]; then
  echo "OK: MultiQC ya existe. Para regenerarlo, mueva primero ${MULTIQC_DIR}."
  exit 0
fi
log "FastQC: inspección de calidad de ${EXPECTED_FASTQ_COUNT} FASTQ (no modifica lecturas)"
fastqc --threads "${THREADS}" --outdir "${FASTQC_DIR}" "${FASTQ_DIR}"/*.fastq.gz
log "MultiQC: resumen conjunto de todos los informes FastQC"
multiqc --force --outdir "${MULTIQC_DIR}" "${FASTQC_DIR}"
echo "Informe: ${MULTIQC_DIR}/multiqc_report.html"
