#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
REPO_ROOT="$(cd "${PIPELINE_ROOT}/.." && pwd)"
CONFIG_FILE="${PEIPERS_CONFIG:-${PIPELINE_ROOT}/config/PRJNA1347747.env}"

if [[ ! -f "${CONFIG_FILE}" ]]; then
  echo "ERROR: no existe la configuración: ${CONFIG_FILE}" >&2
  exit 2
fi
# shellcheck disable=SC1090
source "${CONFIG_FILE}"

RUN_MANIFEST="${PIPELINE_ROOT}/manifests/${PROJECT_ID}_runs.tsv"
SAMPLE_METADATA="${PIPELINE_ROOT}/manifests/${PROJECT_ID}_samples.tsv"
REFERENCE_CHECKSUMS="${PIPELINE_ROOT}/manifests/reference.sha256"
REFERENCE_DIR="${WORK_DIR}/reference/${REFERENCE_RELEASE}"
CDNA_FASTA="${REFERENCE_DIR}/ITAG4.0_cDNA.fasta"
GFF_FILE="${REFERENCE_DIR}/ITAG4.0_gene_models.gff"
SALMON_INDEX="${REFERENCE_DIR}/salmon_index_k${SALMON_K}"
FASTQ_DIR="${WORK_DIR}/fastq"
FASTQC_DIR="${WORK_DIR}/qc/fastqc"
MULTIQC_DIR="${WORK_DIR}/qc/multiqc"
QUANT_DIR="${WORK_DIR}/quant"
MATRIX_DIR="${WORK_DIR}/matrices"
BUNDLE_DIR="${WORK_DIR}/bundle/${PROJECT_ID}"
CONSULTOR_DIR="${WORK_DIR}/consultor_standalone/${PROJECT_ID}"
DESEQ2_DIR="${WORK_DIR}/deseq2_results"
SAMPLE_COUNT="$(awk 'END {print NR - 1}' "${RUN_MANIFEST}")"
EXPECTED_FASTQ_COUNT="$((SAMPLE_COUNT * 2))"

log() {
  printf '\n[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"
}

die() {
  echo "ERROR: $*" >&2
  exit 1
}

require_command() {
  command -v "$1" >/dev/null 2>&1 || die "falta el comando '$1'. Activa el ambiente Conda indicado en la guía."
}

ensure_work_dirs() {
  mkdir -p "${WORK_DIR}" "${REFERENCE_DIR}" "${FASTQ_DIR}" "${FASTQC_DIR}" \
    "${MULTIQC_DIR}" "${QUANT_DIR}" "${MATRIX_DIR}"
}
