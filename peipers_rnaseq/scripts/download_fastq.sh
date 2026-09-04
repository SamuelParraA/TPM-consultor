#!/usr/bin/env bash
set -Eeuo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/common.sh"
require_command wget
require_command md5sum
ensure_work_dirs
python "${SCRIPT_DIR}/validate_project.py" --runs "${RUN_MANIFEST}" --samples "${SAMPLE_METADATA}" --expected-project "${PROJECT_ID}"

download_mate() {
  local run="$1" mate="$2" url="$3" expected_md5="$4" expected_bytes="$5"
  local destination="${FASTQ_DIR}/${run}_${mate}.fastq.gz" partial="${FASTQ_DIR}/${run}_${mate}.fastq.gz.part"
  local actual_md5 actual_bytes
  if [[ -f "${destination}" ]]; then
    actual_md5="$(md5sum "${destination}" | awk '{print $1}')"
    actual_bytes="$(stat -c '%s' "${destination}")"
    [[ "${actual_md5}" == "${expected_md5}" && "${actual_bytes}" == "${expected_bytes}" ]] ||
      die "${destination} existe pero no coincide con MD5/tamaño; no se sobrescribirá"
    echo "  OK existente: $(basename "${destination}")"
    return
  fi
  echo "  Descargando mate ${mate}: ${run}"
  wget --continue --tries=30 --timeout=60 --retry-connrefused --output-document="${partial}" "${url}"
  actual_md5="$(md5sum "${partial}" | awk '{print $1}')"
  actual_bytes="$(stat -c '%s' "${partial}")"
  [[ "${actual_md5}" == "${expected_md5}" ]] || die "MD5 incorrecto para ${run}, mate ${mate}"
  [[ "${actual_bytes}" == "${expected_bytes}" ]] || die "tamaño incorrecto para ${run}, mate ${mate}"
  mv "${partial}" "${destination}"
}

log "Descarga reanudable de FASTQ desde ENA"
{
  read -r _header
  while IFS=$'\t' read -r project run sample experiment sample_alias experiment_alias genotype treatment salinity replicate tissue layout platform read_length read_count base_count url1 md51 bytes1 url2 md52 bytes2 order; do
    [[ -n "${run}" ]] || continue
    echo "Muestra ${order}/${SAMPLE_COUNT}: ${run} (${genotype}, ${treatment}, ${replicate})"
    download_mate "${run}" 1 "${url1}" "${md51}" "${bytes1}"
    download_mate "${run}" 2 "${url2}" "${md52}" "${bytes2}"
  done
} < "${RUN_MANIFEST}"
touch "${FASTQ_DIR}/.download_complete"
log "Descarga completa: ${EXPECTED_FASTQ_COUNT} archivos verificados por MD5 y tamaño"
