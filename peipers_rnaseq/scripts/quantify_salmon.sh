#!/usr/bin/env bash
set -Eeuo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/common.sh"
require_command salmon
ensure_work_dirs
[[ -f "${SALMON_INDEX}/info.json" ]] || die "falta el índice; ejecute bash run_pipeline.sh reference"
python "${SCRIPT_DIR}/validate_salmon_index.py" --info "${SALMON_INDEX}/info.json" --version "${SALMON_VERSION_REQUIRED}" --k "${SALMON_K}"

quantify_one() {
  local run="$1" order="$2"
  local mate1="${FASTQ_DIR}/${run}_1.fastq.gz" mate2="${FASTQ_DIR}/${run}_2.fastq.gz"
  local output="${QUANT_DIR}/${run}" temporary="${QUANT_DIR}/${run}.building.$$"
  [[ -f "${mate1}" && -f "${mate2}" ]] || die "faltan FASTQ para ${run}"
  if [[ -f "${output}/quant.sf" && -f "${output}/aux_info/meta_info.json" ]]; then
    python "${SCRIPT_DIR}/validate_quant.py" --quant-dir "${output}" --sample "${run}"
    return
  fi
  [[ ! -e "${output}" ]] || die "salida incompleta en ${output}; muévala antes de reintentar"
  echo "  [${order}/${SAMPLE_COUNT}] Cuantificando ${run}"
  # Compatibilidad: mismos parámetros históricos. No activar seqBias/gcBias
  # ni recortar lecturas por defecto.
  salmon quant -i "${SALMON_INDEX}" -l A -1 "${mate1}" -2 "${mate2}" -p "${THREADS}" -o "${temporary}"
  python "${SCRIPT_DIR}/validate_quant.py" --quant-dir "${temporary}" --sample "${run}"
  mv "${temporary}" "${output}"
}

log "Cuantificación Salmon compatible con TPM-consultor"
{
  read -r _header
  while IFS=$'\t' read -r project run sample experiment sample_alias experiment_alias genotype treatment salinity replicate tissue layout platform read_length read_count base_count url1 md51 bytes1 url2 md52 bytes2 order; do
    [[ -n "${run}" ]] || continue
    quantify_one "${run}" "${order}"
  done
} < "${RUN_MANIFEST}"
touch "${QUANT_DIR}/.quant_complete"
log "Las ${SAMPLE_COUNT} muestras tienen quant.sf y métricas Salmon válidas"
