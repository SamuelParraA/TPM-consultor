#!/usr/bin/env bash
set -Eeuo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/common.sh"
require_command salmon
[[ -f "${CDNA_FASTA}" ]] || die "falta ${CDNA_FASTA}; ejecute primero la etapa reference"

if [[ -f "${SALMON_INDEX}/info.json" ]]; then
  python "${SCRIPT_DIR}/validate_salmon_index.py" --info "${SALMON_INDEX}/info.json" --version "${SALMON_VERSION_REQUIRED}" --k "${SALMON_K}"
  echo "OK: el índice Salmon ya existe y es compatible."
  exit 0
fi
[[ ! -e "${SALMON_INDEX}" ]] || die "existe un índice incompleto en ${SALMON_INDEX}; revíselo y muévalo"
temporary="${SALMON_INDEX}.building.$$"
log "Construyendo índice Salmon k=${SALMON_K}"
salmon index -t "${CDNA_FASTA}" -i "${temporary}" -k "${SALMON_K}" -p "${THREADS}"
python "${SCRIPT_DIR}/validate_salmon_index.py" --info "${temporary}/info.json" --version "${SALMON_VERSION_REQUIRED}" --k "${SALMON_K}"
mv "${temporary}" "${SALMON_INDEX}"
echo "OK: índice construido en ${SALMON_INDEX}"
