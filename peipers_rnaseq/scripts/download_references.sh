#!/usr/bin/env bash
set -Eeuo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/common.sh"
require_command wget
require_command sha256sum
ensure_work_dirs
BASE_URL="https://solgenomics.net/ftp/tomato_genome/annotation/ITAG4.0_release"

expected_hash() {
  awk -v name="$1" '$2 == name {print $1}' "${REFERENCE_CHECKSUMS}"
}

download_one() {
  local name="$1" destination="${REFERENCE_DIR}/$1" partial="${REFERENCE_DIR}/$1.part"
  local expected actual
  expected="$(expected_hash "${name}")"
  [[ -n "${expected}" ]] || die "no hay SHA-256 registrado para ${name}"
  if [[ -f "${destination}" ]]; then
    actual="$(sha256sum "${destination}" | awk '{print $1}')"
    [[ "${actual}" == "${expected}" ]] || die "${destination} existe pero su SHA-256 no coincide"
    echo "OK existente: ${name}"
    return
  fi
  log "Descargando ${name} desde Sol Genomics Network"
  wget --continue --tries=20 --timeout=60 --output-document="${partial}" "${BASE_URL}/${name}"
  actual="$(sha256sum "${partial}" | awk '{print $1}')"
  [[ "${actual}" == "${expected}" ]] || die "SHA-256 incorrecto en ${name}; no se instalará"
  mv "${partial}" "${destination}"
}

download_one "ITAG4.0_cDNA.fasta"
download_one "ITAG4.0_gene_models.gff"
python "${SCRIPT_DIR}/validate_reference_pair.py" --cdna "${CDNA_FASTA}" --gff "${GFF_FILE}"
log "Las referencias coinciden bit a bit con TPM-consultor"
