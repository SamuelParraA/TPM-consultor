#!/usr/bin/env bash
# Peipers-RNAseq fue creado por Samuel Parra A., PhD.
# sa.parra@uandresbello.edu | https://orcid.org/0000-0002-9129-4133
# Copyright (c) 2026 Samuel Parra A., PhD. SPDX-License-Identifier: MIT
set -Eeuo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
COMMAND="${1:-help}"
if [[ $# -ge 2 ]]; then
  export PEIPERS_CONFIG="$2"
fi

usage() {
  cat <<'EOF'
Uso:
  bash run_pipeline.sh ETAPA [archivo_config]

Etapas:
  preflight   Comprueba WSL, programas, versiones, metadata y espacio.
  reference   Descarga ITAG4.0 y construye el índice Salmon.
  download    Descarga los 36 FASTQ con reanudación y verifica MD5.
  qc          Ejecuta FastQC y MultiQC sin modificar las lecturas.
  quant       Cuantifica las 18 muestras con Salmon 2.0.1.
  matrices    Construye matrices de TPM/NumReads y tabla tx2gene.
  bundle      Crea el paquete pequeño y auditable apto para GitHub.
  consultor   Construye el consultor standalone y su SQLite local.
  verify      Verifica el paquete generado.
  deseq2      Ejecuta el módulo opcional de expresión diferencial.
  all         Ejecuta desde preflight hasta consultor (descarga ~100 GB).
EOF
}

run_stage() {
  bash "${HERE}/scripts/$1.sh"
}

case "${COMMAND}" in
  preflight) run_stage preflight ;;
  reference) run_stage download_references; run_stage build_salmon_index ;;
  download) run_stage download_fastq ;;
  qc) run_stage run_fastqc ;;
  quant) run_stage quantify_salmon ;;
  matrices) run_stage collect_matrices ;;
  bundle) run_stage build_bundle ;;
  consultor) run_stage build_consultor ;;
  verify) run_stage verify_bundle ;;
  deseq2) run_stage run_deseq2 ;;
  all)
    run_stage preflight
    run_stage download_references
    run_stage build_salmon_index
    run_stage download_fastq
    run_stage run_fastqc
    run_stage quantify_salmon
    run_stage collect_matrices
    run_stage build_bundle
    run_stage verify_bundle
    run_stage build_consultor
    ;;
  help|-h|--help) usage ;;
  *) usage; echo; echo "Etapa desconocida: ${COMMAND}" >&2; exit 2 ;;
esac
