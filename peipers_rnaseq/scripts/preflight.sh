#!/usr/bin/env bash
set -Eeuo pipefail
# shellcheck source=common.sh
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/common.sh"

ensure_work_dirs
log "Comprobación previa de ${PROJECT_ID}"

if grep -qi microsoft /proc/version 2>/dev/null; then
  echo "OK: se detectó Ubuntu/WSL."
else
  echo "AVISO: no se detectó WSL. El pipeline también funciona en Linux nativo."
fi

case "${REPO_ROOT}" in
  /mnt/c/*|/mnt/d/*|/mnt/e/*)
    echo "AVISO: el repositorio está en un disco montado desde Windows (${REPO_ROOT})."
    echo "       Para mejor rendimiento, clone el código en ~/TPM-consultor."
    ;;
esac
case "${WORK_DIR}" in
  *OneDrive*|*onedrive*) die "WORK_DIR está dentro de OneDrive. Elija una carpeta local no sincronizada." ;;
esac

for command in python salmon wget md5sum sha256sum fastqc multiqc; do
  require_command "${command}"
done

salmon_version="$(salmon --version 2>&1 | grep -Eo '[0-9]+\.[0-9]+\.[0-9]+' | head -n 1 || true)"
[[ "${salmon_version}" == "${SALMON_VERSION_REQUIRED}" ]] ||
  die "se requiere Salmon ${SALMON_VERSION_REQUIRED}; se encontró '${salmon_version:-desconocido}'."
echo "OK: Salmon ${salmon_version}."

python "${SCRIPT_DIR}/validate_project.py" --runs "${RUN_MANIFEST}" --samples "${SAMPLE_METADATA}" --expected-project "${PROJECT_ID}"

available_kb="$(df -Pk "${WORK_DIR}" | awk 'NR==2 {print $4}')"
available_gb="$((available_kb / 1024 / 1024))"
echo "Espacio libre en WORK_DIR: ${available_gb} GB."
if (( available_gb < MIN_FREE_GB )); then
  die "hay menos de ${MIN_FREE_GB} GB libres. Cambie WORK_DIR a un disco con más espacio."
fi
echo "Hilos configurados: ${THREADS}."
echo "Directorio de trabajo: ${WORK_DIR}"
echo "OK: el equipo está listo para comenzar."
