#!/usr/bin/env python3
import argparse
import hashlib
import json
from pathlib import Path

FORBIDDEN_SUFFIXES = (".fastq", ".fastq.gz", ".fq", ".fq.gz", ".sra")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


parser = argparse.ArgumentParser()
parser.add_argument("--bundle", type=Path, required=True)
args = parser.parse_args()
manifest_path = args.bundle / "dataset_manifest.json"
if not manifest_path.is_file():
    raise SystemExit(f"ERROR: falta {manifest_path}")
manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
errors = []
checksum_path = args.bundle / "SHA256SUMS"
if checksum_path.is_file():
    checksum_members = set()
    for line_number, line in enumerate(checksum_path.read_text(encoding="utf-8").splitlines(), start=1):
        try:
            expected_hash, relative_name = line.split("  ", 1)
        except ValueError:
            errors.append(f"SHA256SUMS, linea {line_number}: formato invalido")
            continue
        relative = Path(relative_name)
        if relative.is_absolute() or ".." in relative.parts:
            errors.append(f"SHA256SUMS, linea {line_number}: ruta insegura {relative_name}")
            continue
        checksum_members.add(relative.as_posix())
        path = args.bundle / relative
        if not path.is_file():
            errors.append(f"SHA256SUMS referencia un archivo ausente: {relative_name}")
        elif sha256(path) != expected_hash:
            errors.append(f"SHA256SUMS no coincide: {relative_name}")
    actual_members = {
        path.relative_to(args.bundle).as_posix()
        for path in args.bundle.rglob("*")
        if path.is_file() and path.name != "SHA256SUMS"
    }
    if checksum_members != actual_members:
        missing = sorted(actual_members - checksum_members)
        extra = sorted(checksum_members - actual_members)
        if missing:
            errors.append("archivos sin entrada en SHA256SUMS: " + ", ".join(missing))
        if extra:
            errors.append("entradas sobrantes en SHA256SUMS: " + ", ".join(extra))
for item in manifest.get("files", []):
    path = args.bundle / item["path"]
    if not path.is_file():
        errors.append(f"falta {item['path']}")
        continue
    if path.stat().st_size != item["bytes"]:
        errors.append(f"tamaño cambió: {item['path']}")
    if sha256(path) != item["sha256"]:
        errors.append(f"SHA-256 cambió: {item['path']}")
forbidden = [
    path.relative_to(args.bundle).as_posix()
    for path in args.bundle.rglob("*")
    if path.is_file() and path.name.lower().endswith(FORBIDDEN_SUFFIXES)
]
if forbidden:
    errors.append("lecturas crudas prohibidas: " + ", ".join(forbidden))
quant_files = list((args.bundle / "quant").glob("SRR*/quant.sf"))
if len(quant_files) != int(manifest.get("sample_count", -1)):
    errors.append(f"hay {len(quant_files)} quant.sf; el manifiesto declara {manifest.get('sample_count')}")
for required in (
    "metadata/run_manifest.tsv", "metadata/sample_metadata.tsv", "metadata/reference.sha256",
    "matrices/Gene_TPM.tsv", "matrices/Gene_NumReads.tsv", "matrices/tx2gene.tsv", "SHA256SUMS",
):
    if not (args.bundle / required).is_file():
        errors.append(f"falta {required}")
if manifest.get("raw_reads_included") is not False:
    errors.append("raw_reads_included debe ser false")
if errors:
    raise SystemExit("ERROR de paquete:\n- " + "\n- ".join(errors))
size = sum(path.stat().st_size for path in args.bundle.rglob("*") if path.is_file())
print(f"OK: paquete {manifest['project_id']} íntegro; {len(quant_files)} muestras; {size / 1024**2:.1f} MiB; sin FASTQ/SRA.")
