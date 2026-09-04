#!/usr/bin/env python3
"""Crea un paquete auditable que contiene resultados, pero nunca lecturas crudas."""

import argparse
import csv
import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path

PIPELINE_VERSION = "1.0.0"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def copy_required(source: Path, destination: Path) -> None:
    if not source.is_file():
        raise SystemExit(f"ERROR: falta archivo requerido: {source}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", required=True)
    parser.add_argument("--quant-dir", type=Path, required=True)
    parser.add_argument("--matrix-dir", type=Path, required=True)
    parser.add_argument("--multiqc-dir", type=Path, required=True)
    parser.add_argument("--runs", type=Path, required=True)
    parser.add_argument("--samples", type=Path, required=True)
    parser.add_argument("--reference-checksums", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise SystemExit(f"ERROR: {args.output} ya existe; no se sobrescribirá")
    temporary = args.output.with_name(args.output.name + ".building")
    if temporary.exists():
        raise SystemExit(f"ERROR: quedó un paquete incompleto en {temporary}; revíselo y muévalo")
    temporary.mkdir(parents=True)

    copy_required(args.runs, temporary / "metadata" / "run_manifest.tsv")
    copy_required(args.samples, temporary / "metadata" / "sample_metadata.tsv")
    copy_required(args.reference_checksums, temporary / "metadata" / "reference.sha256")
    with args.runs.open(encoding="utf-8-sig", newline="") as handle:
        runs = list(csv.DictReader(handle, delimiter="\t"))
    for row in runs:
        run = row["run_accession"]
        source = args.quant_dir / run
        target = temporary / "quant" / run
        copy_required(source / "quant.sf", target / "quant.sf")
        copy_required(source / "cmd_info.json", target / "cmd_info.json")
        copy_required(source / "aux_info" / "meta_info.json", target / "aux_info" / "meta_info.json")
        for optional in ("lib_format_counts.json", "logs/salmon_quant.log"):
            if (source / optional).is_file():
                copy_required(source / optional, target / optional)

    for name in (
        "Transcript_TPM.tsv", "Transcript_NumReads.tsv", "Gene_TPM.tsv",
        "Gene_NumReads.tsv", "tx2gene.tsv", "salmon_qc.tsv", "matrix_summary.json",
    ):
        copy_required(args.matrix_dir / name, temporary / "matrices" / name)
    report = args.multiqc_dir / "multiqc_report.html"
    if report.is_file():
        copy_required(report, temporary / "qc" / "multiqc_report.html")
        data_dir = args.multiqc_dir / "multiqc_data"
        if data_dir.is_dir():
            shutil.copytree(data_dir, temporary / "qc" / "multiqc_data")

    readme = f"""# Dataset {args.project}

Paquete producido por Peipers-RNAseq docente v{PIPELINE_VERSION}.

Incluye metadata curada, salidas Salmon 2.0.1, matrices de expresión y QC.
No incluye FASTQ, archivos SRA ni el índice Salmon.
"""
    (temporary / "README.md").write_text(readme, encoding="utf-8")
    manifest = {
        "schema_version": 1,
        "pipeline": "Peipers-RNAseq docente",
        "pipeline_version": PIPELINE_VERSION,
        "project_id": args.project,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "reference": "ITAG4.0",
        "salmon_version": "2.0.1",
        "sample_count": len(runs),
        "raw_reads_included": False,
        "files": [],
    }
    manifest_path = temporary / "dataset_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    tracked = sorted(
        path for path in temporary.rglob("*")
        if path.is_file() and path.name not in {"dataset_manifest.json", "SHA256SUMS"}
    )
    manifest["files"] = [
        {"path": path.relative_to(temporary).as_posix(), "bytes": path.stat().st_size, "sha256": sha256(path)}
        for path in tracked
    ]
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    checksummed = sorted(path for path in temporary.rglob("*") if path.is_file() and path.name != "SHA256SUMS")
    lines = [f"{sha256(path)}  {path.relative_to(temporary).as_posix()}" for path in checksummed]
    (temporary / "SHA256SUMS").write_text("\n".join(lines) + "\n", encoding="utf-8")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary.rename(args.output)
    print(f"OK: paquete GitHub creado en {args.output}")


if __name__ == "__main__":
    main()
