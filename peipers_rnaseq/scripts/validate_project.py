#!/usr/bin/env python3
"""Valida la metadata antes de descargar cerca de 100 GB."""

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path

MD5 = re.compile(r"^[0-9a-f]{32}$")
RUN = re.compile(r"^SRR\d+$")


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise SystemExit(f"ERROR: no existe {path}")
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def require_columns(rows: list[dict[str, str]], columns: set[str], label: str) -> None:
    if not rows:
        raise SystemExit(f"ERROR: {label} está vacío")
    missing = columns - set(rows[0])
    if missing:
        raise SystemExit(f"ERROR: faltan columnas en {label}: {', '.join(sorted(missing))}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runs", type=Path, required=True)
    parser.add_argument("--samples", type=Path, required=True)
    parser.add_argument("--expected-project", required=True)
    args = parser.parse_args()

    runs = read_tsv(args.runs)
    samples = read_tsv(args.samples)
    require_columns(
        runs,
        {
            "project_id", "run_accession", "sample_accession", "experiment_accession",
            "genotype", "treatment", "replicate", "library_layout", "read_length",
            "read_count", "base_count", "fastq_1_url", "fastq_1_md5", "fastq_1_bytes",
            "fastq_2_url", "fastq_2_md5", "fastq_2_bytes", "order",
        },
        "manifiesto de corridas",
    )
    require_columns(
        samples,
        {"SampleID", "timepoint", "replicate", "genotype", "treatment", "order"},
        "metadata de muestras",
    )

    errors: list[str] = []
    accessions = [row["run_accession"] for row in runs]
    if args.expected_project == "PRJNA1347747" and len(runs) != 18:
        errors.append(f"se esperaban 18 corridas y hay {len(runs)}")
    duplicates = [key for key, count in Counter(accessions).items() if count > 1]
    if duplicates:
        errors.append(f"SRR duplicados: {duplicates}")

    conditions: dict[tuple[str, str], set[str]] = defaultdict(set)
    total_bytes = 0
    for number, row in enumerate(runs, 2):
        run = row["run_accession"]
        if row["project_id"] != args.expected_project:
            errors.append(f"línea {number}: BioProject incorrecto")
        if not RUN.fullmatch(run):
            errors.append(f"línea {number}: SRR inválido: {run}")
        if row["library_layout"] != "PAIRED":
            errors.append(f"{run}: se requiere library_layout=PAIRED")
        for mate in ("1", "2"):
            url = row[f"fastq_{mate}_url"]
            if not url.startswith("https://ftp.sra.ebi.ac.uk/") or not url.endswith(f"_{mate}.fastq.gz"):
                errors.append(f"{run}: URL del mate {mate} no es la esperada")
            if not MD5.fullmatch(row[f"fastq_{mate}_md5"]):
                errors.append(f"{run}: MD5 inválido para mate {mate}")
            try:
                total_bytes += int(row[f"fastq_{mate}_bytes"])
            except ValueError:
                errors.append(f"{run}: tamaño inválido para mate {mate}")
        try:
            expected_bases = int(row["read_count"]) * int(row["read_length"]) * 2
            if expected_bases != int(row["base_count"]):
                errors.append(f"{run}: base_count no coincide con pares × 2 × longitud")
        except ValueError:
            errors.append(f"{run}: read_count/read_length/base_count no son enteros")
        conditions[(row["genotype"], row["treatment"])].add(row["replicate"])

    if args.expected_project == "PRJNA1347747":
        expected_conditions = {
            (genotype, treatment): {"R1", "R2", "R3"}
            for genotype in ("WT", "epi", "ACCD")
            for treatment in ("Control", "Salt")
        }
        if conditions != expected_conditions:
            errors.append(f"diseño experimental incompleto: {dict(conditions)}")
    elif any("REVISAR" in str(value) for row in runs for value in row.values()):
        errors.append("el manifiesto genérico todavía contiene campos REVISAR")

    sample_ids = [row["SampleID"] for row in samples]
    if set(sample_ids) != set(accessions):
        errors.append("los SampleID no coinciden exactamente con los SRR del manifiesto")
    if len(sample_ids) != len(set(sample_ids)):
        errors.append("hay SampleID duplicados")
    try:
        orders = sorted(int(row["order"]) for row in samples)
    except ValueError:
        errors.append("la columna order contiene valores no enteros")
        orders = []
    if orders != list(range(1, len(samples) + 1)):
        errors.append("la columna order de muestras debe ser consecutiva desde 1")

    if errors:
        raise SystemExit("ERROR de metadata:\n- " + "\n- ".join(errors))
    print(f"OK: {len(runs)} corridas y metadata coherente.")
    print(f"Descarga FASTQ esperada: {total_bytes / 1024**3:.1f} GiB.")


if __name__ == "__main__":
    main()
