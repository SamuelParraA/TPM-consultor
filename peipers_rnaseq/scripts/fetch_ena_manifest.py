#!/usr/bin/env python3
"""Genera un manifiesto preliminar para reutilizar el pipeline con otro BioProject."""

import argparse
import csv
import io
import re
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

FIELDS = [
    "run_accession", "sample_accession", "experiment_accession", "sample_alias",
    "experiment_alias", "library_layout", "instrument_model", "nominal_length",
    "read_count", "base_count", "fastq_ftp", "fastq_md5", "fastq_bytes",
]
OUTPUT_FIELDS = [
    "project_id", "run_accession", "sample_accession", "experiment_accession",
    "sample_alias", "experiment_alias", "genotype", "treatment", "salinity_mM",
    "replicate", "tissue", "library_layout", "platform", "read_length", "read_count",
    "base_count", "fastq_1_url", "fastq_1_md5", "fastq_1_bytes", "fastq_2_url",
    "fastq_2_md5", "fastq_2_bytes", "order",
]


def infer_replicate(alias: str) -> str:
    match = re.search(r"(?:^|[-_])R([1-9])(?:$|[-_])", alias, re.IGNORECASE)
    return f"R{match.group(1)}" if match else "REVISAR"


def infer_treatment(alias: str):
    text = alias.casefold()
    if "salt" in text or "sal" in text or "nacl" in text:
        return "Salt", "REVISAR"
    if "control" in text or re.search(r"(?:^|[-_])c(?:$|[-_])", text):
        return "Control", "0"
    return "REVISAR", "REVISAR"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    query = urlencode({
        "accession": args.project,
        "result": "read_run",
        "fields": ",".join(FIELDS),
        "format": "tsv",
    })
    request = Request(
        "https://www.ebi.ac.uk/ena/portal/api/filereport?" + query,
        headers={"User-Agent": "Peipers-RNAseq-teaching/1.0"},
    )
    with urlopen(request, timeout=120) as response:
        rows = list(csv.DictReader(io.StringIO(response.read().decode("utf-8")), delimiter="\t"))
    if not rows:
        raise SystemExit(f"ERROR: ENA no devolvió corridas para {args.project}")

    output = []
    for order, row in enumerate(rows, 1):
        urls = row["fastq_ftp"].split(";")
        md5s = row["fastq_md5"].split(";")
        sizes = row["fastq_bytes"].split(";")
        if row["library_layout"] != "PAIRED" or not (len(urls) == len(md5s) == len(sizes) == 2):
            raise SystemExit(f"ERROR: {row['run_accession']} no es paired-end con dos FASTQ")
        alias = row["experiment_alias"]
        treatment, salinity = infer_treatment(alias)
        read_length = row.get("nominal_length", "")
        if not read_length and row.get("read_count") and row.get("base_count"):
            read_length = str(int(row["base_count"]) // (2 * int(row["read_count"])))
        output.append({
            "project_id": args.project,
            "run_accession": row["run_accession"],
            "sample_accession": row["sample_accession"],
            "experiment_accession": row["experiment_accession"],
            "sample_alias": row["sample_alias"],
            "experiment_alias": alias,
            "genotype": "REVISAR",
            "treatment": treatment,
            "salinity_mM": salinity,
            "replicate": infer_replicate(alias),
            "tissue": "REVISAR",
            "library_layout": row["library_layout"],
            "platform": row["instrument_model"],
            "read_length": read_length,
            "read_count": row["read_count"],
            "base_count": row["base_count"],
            "fastq_1_url": "https://" + urls[0].removeprefix("https://"),
            "fastq_1_md5": md5s[0],
            "fastq_1_bytes": sizes[0],
            "fastq_2_url": "https://" + urls[1].removeprefix("https://"),
            "fastq_2_md5": md5s[1],
            "fastq_2_bytes": sizes[1],
            "order": str(order),
        })
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output)
    print(f"Borrador: {args.output} ({len(output)} corridas).")
    print("IMPORTANTE: reemplace cada REVISAR usando el artículo y los BioSamples.")


if __name__ == "__main__":
    main()
