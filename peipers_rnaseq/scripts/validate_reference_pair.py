#!/usr/bin/env python3
"""Comprueba coherencia entre el cDNA cuantificable y los modelos del GFF."""

import argparse
from pathlib import Path


def clean(value: str) -> str:
    return value.split(":", 1)[-1]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cdna", type=Path, required=True)
    parser.add_argument("--gff", type=Path, required=True)
    args = parser.parse_args()
    cdna = []
    with args.cdna.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                cdna.append(line[1:].split()[0])
    gff_transcripts = set()
    with args.gff.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] not in {"mRNA", "transcript"}:
                continue
            attrs = dict(item.split("=", 1) for item in fields[8].split(";") if "=" in item)
            transcript = clean(attrs.get("ID", attrs.get("Name", "")))
            if transcript:
                gff_transcripts.add(transcript)
    errors = []
    if len(cdna) != 34075 or len(set(cdna)) != 34075:
        errors.append(f"cDNA tiene {len(cdna)} encabezados y {len(set(cdna))} únicos")
    absent = sorted(set(cdna) - gff_transcripts)
    if absent:
        errors.append(f"{len(absent)} transcritos cDNA no aparecen en GFF")
    if errors:
        raise SystemExit("ERROR: " + "; ".join(errors))
    without_cdna = len(gff_transcripts - set(cdna))
    if without_cdna:
        raise SystemExit(f"ERROR: {without_cdna} modelos GFF carecen de encabezado cDNA")
    print("OK: 34.075 cDNA mapeables uno a uno al GFF.")
    print("Nota: Salmon colapsa 99 secuencias idénticas y cuantifica 33.976 referencias retenidas.")


if __name__ == "__main__":
    main()
