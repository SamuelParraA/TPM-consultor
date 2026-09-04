#!/usr/bin/env python3
"""Combina quant.sf y agrega transcritos ITAG4.0 a nivel de gen."""

from __future__ import annotations

import argparse
import csv
import json
from collections import OrderedDict
from pathlib import Path
from urllib.parse import unquote


def read_tsv(path: Path) -> list[dict[str, str]]:
    delimiter = "," if path.suffix.lower() == ".csv" else "\t"
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter=delimiter))
    for row in rows:
        row["SampleID"] = row.get("SampleID") or row.get("sample", "")
    return rows


def parse_attributes(text: str) -> dict[str, str]:
    result = {}
    for item in text.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = unquote(value)
    return result


def clean_feature_id(value: str) -> str:
    return value.split(":", 1)[-1]


def parse_gff(path: Path):
    genes = OrderedDict()
    transcripts = OrderedDict()
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            seqid, _, feature, start, end, _, strand, _, raw = fields
            attrs = parse_attributes(raw)
            if feature == "gene":
                gene_id = clean_feature_id(attrs.get("ID", attrs.get("Name", "")))
                if gene_id:
                    genes[gene_id] = {
                        "GeneID": gene_id,
                        "Alias": attrs.get("Alias", ""),
                        "Name": attrs.get("Name", gene_id),
                        "Description": attrs.get("Note", attrs.get("description", "")),
                        "Chromosome": seqid,
                        "Start": start,
                        "End": end,
                        "Strand": strand,
                    }
            elif feature in {"mRNA", "transcript"}:
                transcript_id = clean_feature_id(attrs.get("ID", attrs.get("Name", "")))
                gene_id = clean_feature_id(attrs.get("Parent", "").split(",", 1)[0])
                if transcript_id and gene_id:
                    transcripts[transcript_id] = gene_id
                    if gene_id in genes and not genes[gene_id]["Description"]:
                        genes[gene_id]["Description"] = attrs.get("Note", "")
    if len(transcripts) < 33976:
        raise SystemExit(f"ERROR: GFF produjo sólo {len(transcripts)} modelos de transcrito")
    return genes, transcripts


def format_value(value: float) -> str:
    return f"{value:.6f}".rstrip("0").rstrip(".") or "0"


def write_matrix(path, id_column, identifiers, sample_ids, values, annotations=None):
    annotation_columns = ["Alias", "Name", "Description", "Chromosome", "Start", "End", "Strand"] if annotations else []
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([id_column, *annotation_columns, *sample_ids])
        for index, identifier in enumerate(identifiers):
            prefix = [identifier]
            if annotations is not None:
                prefix.extend(annotations[identifier].get(column, "") for column in annotation_columns)
            writer.writerow([*prefix, *(format_value(values[sample][index]) for sample in sample_ids)])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quant-dir", type=Path, required=True)
    parser.add_argument("--gff", type=Path, required=True)
    parser.add_argument("--samples", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    metadata = sorted(read_tsv(args.samples), key=lambda row: int(row["order"]))
    sample_ids = [row["SampleID"] for row in metadata]
    genes, transcript_to_gene = parse_gff(args.gff)
    first_quant = args.quant_dir / sample_ids[0] / "quant.sf"
    if not first_quant.is_file():
        raise SystemExit(f"ERROR: falta cuantificación para {sample_ids[0]}")
    with first_quant.open(encoding="utf-8", newline="") as handle:
        transcript_ids = [row["Name"] for row in csv.DictReader(handle, delimiter="\t")]
    if len(transcript_ids) != 33976 or len(set(transcript_ids)) != 33976:
        raise SystemExit(f"ERROR: el primer quant.sf tiene {len(transcript_ids)} transcritos únicos esperados")
    missing_from_gff = [transcript for transcript in transcript_ids if transcript not in transcript_to_gene]
    if missing_from_gff:
        raise SystemExit(f"ERROR: {len(missing_from_gff)} transcritos cuantificados no están en el GFF")
    transcript_index = {transcript: index for index, transcript in enumerate(transcript_ids)}
    gene_ids = list(genes)
    gene_index = {gene: index for index, gene in enumerate(gene_ids)}

    transcript_tpm = {}
    transcript_reads = {}
    gene_tpm = {}
    gene_reads = {}
    qc_rows = []

    for position, sample in enumerate(sample_ids, 1):
        quant_path = args.quant_dir / sample / "quant.sf"
        meta_path = args.quant_dir / sample / "aux_info" / "meta_info.json"
        if not quant_path.is_file() or not meta_path.is_file():
            raise SystemExit(f"ERROR: falta cuantificación completa para {sample}")
        tpm_values = [0.0] * len(transcript_ids)
        read_values = [0.0] * len(transcript_ids)
        seen = set()
        with quant_path.open(encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                transcript = row["Name"]
                if transcript not in transcript_index:
                    raise SystemExit(f"ERROR: {sample} contiene transcrito ajeno a ITAG4.0: {transcript}")
                index = transcript_index[transcript]
                seen.add(transcript)
                tpm_values[index] = float(row["TPM"])
                read_values[index] = float(row["NumReads"])
        if len(seen) != len(transcript_ids):
            raise SystemExit(f"ERROR: {sample} contiene {len(seen)} de {len(transcript_ids)} transcritos")
        transcript_tpm[sample] = tpm_values
        transcript_reads[sample] = read_values

        sample_gene_tpm = [0.0] * len(gene_ids)
        sample_gene_reads = [0.0] * len(gene_ids)
        for transcript, index in transcript_index.items():
            gene = transcript_to_gene[transcript]
            if gene in gene_index:
                target = gene_index[gene]
                sample_gene_tpm[target] += tpm_values[index]
                sample_gene_reads[target] += read_values[index]
        gene_tpm[sample] = sample_gene_tpm
        gene_reads[sample] = sample_gene_reads

        info = json.loads(meta_path.read_text(encoding="utf-8"))
        qc_rows.append({
            "SampleID": sample,
            "SalmonVersion": info.get("salmon_version", ""),
            "LibraryType": ",".join(info.get("library_types", [])),
            "ProcessedFragments": info.get("num_processed", ""),
            "MappedFragments": info.get("num_mapped", ""),
            "MappingRatePercent": info.get("percent_mapped", ""),
            "ValidTargets": info.get("num_valid_targets", ""),
            "SeqBias": info.get("seq_bias_correct", False),
            "GCBias": info.get("gc_bias_correct", False),
        })
        print(f"[{position}/{len(sample_ids)}] {sample}")

    write_matrix(args.output / "Transcript_TPM.tsv", "Transcript", transcript_ids, sample_ids, transcript_tpm)
    write_matrix(args.output / "Transcript_NumReads.tsv", "Transcript", transcript_ids, sample_ids, transcript_reads)
    write_matrix(args.output / "Gene_TPM.tsv", "GeneID", gene_ids, sample_ids, gene_tpm, genes)
    write_matrix(args.output / "Gene_NumReads.tsv", "GeneID", gene_ids, sample_ids, gene_reads, genes)

    with (args.output / "tx2gene.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["TXNAME", "GENEID"])
        writer.writerows((transcript, transcript_to_gene[transcript]) for transcript in transcript_ids)
    with (args.output / "salmon_qc.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(qc_rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(qc_rows)

    summary = {
        "samples": len(sample_ids),
        "transcripts": len(transcript_ids),
        "genes": len(gene_ids),
        "reference": "ITAG4.0",
        "aggregation": "gene TPM and NumReads are sums across transcripts mapped by the ITAG4.0 GFF",
    }
    (args.output / "matrix_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"OK: matrices escritas en {args.output}")


if __name__ == "__main__":
    main()
