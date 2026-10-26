#!/usr/bin/env python3
"""Plataforma local para explorar expresión transcriptómica Salmon + ITAG4.0."""

from __future__ import annotations

import argparse
import csv
from difflib import SequenceMatcher
import json
import os
import re
import sqlite3
import sys
from collections import defaultdict
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse

ROOT = Path(__file__).resolve().parent
DEFAULT_CONFIG = ROOT / "config.example.json"
DEFAULT_DB = Path(os.environ.get("DB_PATH", str(ROOT / "rnaseq_index.sqlite")))
LEGACY_QUANT = Path(os.environ.get("QUANT_DIR", r"I:\transcriptomica\salmon_timecourse"))
LEGACY_GFF = Path(os.environ.get("GFF_PATH", r"I:\transcriptomica\reference\ITAG4.0_gene_models.gff"))


def strip_prefix(value: str) -> str:
    return re.sub(r"^(gene|mRNA|transcript):", "", value or "")


def decode_attr(value: str) -> str:
    return (
        (value or "")
        .replace("%20", " ")
        .replace("%2C", ",")
        .replace("%3B", ";")
        .replace("%3A", ":")
    )


def attributes(text: str) -> dict[str, str]:
    result = {}
    for item in text.rstrip().split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = decode_attr(value)
    return result


def clean_name(value: str) -> str:
    value = (value or "").strip()
    return value or "No especificado"


def read_table(path: Path) -> list[dict[str, str]]:
    delimiter = "\t" if path.suffix.lower() in {".tsv", ".tab"} else ","
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return [{k.strip(): (v or "").strip() for k, v in row.items()} for row in csv.DictReader(handle, delimiter=delimiter)]


def load_config(path: Path) -> dict:
    if path.exists():
        with path.open(encoding="utf-8") as handle:
            config = json.load(handle)
    else:
        config = {}
    data_root = Path(os.environ.get("TRANSCRIPTOMICA_ROOT", config.get("data_root", r"I:\transcriptomica")))
    reference_dir = Path(config.get("reference_dir", str(data_root / "reference")))
    annotation_dir = Path(config.get("annotation_dir", str(data_root / "annotation")))
    return {
        "data_root": data_root,
        "reference_dir": Path(os.environ.get("REFERENCE_DIR", str(reference_dir))),
        "annotation_dir": Path(os.environ.get("ANNOTATION_DIR", str(annotation_dir))),
        "gff": Path(os.environ.get("GFF_PATH", config.get("gff", str(reference_dir / "ITAG4.0_gene_models.gff")))),
        "experiments": config.get("experiments", []),
    }


def db_has_column(con: sqlite3.Connection, table: str, column: str) -> bool:
    return any(row[1] == column for row in con.execute(f"PRAGMA table_info({table})"))


def discover_experiments(config: dict, fallback_quant: Path | None = None) -> list[dict]:
    configured = []
    for item in config.get("experiments") or []:
        root = Path(item["path"])
        configured.append({
            "experiment_id": item.get("id") or root.name,
            "name": item.get("name") or root.name.replace("_", " ").title(),
            "path": root,
            "metadata": Path(item.get("metadata", str(root / "sample_metadata.tsv"))),
        })
    if configured:
        return configured

    data_root = Path(config["data_root"])
    candidates = []
    if fallback_quant:
        candidates.append(fallback_quant)
    if data_root.exists():
        candidates.extend(path for path in data_root.iterdir() if path.is_dir())
    seen, experiments = set(), []
    for root in candidates:
        try:
            root = root.resolve()
        except OSError:
            continue
        if root in seen:
            continue
        seen.add(root)
        if list(root.glob("*/quant.sf")):
            experiments.append({
                "experiment_id": root.name,
                "name": root.name.replace("salmon_", "").replace("_", " ").title(),
                "path": root,
                "metadata": next((root / name for name in ("sample_metadata.tsv", "sample_metadata.csv") if (root / name).exists()), ROOT / "sample_metadata.csv"),
            })
    return experiments


def load_metadata(path: Path, sample_names: list[str]) -> list[dict[str, str]]:
    configured = {}
    if path.exists():
        for row in read_table(path):
            sample = row.get("SampleID") or row.get("sample") or row.get("SRR") or row.get("srr")
            if sample:
                configured[sample.strip()] = row
    rows = []
    for order, sample in enumerate(sample_names, 1):
        row = configured.get(sample, {})
        bio = row.get("Nombre biológico") or row.get("Nombre biologico") or row.get("biological_name") or row.get("name") or row.get("timepoint") or sample
        time = row.get("Tiempo") or row.get("time") or row.get("timepoint") or bio
        replicate = row.get("Réplica") or row.get("Replica") or row.get("replicate") or ""
        treatment = row.get("Tratamiento") or row.get("treatment") or ""
        group = row.get("Grupo") or row.get("group") or ""
        label = " · ".join(x for x in [bio, replicate, group, treatment] if x and x != sample)
        rows.append({
            "sample": sample,
            "label": label or bio or sample,
            "timepoint": clean_name(time),
            "replicate": replicate.strip(),
            "group_name": group.strip(),
            "treatment": treatment.strip(),
            "color": (row.get("Color sugerido") or row.get("color") or "").strip(),
            "comments": (row.get("Comentarios") or row.get("comments") or "").strip(),
            "order": int(row.get("order") or row.get("Orden") or order),
        })
    return sorted(rows, key=lambda row: (row["order"], row["sample"]))


def parse_gff(gff_path: Path) -> tuple[dict[str, dict], dict[str, str]]:
    if not gff_path.exists():
        raise FileNotFoundError(f"No se encontró el GFF: {gff_path}")
    genes: dict[str, dict] = {}
    transcript_to_gene: dict[str, str] = {}
    with gff_path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            seqid, _, feature, start, end, _, strand, _, raw_attrs = fields
            attrs = attributes(raw_attrs)
            length = int(end) - int(start) + 1
            if feature == "gene":
                gene_id = strip_prefix(attrs.get("ID", attrs.get("Name", "")))
                if gene_id:
                    dbxref = attrs.get("Dbxref", "")
                    genes[gene_id] = {
                        "gene_id": gene_id,
                        "alias": attrs.get("Alias", ""),
                        "name": attrs.get("Name", gene_id),
                        "description": attrs.get("Note", attrs.get("description", "")),
                        "go": attrs.get("Ontology_term", ""),
                        "interpro": ";".join(re.findall(r"InterPro:[^,;]+", dbxref)),
                        "pfam": ";".join(re.findall(r"Pfam:[^,;]+", dbxref)),
                        "dbxref": dbxref,
                        "seqid": seqid,
                        "start": int(start),
                        "end": int(end),
                        "strand": strand,
                        "exon_count": 0,
                        "gene_length": length,
                        "mrna_length": 0,
                        "transcripts": [],
                    }
            elif feature in {"mRNA", "transcript"}:
                transcript_id = strip_prefix(attrs.get("ID", attrs.get("Name", "")))
                gene_id = strip_prefix(attrs.get("Parent", "").split(",")[0])
                if transcript_id and gene_id:
                    transcript_to_gene[transcript_id] = gene_id
                    if gene_id in genes:
                        genes[gene_id]["transcripts"].append(transcript_id)
                        genes[gene_id]["mrna_length"] = max(genes[gene_id]["mrna_length"], length)
                        if attrs.get("Note") and not genes[gene_id]["description"]:
                            genes[gene_id]["description"] = attrs["Note"]
            elif feature == "exon":
                gene_id = strip_prefix(attrs.get("Parent", "").split(",")[0])
                gene_id = transcript_to_gene.get(gene_id, gene_id)
                if gene_id in genes:
                    genes[gene_id]["exon_count"] += 1
    return genes, transcript_to_gene


def load_annotation_tables(annotation_dir: Path) -> dict[str, list[dict[str, str]]]:
    annotations: dict[str, list[dict[str, str]]] = defaultdict(list)
    if not annotation_dir.exists():
        return annotations
    for table in sorted(annotation_dir.glob("*.tsv")) + sorted(annotation_dir.glob("*.csv")):
        for row in read_table(table):
            gene_id = row.get("ITAG") or row.get("gene_id") or row.get("GeneID") or row.get("Gene ID")
            if not gene_id:
                continue
            payload = {k: v for k, v in row.items() if k and v and k not in {"ITAG", "gene_id", "GeneID", "Gene ID"}}
            if payload:
                annotations[gene_id].append({"source": table.name, **payload})
    literature = annotation_dir / "literature"
    if literature.exists():
        for table in sorted(literature.glob("*.tsv")) + sorted(literature.glob("*.csv")):
            for row in read_table(table):
                gene_id = row.get("ITAG") or row.get("gene_id") or row.get("GeneID") or row.get("Gene ID")
                if gene_id:
                    payload = {k: v for k, v in row.items() if k and v and k not in {"ITAG", "gene_id", "GeneID", "Gene ID"}}
                    annotations[gene_id].append({"source": f"literature/{table.name}", **payload})
    return annotations


def build_index(config_path: Path, quant_dir: Path | None, gff_path: Path | None, db_path: Path) -> None:
    config = load_config(config_path)
    if gff_path:
        config["gff"] = gff_path
    experiments = discover_experiments(config, quant_dir)
    if not experiments:
        raise FileNotFoundError("No se encontraron experimentos con subcarpetas */quant.sf")

    print(f"Leyendo anotación oficial: {config['gff']}")
    genes, transcript_to_gene = parse_gff(Path(config["gff"]))
    curated_annotations = load_annotation_tables(Path(config["annotation_dir"]))

    db_path.unlink(missing_ok=True)
    con = sqlite3.connect(db_path)
    con.executescript("""
        CREATE TABLE experiments (
          experiment_id TEXT PRIMARY KEY, name TEXT, path TEXT, metadata_path TEXT
        );
        CREATE TABLE genes (
          gene_id TEXT PRIMARY KEY, alias TEXT, name TEXT, description TEXT, go TEXT,
          interpro TEXT, pfam TEXT, dbxref TEXT, seqid TEXT, start INTEGER, end INTEGER,
          strand TEXT, exon_count INTEGER, gene_length INTEGER, mrna_length INTEGER,
          transcripts TEXT, search_text TEXT
        );
        CREATE TABLE gene_annotations (
          gene_id TEXT, source TEXT, field TEXT, value TEXT
        );
        CREATE TABLE samples (
          experiment_id TEXT, sample TEXT, label TEXT, timepoint TEXT, replicate TEXT,
          group_name TEXT, treatment TEXT, color TEXT, comments TEXT, sort_order INTEGER,
          PRIMARY KEY (experiment_id, sample)
        );
        CREATE TABLE expression (
          experiment_id TEXT, gene_id TEXT, sample TEXT, tpm REAL, reads REAL,
          PRIMARY KEY (experiment_id, gene_id, sample)
        );
        CREATE INDEX idx_gene_search ON genes(search_text);
        CREATE INDEX idx_expr_gene ON expression(experiment_id, gene_id);
    """)

    gene_rows = []
    annotation_rows = []
    for gene in genes.values():
        extra = curated_annotations.get(gene["gene_id"], [])
        extra_text = " ".join(" ".join(item.values()) for item in extra)
        tx = ",".join(gene["transcripts"])
        search = " ".join([
            gene["gene_id"], gene["alias"], gene["name"], gene["description"], gene["go"],
            gene["interpro"], gene["pfam"], gene["dbxref"], tx, extra_text
        ]).lower()
        gene_rows.append((
            gene["gene_id"], gene["alias"], gene["name"], gene["description"], gene["go"],
            gene["interpro"], gene["pfam"], gene["dbxref"], gene["seqid"], gene["start"], gene["end"],
            gene["strand"], gene["exon_count"], gene["gene_length"], gene["mrna_length"], tx, search
        ))
        for item in extra:
            source = item.get("source", "annotation")
            for field, value in item.items():
                if field != "source":
                    annotation_rows.append((gene["gene_id"], source, field, value))
    con.executemany("INSERT INTO genes VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", gene_rows)
    con.executemany("INSERT INTO gene_annotations VALUES (?,?,?,?)", annotation_rows)

    for exp in experiments:
        quant_files = sorted(Path(exp["path"]).glob("*/quant.sf"))
        if not quant_files:
            continue
        exp_id = exp["experiment_id"]
        con.execute("INSERT INTO experiments VALUES (?,?,?,?)", (exp_id, exp["name"], str(exp["path"]), str(exp["metadata"])))
        metadata = load_metadata(Path(exp["metadata"]), [path.parent.name for path in quant_files])
        con.executemany("INSERT INTO samples VALUES (?,?,?,?,?,?,?,?,?,?)", [
            (exp_id, m["sample"], m["label"], m["timepoint"], m["replicate"], m["group_name"], m["treatment"], m["color"], m["comments"], m["order"])
            for m in metadata
        ])
        expression: dict[tuple[str, str], list[float]] = defaultdict(lambda: [0.0, 0.0])
        print(f"Indexando experimento: {exp['name']} ({len(quant_files)} muestras)")
        for index, quant_path in enumerate(quant_files, 1):
            sample = quant_path.parent.name
            print(f"  [{index}/{len(quant_files)}] {sample}")
            with quant_path.open(encoding="utf-8", newline="") as handle:
                for row in csv.DictReader(handle, delimiter="\t"):
                    transcript = row["Name"]
                    gene_id = transcript_to_gene.get(transcript)
                    if not gene_id:
                        candidate = transcript.rsplit(".", 1)[0]
                        gene_id = candidate if candidate in genes else None
                    if gene_id:
                        values = expression[(gene_id, sample)]
                        values[0] += float(row["TPM"])
                        values[1] += float(row["NumReads"])
        con.executemany("INSERT INTO expression VALUES (?,?,?,?,?)", [
            (exp_id, gene_id, sample, values[0], values[1])
            for (gene_id, sample), values in expression.items()
        ])
    con.commit()
    con.close()
    print(f"Índice listo: {len(genes):,} genes, {len(experiments)} experimento(s) -> {db_path}")


class Handler(SimpleHTTPRequestHandler):
    db_path: Path

    def __init__(self, *args, **kwargs):
        super().__init__(*args, directory=str(ROOT / "web"), **kwargs)

    def send_json(self, payload, status=200):
        body = json.dumps(payload, ensure_ascii=False).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    def connect(self) -> sqlite3.Connection:
        con = sqlite3.connect(self.db_path)
        con.row_factory = sqlite3.Row
        return con

    def is_multi_experiment(self, con: sqlite3.Connection) -> bool:
        return db_has_column(con, "samples", "experiment_id")

    def default_experiment(self, con: sqlite3.Connection) -> str:
        if self.is_multi_experiment(con):
            row = con.execute("SELECT experiment_id FROM experiments ORDER BY name LIMIT 1").fetchone()
            return row["experiment_id"] if row else ""
        return "default"

    def rank_gene(self, row, query: str) -> float:
        identifiers = [str(row[key] or "").lower() for key in row.keys() if key in {"gene_id", "alias", "name"}]
        description = str(row["description"] or "").lower() if "description" in row.keys() else ""
        if query in identifiers:
            return 1000
        if any(value.startswith(query) for value in identifiers if value):
            return 900
        if any(query in value for value in identifiers if value):
            return 800
        position = description.find(query)
        phrase_score = 500 - min(position, 200) if position >= 0 else 0
        words = re.findall(r"[a-z0-9_-]+", description)
        width = max(1, len(query.split())) + 2
        windows = [" ".join(words[i:i + width]) for i in range(len(words))]
        similarity = max((SequenceMatcher(None, query, value).ratio() for value in windows), default=0)
        return phrase_score + similarity * 100

    def gene_payload(self, con: sqlite3.Connection, gene, experiment_id: str) -> dict:
        item = dict(gene)
        item.pop("search_text", None)
        item["transcripts"] = item.get("transcripts", "").split(",") if item.get("transcripts") else []
        if self.is_multi_experiment(con):
            item["expression"] = [dict(row) for row in con.execute("""SELECT s.sample, s.label, s.timepoint, s.replicate,
                s.group_name, s.treatment, s.color, s.comments, COALESCE(e.tpm,0) AS tpm, COALESCE(e.reads,0) AS reads
                FROM samples s LEFT JOIN expression e
                  ON e.experiment_id=s.experiment_id AND e.sample=s.sample AND e.gene_id=?
                WHERE s.experiment_id=?
                ORDER BY s.sort_order, s.sample""", (gene["gene_id"], experiment_id))]
            item["curated_annotations"] = [dict(row) for row in con.execute(
                "SELECT source, field, value FROM gene_annotations WHERE gene_id=? ORDER BY source, field",
                (gene["gene_id"],)
            )]
        else:
            item["expression"] = [dict(row) for row in con.execute("""SELECT s.sample, s.sample AS label,
                s.timepoint, s.replicate, '' AS group_name, '' AS treatment, '' AS color, '' AS comments,
                COALESCE(e.tpm,0) AS tpm, COALESCE(e.reads,0) AS reads
                FROM samples s LEFT JOIN expression e ON e.sample=s.sample AND e.gene_id=?
                ORDER BY s.sort_order, s.sample""", (gene["gene_id"],))]
            for row in item["expression"]:
                row["label"] = " ".join(x for x in [row.get("timepoint"), row.get("replicate")] if x) or row["sample"]
            item["curated_annotations"] = []
        return item

    def do_GET(self):
        parsed = urlparse(self.path)
        if not parsed.path.startswith("/api/"):
            return super().do_GET()
        con = self.connect()
        try:
            qs = parse_qs(parsed.query)
            experiment_id = qs.get("experiment", [self.default_experiment(con)])[0]
            if parsed.path == "/api/info":
                if self.is_multi_experiment(con):
                    experiments = [dict(row) for row in con.execute("SELECT * FROM experiments ORDER BY name")]
                    samples = [dict(row) for row in con.execute(
                        "SELECT * FROM samples WHERE experiment_id=? ORDER BY sort_order, sample", (experiment_id,)
                    )]
                else:
                    experiments = [{"experiment_id": "default", "name": "Time course salinidad"}]
                    samples = [dict(row) for row in con.execute("SELECT *, sample AS label FROM samples ORDER BY sort_order, sample")]
                    for row in samples:
                        row["label"] = " ".join(x for x in [row.get("timepoint"), row.get("replicate")] if x) or row["sample"]
                genes = con.execute("SELECT COUNT(*) FROM genes").fetchone()[0]
                return self.send_json({"samples": samples, "gene_count": genes, "experiments": experiments, "active_experiment": experiment_id})

            if parsed.path == "/api/search":
                query = qs.get("q", [""])[0].strip().lower()
                if len(query) < 2:
                    return self.send_json([])
                escaped = query.replace("%", "\\%").replace("_", "\\_")
                rows = con.execute("""SELECT gene_id, alias, name, description, seqid, start, end, strand
                    FROM genes WHERE search_text LIKE ? ESCAPE '\\' LIMIT 300""", (f"%{escaped}%",)).fetchall()
                ranked = sorted(rows, key=lambda row: (-self.rank_gene(row, query), row["gene_id"]))[:10]
                return self.send_json([dict(row) | {"match_score": round(self.rank_gene(row, query), 2)} for row in ranked])

            if parsed.path == "/api/batch":
                raw = qs.get("q", [""])[0]
                terms = [term.strip().lower() for term in re.split(r"[,;\n\r\t]+", raw) if term.strip()]
                if not terms:
                    return self.send_json([])
                if len(terms) > 50:
                    return self.send_json({"error": "Máximo 50 genes por consulta"}, 400)
                found, seen = [], set()
                for term in terms:
                    escaped = term.replace("%", "\\%").replace("_", "\\_")
                    rows = con.execute("""SELECT * FROM genes WHERE search_text LIKE ? ESCAPE '\\' LIMIT 100""",
                                       (f"%{escaped}%",)).fetchall()
                    if not rows:
                        continue
                    gene = sorted(rows, key=lambda row: (-self.rank_gene(row, term), row["gene_id"]))[0]
                    if gene["gene_id"] in seen:
                        continue
                    seen.add(gene["gene_id"])
                    found.append(self.gene_payload(con, gene, experiment_id))
                return self.send_json(found)

            if parsed.path.startswith("/api/gene/"):
                gene_id = parsed.path[len("/api/gene/"):]
                gene = con.execute("SELECT * FROM genes WHERE gene_id=?", (gene_id,)).fetchone()
                if not gene:
                    return self.send_json({"error": "Gen no encontrado"}, 404)
                return self.send_json(self.gene_payload(con, gene, experiment_id))
            return self.send_json({"error": "Ruta no encontrada"}, 404)
        finally:
            con.close()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--quant-dir", type=Path, default=Path(os.environ["QUANT_DIR"]) if os.environ.get("QUANT_DIR") else None)
    parser.add_argument("--gff", type=Path, default=Path(os.environ["GFF_PATH"]) if os.environ.get("GFF_PATH") else None)
    parser.add_argument("--db", type=Path, default=DEFAULT_DB)
    parser.add_argument("--reindex", action="store_true")
    parser.add_argument("--index-only", action="store_true")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=int(os.environ.get("PORT", "8765")))
    args = parser.parse_args()
    if args.reindex or not args.db.exists():
        build_index(args.config, args.quant_dir, args.gff, args.db)
    if args.index_only:
        return
    Handler.db_path = args.db
    server = ThreadingHTTPServer((args.host, args.port), Handler)
    print(f"Dashboard: http://{args.host}:{args.port}")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nServidor detenido.")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        raise
