#!/usr/bin/env python3
"""Dashboard local para explorar cuantificaciones Salmon cruzadas con un GFF3."""

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
DEFAULT_QUANT = Path(os.environ.get("QUANT_DIR", r"I:\transcriptomica\salmon_timecourse"))
DEFAULT_GFF = Path(os.environ.get("GFF_PATH", r"I:\transcriptomica\reference\ITAG4.0_gene_models.gff"))
DEFAULT_METADATA = Path(os.environ.get("METADATA_PATH", str(ROOT / "sample_metadata.csv")))
DEFAULT_DB = Path(os.environ.get("DB_PATH", str(ROOT / "rnaseq_index.sqlite")))


def strip_prefix(value: str) -> str:
    return re.sub(r"^(gene|mRNA|transcript):", "", value or "")


def attributes(text: str) -> dict[str, str]:
    result = {}
    for item in text.rstrip().split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = value.replace("%20", " ").replace("%2C", ",")
    return result


def load_metadata(path: Path, sample_names: list[str]) -> list[dict[str, str]]:
    configured = {}
    if path.exists():
        with path.open(encoding="utf-8-sig", newline="") as handle:
            configured = {row["sample"].strip(): row for row in csv.DictReader(handle)}
    rows = []
    for order, name in enumerate(sample_names, 1):
        row = configured.get(name, {})
        rows.append({
            "sample": name,
            "timepoint": row.get("timepoint", name).strip() or name,
            "replicate": row.get("replicate", "").strip(),
            "order": int(row.get("order", order) or order),
        })
    return sorted(rows, key=lambda row: (row["order"], row["sample"]))


def build_index(quant_dir: Path, gff_path: Path, metadata_path: Path, db_path: Path) -> None:
    quant_files = sorted(quant_dir.glob("*/quant.sf"))
    if not quant_files:
        raise FileNotFoundError(f"No se encontraron archivos */quant.sf en {quant_dir}")
    if not gff_path.exists():
        raise FileNotFoundError(f"No se encontró el GFF: {gff_path}")

    print(f"Leyendo anotación: {gff_path}")
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
            if feature == "gene":
                gene_id = strip_prefix(attrs.get("ID", attrs.get("Name", "")))
                if gene_id:
                    genes[gene_id] = {
                        "gene_id": gene_id,
                        "alias": attrs.get("Alias", ""),
                        "name": attrs.get("Name", gene_id),
                        "description": attrs.get("Note", attrs.get("description", "")),
                        "seqid": seqid, "start": int(start), "end": int(end), "strand": strand,
                        "transcripts": [],
                    }
            elif feature in {"mRNA", "transcript"}:
                transcript_id = strip_prefix(attrs.get("ID", attrs.get("Name", "")))
                gene_id = strip_prefix(attrs.get("Parent", "").split(",")[0])
                if transcript_id and gene_id:
                    transcript_to_gene[transcript_id] = gene_id
                    if gene_id in genes:
                        genes[gene_id]["transcripts"].append(transcript_id)
                        if attrs.get("Note") and not genes[gene_id]["description"]:
                            genes[gene_id]["description"] = attrs["Note"]

    metadata = load_metadata(metadata_path, [path.parent.name for path in quant_files])
    meta_by_sample = {row["sample"]: row for row in metadata}
    expression: dict[tuple[str, str], list[float]] = defaultdict(lambda: [0.0, 0.0])
    for index, quant_path in enumerate(quant_files, 1):
        sample = quant_path.parent.name
        print(f"[{index}/{len(quant_files)}] {sample}")
        with quant_path.open(encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                transcript = row["Name"]
                gene_id = transcript_to_gene.get(transcript)
                if not gene_id:
                    # ITAG: Solyc01g000010.4.1 -> Solyc01g000010.4
                    candidate = transcript.rsplit(".", 1)[0]
                    gene_id = candidate if candidate in genes else None
                if gene_id:
                    values = expression[(gene_id, sample)]
                    values[0] += float(row["TPM"])
                    values[1] += float(row["NumReads"])

    db_path.unlink(missing_ok=True)
    con = sqlite3.connect(db_path)
    con.executescript("""
        CREATE TABLE genes (gene_id TEXT PRIMARY KEY, alias TEXT, name TEXT, description TEXT,
          seqid TEXT, start INTEGER, end INTEGER, strand TEXT, transcripts TEXT, search_text TEXT);
        CREATE TABLE samples (sample TEXT PRIMARY KEY, timepoint TEXT, replicate TEXT, sort_order INTEGER);
        CREATE TABLE expression (gene_id TEXT, sample TEXT, tpm REAL, reads REAL,
          PRIMARY KEY (gene_id, sample));
        CREATE INDEX idx_gene_search ON genes(search_text);
        CREATE INDEX idx_expr_gene ON expression(gene_id);
    """)
    con.executemany("INSERT INTO samples VALUES (?,?,?,?)", [
        (m["sample"], m["timepoint"], m["replicate"], m["order"]) for m in metadata
    ])
    gene_rows = []
    for gene in genes.values():
        tx = ",".join(gene["transcripts"])
        search = " ".join([gene["gene_id"], gene["alias"], gene["name"], gene["description"], tx]).lower()
        gene_rows.append((gene["gene_id"], gene["alias"], gene["name"], gene["description"],
                          gene["seqid"], gene["start"], gene["end"], gene["strand"], tx, search))
    con.executemany("INSERT INTO genes VALUES (?,?,?,?,?,?,?,?,?,?)", gene_rows)
    con.executemany("INSERT INTO expression VALUES (?,?,?,?)", [
        (gene_id, sample, values[0], values[1]) for (gene_id, sample), values in expression.items()
    ])
    con.commit()
    con.close()
    print(f"Índice listo: {len(genes):,} genes, {len(metadata)} muestras -> {db_path}")


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

    def do_GET(self):
        parsed = urlparse(self.path)
        if not parsed.path.startswith("/api/"):
            return super().do_GET()
        con = sqlite3.connect(self.db_path)
        con.row_factory = sqlite3.Row
        try:
            if parsed.path == "/api/info":
                samples = [dict(row) for row in con.execute("SELECT * FROM samples ORDER BY sort_order, sample")]
                genes = con.execute("SELECT COUNT(*) FROM genes").fetchone()[0]
                return self.send_json({"samples": samples, "gene_count": genes})
            if parsed.path == "/api/search":
                query = parse_qs(parsed.query).get("q", [""])[0].strip().lower()
                if len(query) < 2:
                    return self.send_json([])
                escaped = query.replace("%", "\\%").replace("_", "\\_")
                rows = con.execute("""SELECT gene_id, alias, name, description, seqid, start, end, strand
                    FROM genes WHERE search_text LIKE ? ESCAPE '\\'
                    LIMIT 250""", (f"%{escaped}%",)).fetchall()

                def relevance(row):
                    gene_id = row["gene_id"].lower()
                    alias = (row["alias"] or "").lower()
                    name = (row["name"] or "").lower()
                    description = (row["description"] or "").lower()
                    identifiers = (gene_id, alias, name)
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
                    similarity = max((SequenceMatcher(None, query, value).ratio()
                                      for value in windows), default=0)
                    return phrase_score + similarity * 100

                ranked = sorted(rows, key=lambda row: (-relevance(row), row["gene_id"]))[:10]
                return self.send_json([dict(row) | {"match_score": round(relevance(row), 2)}
                                       for row in ranked])
            if parsed.path == "/api/batch":
                raw = parse_qs(parsed.query).get("q", [""])[0]
                # Los espacios forman parte de descripciones; solo coma, punto y coma o salto separan genes.
                terms = [term.strip().lower() for term in re.split(r"[,;\n\r\t]+", raw) if term.strip()]
                if not terms:
                    return self.send_json([])
                if len(terms) > 50:
                    return self.send_json({"error": "Máximo 50 genes por consulta"}, 400)
                found, seen = [], set()
                for term in terms:
                    escaped = term.replace("%", "\\%").replace("_", "\\_")
                    gene = con.execute("""SELECT * FROM genes WHERE search_text LIKE ? ESCAPE '\\'
                        ORDER BY CASE WHEN lower(gene_id)=? THEN 0 WHEN lower(alias)=? THEN 1
                        WHEN lower(name)=? THEN 2 ELSE 3 END, gene_id LIMIT 1""",
                        (f"%{escaped}%", term, term, term)).fetchone()
                    if not gene or gene["gene_id"] in seen:
                        continue
                    seen.add(gene["gene_id"])
                    item = dict(gene)
                    item.pop("search_text", None)
                    item["transcripts"] = item["transcripts"].split(",") if item["transcripts"] else []
                    item["expression"] = [dict(row) for row in con.execute("""SELECT s.sample,
                        s.timepoint, s.replicate, COALESCE(e.tpm,0) AS tpm, COALESCE(e.reads,0) AS reads
                        FROM samples s LEFT JOIN expression e ON e.sample=s.sample AND e.gene_id=?
                        ORDER BY s.sort_order, s.sample""", (gene["gene_id"],))]
                    found.append(item)
                return self.send_json(found)
            if parsed.path.startswith("/api/gene/"):
                gene_id = parsed.path[len("/api/gene/"):]
                gene = con.execute("SELECT * FROM genes WHERE gene_id=?", (gene_id,)).fetchone()
                if not gene:
                    return self.send_json({"error": "Gen no encontrado"}, 404)
                values = [dict(row) for row in con.execute("""SELECT s.sample, s.timepoint, s.replicate,
                    COALESCE(e.tpm,0) AS tpm, COALESCE(e.reads,0) AS reads
                    FROM samples s LEFT JOIN expression e ON e.sample=s.sample AND e.gene_id=?
                    ORDER BY s.sort_order, s.sample""", (gene_id,))]
                result = dict(gene)
                result["transcripts"] = result["transcripts"].split(",") if result["transcripts"] else []
                result.pop("search_text", None)
                result["expression"] = values
                return self.send_json(result)
            return self.send_json({"error": "Ruta no encontrada"}, 404)
        finally:
            con.close()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quant-dir", type=Path, default=DEFAULT_QUANT)
    parser.add_argument("--gff", type=Path, default=DEFAULT_GFF)
    parser.add_argument("--metadata", type=Path, default=DEFAULT_METADATA)
    parser.add_argument("--db", type=Path, default=DEFAULT_DB)
    parser.add_argument("--reindex", action="store_true")
    parser.add_argument("--index-only", action="store_true")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8765)
    args = parser.parse_args()
    if args.reindex or not args.db.exists():
        build_index(args.quant_dir, args.gff, args.metadata, args.db)
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
