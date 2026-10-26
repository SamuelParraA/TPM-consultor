#!/usr/bin/env python3
import argparse
import sqlite3
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--directory", type=Path, required=True)
parser.add_argument("--experiment", required=True)
parser.add_argument("--samples", type=int, required=True)
args = parser.parse_args()
database = args.directory / "rnaseq_index.sqlite"
if not database.is_file():
    raise SystemExit(f"ERROR: falta {database}")
with sqlite3.connect(database) as connection:
    tables = {row[0] for row in connection.execute("SELECT name FROM sqlite_master WHERE type='table'")}
    required = {"experiments", "genes", "samples", "expression"}
    if not required <= tables:
        raise SystemExit(f"ERROR: esquema SQLite incompleto: {tables}")
    experiment = connection.execute("SELECT name FROM experiments WHERE experiment_id = ?", (args.experiment,)).fetchone()
    sample_count = connection.execute("SELECT COUNT(*) FROM samples WHERE experiment_id = ?", (args.experiment,)).fetchone()[0]
    gene_count = connection.execute("SELECT COUNT(*) FROM genes").fetchone()[0]
    expression_count = connection.execute("SELECT COUNT(*) FROM expression WHERE experiment_id = ?", (args.experiment,)).fetchone()[0]
errors = []
if experiment is None:
    errors.append(f"falta experimento {args.experiment}")
if sample_count != args.samples:
    errors.append(f"hay {sample_count} muestras; se esperaban {args.samples}")
if not 33_000 <= gene_count <= 35_000:
    errors.append(f"número de genes inesperado: {gene_count}")
if expression_count < 33_000 * args.samples:
    errors.append(f"muy pocas filas de expresión: {expression_count}")
if errors:
    raise SystemExit("ERROR: " + "; ".join(errors))
size = database.stat().st_size / 1024**2
print(f"OK: consultor; {sample_count} muestras; {gene_count} genes; {expression_count} TPM; SQLite {size:.1f} MiB.")
