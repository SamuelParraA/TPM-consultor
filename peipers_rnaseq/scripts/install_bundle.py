#!/usr/bin/env python3
"""Instala un paquete validado en una raíz de datos del mantenedor."""

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--bundle", type=Path, required=True)
parser.add_argument("--destination", type=Path, required=True)
parser.add_argument("--experiment-id", default="ethylene_rootstocks_salt")
parser.add_argument("--experiment-name", default="Etileno del portainjerto y salinidad (PRJNA1347747)")
args = parser.parse_args()
subprocess.run([sys.executable, str(Path(__file__).with_name("verify_bundle.py")), "--bundle", str(args.bundle)], check=True)
if args.destination.exists():
    raise SystemExit(f"ERROR: {args.destination} ya existe; no se sobrescribirá")
temporary = args.destination.with_name(args.destination.name + ".building")
if temporary.exists():
    raise SystemExit(f"ERROR: existe instalación incompleta: {temporary}")
temporary.mkdir(parents=True)
shutil.copytree(args.bundle / "quant", temporary / "quant")
shutil.copy2(args.bundle / "metadata" / "sample_metadata.tsv", temporary / "quant" / "sample_metadata.tsv")
shutil.copy2(args.bundle / "metadata" / "run_manifest.tsv", temporary / "run_manifest.tsv")
shutil.copy2(args.bundle / "dataset_manifest.json", temporary / "dataset_manifest.json")
args.destination.parent.mkdir(parents=True, exist_ok=True)
temporary.rename(args.destination)
snippet = {
    "id": args.experiment_id,
    "name": args.experiment_name,
    "path": str(args.destination / "quant"),
    "metadata": str(args.destination / "quant" / "sample_metadata.tsv"),
}
print(f"OK: dataset instalado en {args.destination}")
print("Añada este objeto a experiments en su config.local.json:")
print(json.dumps(snippet, indent=2, ensure_ascii=False))
