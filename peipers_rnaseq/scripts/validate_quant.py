#!/usr/bin/env python3
import argparse
import csv
import json
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--quant-dir", type=Path, required=True)
parser.add_argument("--sample", required=True)
args = parser.parse_args()
quant = args.quant_dir / "quant.sf"
meta = args.quant_dir / "aux_info" / "meta_info.json"
if not quant.is_file() or not meta.is_file():
    raise SystemExit(f"ERROR: cuantificación incompleta de {args.sample}")

rows = 0
tpm_sum = 0.0
with quant.open(encoding="utf-8", newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    if reader.fieldnames != ["Name", "Length", "EffectiveLength", "TPM", "NumReads"]:
        raise SystemExit(f"ERROR: columnas inesperadas en {quant}")
    for row in reader:
        rows += 1
        tpm_sum += float(row["TPM"])
information = json.loads(meta.read_text(encoding="utf-8"))
rate = float(information.get("percent_mapped", -1))
version = str(information.get("salmon_version", ""))
errors = []
if rows != 33976:
    errors.append(f"{rows} transcritos (esperados: 33976)")
if not 999_000 <= tpm_sum <= 1_001_000:
    errors.append(f"suma TPM={tpm_sum:.2f}, fuera del rango esperado")
if version != "2.0.1":
    errors.append(f"Salmon {version}, esperado 2.0.1")
if rate < 40:
    errors.append(f"tasa de mapeo anormalmente baja: {rate:.2f}%")
if errors:
    raise SystemExit(f"ERROR en {args.sample}: " + "; ".join(errors))
print(f"OK {args.sample}: {rows} transcritos; TPM={tpm_sum:.1f}; mapeo={rate:.2f}%")
