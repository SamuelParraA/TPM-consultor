#!/usr/bin/env python3
import argparse
import json
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--info", type=Path, required=True)
parser.add_argument("--version", required=True)
parser.add_argument("--k", type=int, required=True)
args = parser.parse_args()
if not args.info.is_file():
    raise SystemExit(f"ERROR: falta {args.info}")
info = json.loads(args.info.read_text(encoding="utf-8"))
errors = []
expected_seq_hash = "4952695c9309330ef494ef455db7371ffce898b39fd3530fd9b828954af29985"
expected_name_hash = "86755ceaf8df52514513f9218406de85f55979d49808bdf1b395230856f4adcb"
if str(info.get("salmon_version")) != args.version:
    errors.append(f"versión {info.get('salmon_version')} != {args.version}")
if int(info.get("k", -1)) != args.k:
    errors.append(f"k={info.get('k')} != {args.k}")
if int(info.get("num_refs", -1)) != 33976:
    errors.append(f"se esperaban 33.976 transcritos y hay {info.get('num_refs')}")
if info.get("seq_hash") != expected_seq_hash:
    errors.append("hash interno de secuencias distinto al índice histórico")
if info.get("name_hash") != expected_name_hash:
    errors.append("hash interno de nombres distinto al índice histórico")
if errors:
    raise SystemExit("ERROR: índice incompatible: " + "; ".join(errors))
print("OK: índice Salmon validado (versión, k y número de transcritos).")
