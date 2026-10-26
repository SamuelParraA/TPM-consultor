#!/usr/bin/env python3
"""Construye una copia autocontenida de TPM-consultor para un experimento."""

import argparse
import csv
import json
import shutil
import subprocess
import sys
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, required=True)
    parser.add_argument("--quant-dir", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--gff", type=Path, required=True)
    parser.add_argument("--experiment-id", required=True)
    parser.add_argument("--experiment-name", required=True)
    parser.add_argument("--sample-count", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    app = args.repo_root / "app.py"
    web = args.repo_root / "web"
    if not app.is_file() or not web.is_dir():
        raise SystemExit("ERROR: no se encontró app.py/web en la raíz del fork TPM-consultor")
    if not args.gff.is_file():
        raise SystemExit(f"ERROR: falta {args.gff}")
    if len(list(args.quant_dir.glob("SRR*/quant.sf"))) != args.sample_count:
        raise SystemExit(f"ERROR: se requieren {args.sample_count} quant.sf antes de construir el consultor")
    if args.output.exists():
        raise SystemExit(f"ERROR: {args.output} ya existe; no se sobrescribirá")

    temporary = args.output.with_name(args.output.name + ".building")
    if temporary.exists():
        raise SystemExit(f"ERROR: quedó una construcción incompleta en {temporary}; revísela y muévala")
    temporary.mkdir(parents=True)
    database = temporary / "rnaseq_index.sqlite"
    config = {
        "data_root": str(args.quant_dir.parent),
        "reference_dir": str(args.gff.parent),
        "annotation_dir": str(args.gff.parent / "annotation_optional"),
        "gff": str(args.gff),
        "experiments": [{
            "id": args.experiment_id,
            "name": args.experiment_name,
            "path": str(args.quant_dir),
            "metadata": str(args.metadata),
        }],
    }
    config_path = temporary / "build_config.json"
    config_path.write_text(json.dumps(config, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    command = [
        sys.executable, str(app), "--config", str(config_path), "--db", str(database),
        "--reindex", "--index-only",
    ]
    print("Construyendo SQLite del consultor...")
    subprocess.run(command, check=True)

    shutil.copy2(app, temporary / "app.py")
    shutil.copytree(web, temporary / "web")
    for name in (
        "requirements.txt", "Dockerfile.standalone", "docker-compose.standalone.yml",
        "iniciar_linux.sh", "iniciar_windows.bat", "LICENSE",
    ):
        source = args.repo_root / name
        if source.is_file():
            shutil.copy2(source, temporary / name)

    with args.metadata.open(encoding="utf-8-sig", newline="") as source, \
            (temporary / "sample_metadata.csv").open("w", encoding="utf-8", newline="") as target:
        reader = csv.DictReader(source, delimiter="\t")
        writer = csv.DictWriter(target, fieldnames=reader.fieldnames or [], lineterminator="\n")
        writer.writeheader()
        writer.writerows(reader)
    shutil.copy2(args.metadata, temporary / "sample_metadata.tsv")
    readme = f"""# Consultor TPM - {args.experiment_name}

Consultor standalone generado a partir de {args.sample_count} bibliotecas de {args.experiment_name},
cuantificadas con Salmon 2.0.1 contra ITAG4.0.

## Abrir localmente

```bash
python app.py
```

Abra http://127.0.0.1:8765. El SQLite ya contiene la expresión.
"""
    (temporary / "README.md").write_text(readme, encoding="utf-8")
    config_path.unlink()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary.rename(args.output)
    print(f"OK: consultor creado en {args.output}")


if __name__ == "__main__":
    main()
