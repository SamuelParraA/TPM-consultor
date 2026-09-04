# Perfil de expresión transcriptómica en tomate

Dashboard web en Python para explorar resultados RNA-seq cuantificados con Salmon y anotados con ITAG4.0. La idea ya no es solo ver TPM: el proyecto queda preparado como una plataforma extensible de exploración transcriptómica.

## Pipeline docente Peipers-RNAseq

La carpeta [`peipers_rnaseq/`](peipers_rnaseq/) contiene el pipeline reproducible
para que un alumno procese RNA-seq de tomate desde ENA/SRA hasta matrices TPM,
un paquete integrable en este consultor y, opcionalmente, DESeq2.

Para comenzar, consulte la
[`guía paso a paso`](peipers_rnaseq/README.md) o su
[`versión PDF`](documentacion/Guia_Peipers_RNAseq_Paso_a_Paso.pdf).

La guía y el pipeline Peipers-RNAseq fueron creados por **Samuel Parra A.,
PhD**, investigador postdoctoral en el laboratorio de la Dra. Claudia Stange,
Centro de Biotecnología Molecular Vegetal, Departamento de Ciencias,
Universidad de Chile. Contacto: [sa.parra@uandresbello.edu](mailto:sa.parra@uandresbello.edu).
ORCID: [0000-0002-9129-4133](https://orcid.org/0000-0002-9129-4133).

Copyright (c) 2026 Samuel Parra A., PhD. Licencia MIT.

## Qué hace

- Detecta experimentos con subcarpetas que contienen `quant.sf`.
- Usa `ITAG4.0_gene_models.gff` como fuente oficial de anotación.
- Busca genes por ITAG, alias, nombre, descripción, GO, familia o pathway.
- Acepta varios genes en una misma consulta, separados por coma, punto y coma o salto de línea.
- Muestra sugerencias desplegables con las 10 coincidencias más similares.
- Usa metadatos para mostrar nombres biológicos de muestras en vez de SRR.
- Resume TPM por condición: promedio, réplicas, desviación estándar y número de muestras.
- Genera gráfico individual, perfiles agrupados, leyenda visible y heatmap.
- Exporta CSV, Excel, SVG, PNG y PDF.
- Usa Verdana y tema claro.

## Estructura esperada

```text
I:\transcriptomica

reference/
    ITAG4.0_gene_models.gff
    ITAG4.0_cDNA.fasta
    salmon_index/

annotation/
    alias.tsv
    gene_family.tsv
    pathways.tsv
    genes_interest.tsv
    go_annotations.tsv
    literature/

salmon_timecourse/
salmon_grafted/
salmon_test/
analysis/
```

Cualquier carpeta dentro de la raíz que tenga subcarpetas con `quant.sf` se considera un experimento.

## Configuración

Edita `config.example.json` o crea tu propio archivo:

```json
{
  "data_root": "I:\\transcriptomica",
  "reference_dir": "I:\\transcriptomica\\reference",
  "annotation_dir": "I:\\transcriptomica\\annotation",
  "gff": "I:\\transcriptomica\\reference\\ITAG4.0_gene_models.gff",
  "experiments": []
}
```

Si `experiments` queda vacío, la app detecta experimentos automáticamente.

## Metadatos por experimento

Cada experimento puede tener un `sample_metadata.tsv` o `sample_metadata.csv` con columnas como:

```text
SampleID    Nombre biológico    Grupo    Tratamiento    Tiempo    Réplica    Color sugerido    Comentarios
```

Si no existe, la app usa el nombre de la carpeta de la muestra como respaldo.

## Anotaciones opcionales

Los TSV/CSV dentro de `annotation/` enriquecen la información, pero no reemplazan el GFF. Deben incluir una columna `ITAG`, `gene_id`, `GeneID` o `Gene ID`.

Ejemplos:

- `alias.tsv`: `ITAG`, `Alias`
- `gene_family.tsv`: `ITAG`, `Family`
- `pathways.tsv`: `ITAG`, `Pathway`
- `genes_interest.tsv`: `ITAG`, `Category`
- `go_annotations.tsv`: `ITAG`, `GO`, `Description`

## Ejecutar localmente

```powershell
python app.py --reindex
python app.py
```

Abre:

<http://127.0.0.1:8765>

También puedes usar una configuración explícita:

```powershell
python app.py --config config.example.json --reindex
```

## Versión standalone

La carpeta `standalone_release/` incluye `rnaseq_index.sqlite`, por lo que no necesita consultar `I:\transcriptomica`. Sirve para publicar o probar el buscador ya indexado.

```powershell
cd standalone_release
python app.py
```

## Docker

```bash
docker compose up --build
```

La imagen espera una raíz `/data` con `reference/`, `annotation/` y carpetas de experimentos.

## Publicación

GitHub Pages solo publica sitios estáticos; esta app necesita Python + SQLite. Para tener un enlace público, úsala en Render, Railway, Fly.io o un VPS con Docker. El `Dockerfile.standalone` ya usa la variable `PORT`, compatible con Render.

## Pruebas

```bash
python -m unittest discover -s tests
```

Licencia MIT.
