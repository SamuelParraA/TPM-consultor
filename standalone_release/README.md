# Perfil de expresión en 75 mM de salinidad

Aplicación web portable para explorar TPM de una serie temporal de RNA-seq de tomate. Cruza archivos `quant.sf` de Salmon con `ITAG4.0_gene_models.gff`, agrega isoformas a nivel de gen y resume tres réplicas biológicas en `0h`, `6h`, `12h`, `24h`, `48h` y `96h`.

## Funciones

- Búsqueda parcial por ID ITAG4, `Name`, `Alias`, descripción `Note` o transcrito.
- Consulta simultánea de hasta 50 términos separados por comas.
- TPM promedio, desviación estándar y valores de réplicas, sin mostrar accesiones SRR.
- Gráfico individual con puntos, promedio y barras de error.
- Heatmap, perfiles múltiples y orden por máximo, fold-change o similitud temporal.
- Exportación CSV, Excel, SVG, PNG y PDF mediante impresión del navegador.
- Interfaz íntegramente en Verdana.

## Ejecución en Windows

Requiere Python 3.10 o posterior y no instala dependencias:

```powershell
python app.py --reindex
```

Después abre <http://127.0.0.1:8765>. En ejecuciones posteriores basta `python app.py`.

Rutas predeterminadas:

- Salmon: `I:\transcriptomica\salmon_timecourse\*/quant.sf`
- GFF: `I:\transcriptomica\reference\ITAG4.0_gene_models.gff`

También se aceptan `--quant-dir`, `--gff`, `--metadata`, `--db`, o las variables `QUANT_DIR`, `GFF_PATH`, `METADATA_PATH` y `DB_PATH`.

## Metadata

`sample_metadata.csv` es la configuración editable. Columnas: `sample`, `timepoint`, `replicate`, `order`. Tras cambiarla hay que ejecutar `python app.py --reindex`. Si faltan temporalidades en el directorio de Salmon, se presentan como datos ausentes y no como TPM cero.

## Docker Compose

Crea un archivo `.env`:

```text
SALMON_DIR=/ruta/a/salmon_timecourse
REFERENCE_DIR=/ruta/a/reference
```

Luego:

```bash
docker compose up --build
```

El índice persistirá en `./data`. Compatible con Linux, Windows/WSL, DockerHub y GitHub Container Registry.

## Docker directo

```bash
docker build -t expresion-tomate .
docker run --rm -p 8765:8765 \
  -v /ruta/salmon:/data/salmon:ro \
  -v /ruta/reference:/data/reference:ro \
  -v "$(pwd)/data:/app/data" expresion-tomate
```

## Pruebas

```bash
python -m unittest discover -s tests
```

## Nota sobre nombres comunes

La búsqueda solo puede resolver símbolos como `ACS7` si aparecen en `Name`, `Alias`, `Note` u otro texto del GFF proporcionado. ITAG4.0 suele usar identificadores `Solyc...` como nombre; no se inventan correspondencias ausentes.

Licencia MIT.

## Versión autónoma para GitHub

El archivo `rnaseq_index.sqlite` incluido contiene los datos necesarios para buscar y visualizar genes sin acceder al GFF ni a los directorios Salmon originales. Consulta [STANDALONE.md](STANDALONE.md) para publicación y despliegue con `Dockerfile.standalone`.
