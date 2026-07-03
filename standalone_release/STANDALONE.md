# Publicación autónoma

Este repositorio funciona sin acceso al disco original porque `rnaseq_index.sqlite` contiene las anotaciones ITAG4.0 y los TPM ya procesados. Los archivos `quant.sf` y el GFF no son necesarios durante la consulta.

## GitHub

Publica como mínimo:

- `app.py`
- `rnaseq_index.sqlite`
- `sample_metadata.csv`
- carpeta `web/`
- `Dockerfile.standalone`
- `docker-compose.standalone.yml`
- `README.md`, `LICENSE` y `.github/`

El índice pesa aproximadamente 62 MB y cabe bajo el límite de 100 MB por archivo de GitHub. Para históricos grandes conviene usar Git LFS o adjuntarlo a una GitHub Release.

## Ejecución directa

Windows: doble clic en `iniciar_windows.bat`.

Linux/WSL:

```bash
chmod +x iniciar_linux.sh
./iniciar_linux.sh
```

## Docker autónomo

```bash
docker compose -f docker-compose.standalone.yml up --build
```

Abre <http://localhost:8765>. El contenedor no monta ni consulta directorios del equipo anfitrión.

## Actualización futura de datos

La reconstrucción sí requiere los `quant.sf` y el GFF:

```bash
python app.py --reindex --quant-dir /ruta/salmon --gff /ruta/ITAG4.0_gene_models.gff
```

Después reemplaza `rnaseq_index.sqlite` en el repositorio.
