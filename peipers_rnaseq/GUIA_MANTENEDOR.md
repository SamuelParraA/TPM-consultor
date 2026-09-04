# Guía del mantenedor: integrar un dataset al consultor

El alumno entrega un pull request con un paquete en **datasets/PRJNA...**;
nunca entrega FASTQ.

## 1. Revisar el pull request

~~~bash
python peipers_rnaseq/scripts/verify_bundle.py --bundle datasets/PRJNA1347747
python -m unittest discover -s tests
~~~

El verificador comprueba número de quant.sf, SHA-256 y tamaños, presencia de
metadata/matrices/tx2gene, ausencia de datos crudos y la declaración
raw_reads_included=false.

Revise además:

- BioProject y artículo;
- condición y réplica de cada SRR;
- Salmon 2.0.1 en meta_info.json;
- tasas de mapeo en salmon_qc.tsv;
- MultiQC;
- SHA-256 de la referencia.

## 2. Instalar en la raíz de datos

En Windows, desde la raíz del repositorio:

~~~powershell
python peipers_rnaseq/scripts/install_bundle.py --bundle datasets/PRJNA1347747 --destination "I:\transcriptomica\PRJNA1347747"
~~~

El instalador valida primero y se niega a sobrescribir un destino existente.

~~~text
I:\transcriptomica\PRJNA1347747\
    dataset_manifest.json
    run_manifest.tsv
    quant\
        sample_metadata.tsv
        SRR35859425\quant.sf
        ...
~~~

## 3. Configurar el consultor unificado

No edite config.example.json con rutas privadas. Cree config.local.json,
ignorado por Git, con todos los experimentos:

~~~json
{
  "data_root": "I:\\transcriptomica",
  "reference_dir": "I:\\transcriptomica\\reference",
  "annotation_dir": "I:\\transcriptomica\\annotation",
  "gff": "I:\\transcriptomica\\reference\\ITAG4.0_gene_models.gff",
  "experiments": [
    {
      "id": "salt_timecourse_selfrooted",
      "name": "Curso temporal de salinidad - raíz propia",
      "path": "I:\\transcriptomica\\salmon_timecourse",
      "metadata": "I:\\transcriptomica\\salmon_timecourse\\sample_metadata.tsv"
    },
    {
      "id": "salt_timecourse_grafted",
      "name": "Curso temporal de salinidad - injertado",
      "path": "I:\\transcriptomica\\salmon_grafted",
      "metadata": "I:\\transcriptomica\\salmon_grafted\\sample_metadata.tsv"
    },
    {
      "id": "ethylene_rootstocks_salt",
      "name": "Etileno del portainjerto y salinidad (PRJNA1347747)",
      "path": "I:\\transcriptomica\\PRJNA1347747\\quant",
      "metadata": "I:\\transcriptomica\\PRJNA1347747\\quant\\sample_metadata.tsv"
    }
  ]
}
~~~

## 4. Reconstruir SQLite

~~~powershell
python app.py --config config.local.json --db rnaseq_index.sqlite --reindex --index-only
python app.py --config config.local.json --db rnaseq_index.sqlite
~~~

Abra [http://127.0.0.1:8765](http://127.0.0.1:8765) y confirme:

- selector con todos los experimentos;
- 18 muestras del proyecto nuevo;
- seis condiciones y tres réplicas;
- búsquedas de Solyc07g049530 y Solyc12g005940;
- exportaciones y gráficos.

## 5. Publicar el SQLite

El SQLite unificado puede superar 100 MiB. No lo divida ni lo añada
repetidamente al historial Git.

Opciones:

1. **GitHub Release** como artefacto versionado.
2. **Git LFS** sólo si se acepta su cuota.
3. Descargar un artefacto de Release durante el despliegue.

Registre fecha, commit, BioProjects, SHA-256 y versión del esquema:

~~~powershell
Get-FileHash -Algorithm SHA256 rnaseq_index.sqlite
~~~

## 6. Consultor standalone

El alumno genera WORK_DIR/consultor_standalone/PRJNA1347747. Puede publicarse
separadamente; es preferible adjuntar SQLite a un Release.

~~~bash
python peipers_rnaseq/scripts/verify_consultor.py --directory RUTA_CONSULTOR --experiment ethylene_rootstocks_salt --samples 18
~~~

## 7. Política para futuros pull requests

~~~text
datasets/PRJNA.../
    README.md
    SHA256SUMS
    dataset_manifest.json
    metadata/
    quant/
    matrices/
    qc/
~~~

Requisitos:

- tomate y referencia ITAG4.0;
- Salmon 2.0.1 y parámetros de compatibilidad;
- metadata revisada manualmente;
- checksums;
- al menos dos réplicas por condición si habrá DE;
- ningún FASTQ, SRA, índice o SQLite;
- verificador y pruebas en verde.

## 8. Interpretación entre estudios

El consultor permite comparar visualmente TPM procesados con la misma
referencia y cuantificador. No convierte estudios independientes en un único
experimento.

- Use TPM para presencia, perfiles y candidatos.
- Mantenga visible experiment_id.
- No interprete diferencias pequeñas entre estudios como regulación.
- Ejecute DESeq2 dentro de cada estudio balanceado.
- Para metaanálisis, modele estudio/lote y compatibilidad biológica.

## 9. Recuperación

Los constructores no sobrescriben destinos. Si una etapa falla queda una
carpeta .building para inspección. No la borre hasta identificar la causa;
muévala a respaldo antes de reintentar.
