# Peipers-RNAseq docente para tomate

## Autoría y contacto

Esta guía y el pipeline Peipers-RNAseq fueron creados por **Samuel Parra A.,
PhD**.

- Correo institucional: [sa.parra@uandresbello.edu](mailto:sa.parra@uandresbello.edu)
- Filiación actual: Centro de Biotecnología Molecular Vegetal, Departamento de
  Ciencias, Universidad de Chile.
- Investigador postdoctoral en el laboratorio de la Dra. Claudia Stange.
- ORCID: [0000-0002-9129-4133](https://orcid.org/0000-0002-9129-4133)

Copyright (c) 2026 Samuel Parra A., PhD. Código y documentación distribuidos
bajo licencia MIT. Consulte [AUTHORSHIP.md](AUTHORSHIP.md),
[CITATION.cff](CITATION.cff) y [LICENSE](LICENSE).

Esta carpeta permite que un alumno con poca experiencia procese datos RNA-seq
de tomate desde ENA/SRA hasta un consultor TPM, aprenda qué hace cada etapa y
entregue un resultado incorporable a
[TPM-consultor](https://github.com/SamuelParraA/TPM-consultor).

El caso preparado es
[PRJNA1347747](https://www.ebi.ac.uk/ena/browser/view/PRJNA1347747), asociado al
trabajo *Low ethylene production in the root(stock) alleviates salt stress in
tomato*. También se incluye una herramienta para iniciar otros BioProjects.

## Productos finales

- 18 cuantificaciones Salmon, una por corrida SRR.
- Matrices de TPM y lecturas estimadas a nivel de transcrito y gen.
- Un informe MultiQC.
- Un paquete verificable apto para un pull request.
- Un consultor standalone con su propio SQLite.
- Opcionalmente, contrastes de expresión diferencial con DESeq2.

Los FASTQ ocupan aproximadamente 92,5 GiB y quedan fuera de GitHub.

## Conceptos que debe aprender el alumno

- **BioProject**: agrupación pública de muestras y archivos de un estudio.
- **SRR**: identificador de una corrida de secuenciación.
- **FASTQ**: secuencias leídas por el equipo y calidad de cada base.
- **Referencia transcriptómica**: transcritos conocidos contra los que se
  cuantifican las lecturas.
- **Salmon**: estima cuántas lecturas proceden de cada transcrito.
- **TPM**: abundancia normalizada dentro de una biblioteca. Sirve para explorar
  perfiles; no reemplaza conteos ni DESeq2.
- **MD5/SHA-256**: huellas digitales para detectar archivos incompletos o
  diferentes.
- **FastQC/MultiQC**: diagnóstico de calidad; no cambian las lecturas.
- **DESeq2**: modelo estadístico de conteos para condiciones con réplicas.

## Diseño experimental curado

Todas las bibliotecas son RNA-seq paired-end de hoja, 100 nt, BGISEQ-500.

| Vástago/portainjerto | Tratamiento | Réplicas SRR |
|---|---|---|
| UD/WT | Control | SRR35859425, SRR35859424, SRR35859415 |
| UD/WT | 75 mM NaCl | SRR35859408, SRR35859423, SRR35859422 |
| UD/epi | Control | SRR35859414, SRR35859413, SRR35859412 |
| UD/epi | 75 mM NaCl | SRR35859421, SRR35859420, SRR35859419 |
| UD/ACCD | Control | SRR35859411, SRR35859410, SRR35859409 |
| UD/ACCD | 75 mM NaCl | SRR35859418, SRR35859417, SRR35859416 |

El manifiesto completo incluye BioSample, experimento, URL, MD5 y tamaño:
**manifests/PRJNA1347747_runs.tsv**.

## Por qué estos TPM son compatibles

El artículo original usó STAR y featureCounts. Para comparar con los TPM
históricos del consultor, esta ruta vuelve a cuantificar con:

- ITAG4.0_cDNA.fasta e ITAG4.0_gene_models.gff;
- Salmon 2.0.1;
- índice k=31;
- detección de orientación -l A;
- lecturas paired-end comprimidas de ENA;
- sin seqBias, gcBias ni recorte por defecto.

Las referencias se descargan del
[release oficial ITAG4.0 de SGN](https://solgenomics.net/ftp/tomato_genome/annotation/ITAG4.0_release/)
y se exigen estas huellas:

~~~text
ITAG4.0_cDNA.fasta       2956f6f4c76c09eaef65d8732d9fa02d3682a0bff35a3b5027473d224f7beb32
ITAG4.0_gene_models.gff  61bdf7c0a607c723fa75c3b3139690489aedf23c0dcc230e15050b6edfc9f820
~~~

> Comparable técnicamente no significa libre de efectos de lote. Cultivar,
> edad, tejido, plataforma y laboratorio aún influyen. No use DESeq2 mezclando
> BioProjects no balanceados.

## Flujo

~~~mermaid
flowchart LR
    A["ENA: pares FASTQ"] --> B["MD5 + FastQC/MultiQC"]
    R["ITAG4.0 con SHA-256 fijo"] --> I["Índice Salmon 2.0.1"]
    B --> Q["Salmon: quant.sf"]
    I --> Q
    Q --> M["Matrices transcrito/gen"]
    M --> G["Paquete para pull request"]
    Q --> C["Consultor standalone"]
    Q --> D["DESeq2 opcional"]
~~~

## 1. Preparar Windows y Ubuntu/WSL2

### 1.1 Revisar recursos

En Administrador de tareas > Rendimiento anote procesadores lógicos, RAM y
espacio libre.

- 8 procesadores lógicos o más.
- 32 GB RAM ideal; con 16 GB use 4 hilos.
- 130 GB libres como mínimo; 150-200 GB dan margen.

No use OneDrive para los datos grandes.

### 1.2 Instalar WSL2

Abra PowerShell como administrador:

~~~powershell
wsl --install -d Ubuntu-24.04
~~~

Reinicie si se solicita. Abra Ubuntu, cree usuario y contraseña. Al escribir
la contraseña no se muestran caracteres; es normal.

~~~bash
sudo apt update
sudo apt install -y git wget bzip2 ca-certificates
nproc
free -h
df -h
~~~

### 1.3 Instalar Miniforge

Conda crea ambientes aislados y evita mezclar versiones.

~~~bash
cd ~
wget -O Miniforge3.sh https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3.sh
~~~

Acepte la licencia, use la ruta sugerida y permita inicializar el shell. Cierre
y vuelva a abrir Ubuntu:

~~~bash
conda --version
~~~

## 2. Crear fork y rama

1. Abra TPM-consultor en GitHub.
2. Pulse **Fork**.
3. Reemplace USUARIO_GITHUB:

~~~bash
cd ~
git clone https://github.com/USUARIO_GITHUB/TPM-consultor.git
cd TPM-consultor
git switch -c dataset/prjna1347747
git config --global user.name "Nombre Apellido"
git config --global user.email "correo-asociado-a-github"
~~~

Clone dentro de la carpeta Linux, no en /mnt/c/.../OneDrive.

## 3. Crear el ambiente

~~~bash
cd ~/TPM-consultor
conda config --set channel_priority strict
conda env create -f peipers_rnaseq/environment.yml
conda activate peipers-rnaseq-tomate
python --version
salmon --version
fastqc --version
multiqc --version
~~~

Salmon debe informar 2.0.1. No actualice programas durante el análisis.

## 4. Configurar sin cambiar archivos versionados

~~~bash
cp peipers_rnaseq/config/PRJNA1347747.env peipers_rnaseq/config/local.env
nano peipers_rnaseq/config/local.env
~~~

Use THREADS=4 con 16 GB RAM y THREADS=8 con 32 GB o más. Cambie WORK_DIR si
otro disco tiene más espacio. Ejemplo:

~~~bash
WORK_DIR="/mnt/d/rnaseq-work/PRJNA1347747"
~~~

El disco de Windows puede ser más lento, pero es preferible a llenar C:. No
use una ruta de OneDrive.

~~~bash
cd ~/TPM-consultor/peipers_rnaseq
bash run_pipeline.sh preflight config/local.env
~~~

El preflight se detiene si falta software, la metadata es inconsistente,
Salmon no es 2.0.1 o falta espacio.

## 5. Referencia e índice

~~~bash
bash run_pipeline.sh reference config/local.env
~~~

La etapa descarga cDNA/GFF, verifica SHA-256, construye el índice y comprueba
versión, k y 33.976 transcritos. Es reanudable.

## 6. Descargar FASTQ

~~~bash
bash run_pipeline.sh download config/local.env
~~~

Hay 36 archivos. La descarga se reanuda si se corta Internet. Cada archivo se
acepta sólo si tamaño y MD5 coinciden con ENA. Si Windows reinicia, ejecute el
mismo comando.

No edite ni copie FASTQ al repositorio.

## 7. Control de calidad

~~~bash
bash run_pipeline.sh qc config/local.env
~~~

Abra WORK_DIR/qc/multiqc/multiqc_report.html. Revise:

- que estén todos los SRR;
- calidad por posición;
- adaptadores intensos;
- contenido GC anómalo;
- una biblioteca claramente distinta en calidad o profundidad.

No se recorta automáticamente para conservar comparabilidad. Si el QC muestra
un problema grave, detenga el análisis y consulte al supervisor.

## 8. Cuantificar

~~~bash
bash run_pipeline.sh quant config/local.env
~~~

El comando conceptual es:

~~~bash
salmon quant -i INDICE -l A -1 MATE_1.fastq.gz -2 MATE_2.fastq.gz -p HILOS -o SALIDA
~~~

Cada resultado debe tener 33.976 transcritos, suma TPM cercana a un millón,
Salmon 2.0.1 y al menos 40 % de mapeo. Una salida incompleta nunca se
sobrescribe silenciosamente.

## 9. Construir matrices

~~~bash
bash run_pipeline.sh matrices config/local.env
~~~

Se generan:

- Transcript_TPM.tsv y Transcript_NumReads.tsv.
- Gene_TPM.tsv y Gene_NumReads.tsv.
- tx2gene.tsv.
- salmon_qc.tsv.

Los TPM de gen son sumas de sus transcritos según ITAG4.0. tx2gene registra la
agregación y permite usar tximport.

## 10. Consultor standalone

~~~bash
bash run_pipeline.sh consultor config/local.env
source config/local.env
cd "$WORK_DIR/consultor_standalone/PRJNA1347747"
python app.py
~~~

Abra [http://127.0.0.1:8765](http://127.0.0.1:8765). Pruebe:

~~~text
Solyc07g049530
Solyc12g005940
Solyc02g084850
~~~

Debe mostrar seis condiciones con tres réplicas.

## 11. Paquete y pull request

~~~bash
cd ~/TPM-consultor/peipers_rnaseq
bash run_pipeline.sh bundle config/local.env
bash run_pipeline.sh verify config/local.env
source config/local.env
cd ~/TPM-consultor
mkdir -p datasets
cp -a "$WORK_DIR/bundle/PRJNA1347747" datasets/
python peipers_rnaseq/scripts/verify_bundle.py --bundle datasets/PRJNA1347747
git status --short
~~~

El paquete incluye metadata, quant.sf, matrices, QC y checksums. Excluye
FASTQ, SRA, referencia, índice y SQLite.

Sólo debe aparecer datasets/PRJNA1347747. Si aparecen datos grandes, deténgase.

~~~bash
git add datasets/PRJNA1347747
git diff --cached --stat
git commit -m "Agregar PRJNA1347747 procesado con Peipers-RNAseq"
git push -u origin dataset/prjna1347747
~~~

En GitHub pulse **Compare & pull request**. Incluya BioProject, diseño, ITAG4.0,
Salmon 2.0.1, confirmación MD5, rango de mapeo, MultiQC y resultado del
verificador. No haga push directo a main.

## 12. DESeq2 opcional

DESeq2 no analiza TPM; importa las estimaciones Salmon con tximport.

~~~bash
cd ~/TPM-consultor
conda env create -f peipers_rnaseq/optional_deseq2/environment-deseq2.yml
conda activate peipers-rnaseq-deseq2
cd peipers_rnaseq
bash run_pipeline.sh deseq2 config/local.env
~~~

El diseño es genotype, treatment y su interacción. Produce sal vs control
dentro de WT/epi/ACCD, comparaciones de portainjertos, interacciones, conteos
normalizados, PCA, distancias y sessionInfo.txt.

Interprete padj < 0.05 junto con log2FoldChange. Una interacción significativa
indica que la respuesta a sal difiere respecto de WT. Use estos contrastes sólo
dentro de PRJNA1347747.

## 13. Otro BioProject de tomate

~~~bash
conda activate peipers-rnaseq-tomate
python scripts/fetch_ena_manifest.py \
  --project PRJNA_OTRO \
  --output manifests/PRJNA_OTRO_runs.draft.tsv
~~~

Después:

1. Lea el artículo y los BioSamples.
2. Reemplace todos los campos REVISAR.
3. Cree PRJNA_OTRO_samples.tsv.
4. Copie la configuración y cambie proyecto, nombres y WORK_DIR.
5. Ejecute validate_project.py.
6. Procese con la misma referencia y Salmon.

Este camino integra datasets de tomate cuantificables contra ITAG4.0. Otra
especie requiere otra referencia y un consultor separado.

## Problemas frecuentes

### Espacio insuficiente

Cambie WORK_DIR. No reduzca MIN_FREE_GB para ignorar el problema.

### MD5 incorrecto

No use el archivo. Conserve el mensaje y mueva el .part sospechoso a otra
carpeta antes de reintentar.

### Salmon no es 2.0.1

~~~bash
conda activate peipers-rnaseq-tomate
conda list salmon
~~~

Si el ambiente fue modificado, recréelo desde environment.yml.

### Windows reinició

Active el ambiente y repita la etapa. El pipeline reutiliza resultados válidos.

### Mapeo bajo 40 %

No cambie el umbral. Revise especie, pares, FastQC y logs con el supervisor.

## Liberar espacio

Primero valide:

~~~bash
bash run_pipeline.sh verify config/local.env
source config/local.env
du -sh "$WORK_DIR/fastq" "$WORK_DIR/quant" "$WORK_DIR/bundle"
readlink -f "$WORK_DIR/fastq"
~~~

Sólo si la ruta resuelta es exactamente la carpeta FASTQ del proyecto:

~~~bash
rm -r -- "$WORK_DIR/fastq"
~~~

Esto libera cerca de 100 GB y no se puede deshacer localmente. Conserve quant,
matrices, bundle y consultor.

## Lista de éxito

- [ ] Trabajé en una rama de mi fork.
- [ ] Preflight terminó sin errores.
- [ ] Las referencias pasaron SHA-256.
- [ ] Los 36 FASTQ pasaron MD5 y tamaño.
- [ ] Revisé MultiQC.
- [ ] Las 18 cuantificaciones pasaron validación.
- [ ] El consultor muestra seis condiciones con tres réplicas.
- [ ] verify_bundle.py terminó con OK.
- [ ] Git no muestra FASTQ, SRA, índices ni SQLite.
- [ ] Abrí un pull request.

La incorporación al consultor unificado está en
[GUIA_MANTENEDOR.md](GUIA_MANTENEDOR.md).
