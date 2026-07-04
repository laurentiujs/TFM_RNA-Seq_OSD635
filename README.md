# Análisis transcriptómico del dataset OSD-635

Repositorio con el código utilizado y los resultados correspondientes al Trabajo Fin de Máster titulado *Microgravedad y remodelado vascular: Análisis transcriptómico de la plasticidad fenotípica
de las células musculares lisas vasculares mediante RNA-Seq* realizado para la **Universidad Internacional de Valencia**. 

Este proyecto contiene un pipeline bioinformático automatizado para el análisis de datos de RNA-seq (bulk) procedentes del repositorio Open Science Data Repository (OSDR) de la NASA
(código de acceso del conjunto de datos: **OSD-635**). El objetivo es identificar diferencias en la expresión génica y rutas metabólicas alteradas bajo condiciones de microgravedad frente a controles en Tierra.

---

## ESTRUCTURA DE DIRECTORIOS

Para garantizar la correcta ejecución de los scripts y la reproducibilidad del entorno, el proyecto debe seguir estrictamente la siguiente jerarquía de directorios. 

*Nota: Los archivos FASTQ crudos y el transcriptoma de referencia no se incluyen en este repositorio web debido a su gran tamaño. El script 01_analisis_primario.sh descarga los archivos FASTQ de forma automática.
No obstante, Salmon requiere que el transcriptoma de referencia esté previamente indexado y situado en la carpeta correspondiente, creado manualmente por el usuario*

**Leyenda del repositorio:**
* **[GitHub]**: Archivo o directorio incluido por defecto al descargar este repositorio.
* **[Auto]**: Directorio generado automáticamente por los scripts durante el análisis.
* **[Manual]**: Elemento que el usuario debe aportar o crear antes de iniciar el pipeline.

```text
rna_seq_analysis/
├── code/									<-- [GitHub]
│   ├── 01_analisis_primario.sh               
│   └── 02_analisis_transcriptomico.R         
├── data/
│   ├── metadata/							<-- [GitHub]
│   │   └── OSD-635_metadata_OSD-635-ISA/
│   │       ├── s_OSD-635.txt                 
│   │       └── (resto de archivos ISA) 
│   ├── raw/								<-- [Auto] Descargas temporales de FASTQ crudos
│   ├── processed/							<-- [Auto] Directorio de trabajo intermedio
│   │   ├── 1_quality_control_pre/			<-- [Auto] Outputs de FastQC pre-trimming
│   │   ├── 1_quality_control_post/			<-- [Auto] Outputs de FastQC post-trimming
│   │   ├── 2_trimming/						<-- [Auto] Outputs y logs de fastp
│   │   ├── 3_pseudoalignment/				<-- [Auto] Archivos de Salmon
│   │   └── 4_multiqc/						<-- [Auto] Archivos MultiQC
│   └── reference_transcriptome/			<-- [Manual] Directorio a crear por el usuario
│       └── index_transcripts_hg38_v49/		<-- [Manual] Transcriptoma previamente indexado
└── results/								<-- [GitHub] Al volver a ejecutar los scripts, se sobreescriben
    ├── multiqc_report_OSD635.html			<-- [GitHub] Reporte interactivo global de MultiQC
    ├── *_quant.sf							<-- [Auto] Generados por el script 01_analisis_primario.sh
    └── R_results/							<-- [GitHub] Resultados de R
        ├── Grafico_MDS.png                   
        ├── Grafico_BCV.png                   
        ├── Grafico_QLDisp.png                
        ├── Grafico_Volcano.png               
        ├── Grafico_Correlacion.png           
        ├── Grafico_Heatmap_Top50.png         
        ├── Grafico_GO_*.png             
        ├── Grafico_KEGG_*.png                
        ├── Grafico_Pathview_*.png            
        ├── Grafico_GSEA_*.png                
        ├── Resultados_*.csv				<-- Tablas 01 a 12 con datos de expresión
        └── Resultados_sessionInfo.txt		<-- Información de la sesión de R

```

## PREPARACIÓN PREVIA

1. Descargar este repositorio
2. Crear manualmente el directorio data/reference_transcriptome/ y alojar aquí el transcriptoma de referencia previamente indexado con el nombre index_transcripts_hg38_v49


## PROGRAMAS NECESARIOS

1. Herramientas para el análisis primario en bash: wget, FastQC, fastp, Salmon, MultiQC
2. Herramientas para el análisis secundario en R: Rstudio y los paquetes de R tximport, edgeR, EnsDb.Hsapiens.v86, org.Hs.eg.db, ggplot2, pheatmap, clusterProfiler, ggrepel, enrichplot, pathview.

## INSTRUCCIONES DE EJECUCIÓN

1. Análisis primario en bash (01_analisis_primario.sh)
	Otorgar permisos de ejecución al script: chmod +x 01_analisis_primario.sh
	Ejecutar el script: ./01_analisis_primario.sh
2. Análisis de expresión diferencial y enriquecimeinto funcional en R (02_analisis_transcriptomico.R)
	Ejecutar el script desde la consola de R: source("code/02_analisis_transcriptomico.R")

## AUTORÍA

Laurentiu Jalba Staver
Máster Universitario en Bioinformática por la Universidad Internacional de Valencia
Julio 2026
