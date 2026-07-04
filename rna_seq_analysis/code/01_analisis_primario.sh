#!/usr/bin/bash

#-----------------------------------------------------------------------------------------------------------------------------------------------
# Script name:  01_analisis_primario.sh

# Description:  Este script consiste en un pipeline automatizado de análisis primario de RNA-seq a partir de un dataset obtenido
#               del repositorio OSDR (Open Science Data Repository) de la NASA. El código de identificación del conjunto de datos es
#               OSD-635. El flujo que realiza este código consiste en la descarga de los archivos crudos de secuenciación R1 y R2
#               de cada muestra para luego hacer un recorte y filtrado por calidad, pseudo-alineamiento y cuantificación. 
#               Finalmente los archivos crudos se eliminan para ahorrar espacio en el disco debido a las limitaciones de la máquina virtual. 
#               Una vez analizados los datos de la primera muestra, se itera este proceso para todas las demás.

# Author:       Laurentiu Jalba Staver
# Institution:  Universidad Internacional de Valencia
# Date:         07/2026
# Version:      1.0
# Usage:        chmod +x 01_analisis_primario.sh
#               ./01_analisis_primario.sh
# Dependencies: wget, fastqc, fastp, salmon, multiqc

# Tener creada previamente la siguiente estructura de directorios:
# rna_seq_analysis/code (Alojar el script aquí)
# rna_seq_analysis/data/reference_transcriptome (Alojar el transcriptoma previamente indexado con el siguiente nombre: index_transcripts_hg38_v49)

#-------------------------------------------------------------------------------------------------------------------------------------------------

# PARADA DE SCRIPT EN CASO DE ERROR
set -e

# LÍNEA DE ERROR
trap 'echo "ERROR: El script ha fallado en la LÍNEA $LINENO"' ERR

echo "===================================================================================================================="
echo "                                 INICIO DEL ANÁLISIS DEL DATASET OSD-635 "
echo "===================================================================================================================="
echo ""

# 1. DEFINICIÓN DE RUTAS DINÁMICAS
script_directory="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_directory="$(dirname "$script_directory")"

# 2. CREACIÓN DE ESTRUCTURA DE DIRECTORIOS NECESARIOS
mkdir -p "$project_directory"/{data/{raw,processed/{1_quality_control_pre,1_quality_control_post,2_trimming,3_pseudoalignment,4_multiqc}},results}

# 3. DEFINICION DE LA RUTA DONDE SE ENCUENTRA EL TRANSCRIPTOMA DE REFERENCIA PREVIAMENTE INDEXADO (hg38_v49)
transcriptome="$project_directory/data/reference_transcriptome/index_transcripts_hg38_v49"

# 4. DEFINICION DE LOS NOMBRES DE LAS 6 MUESTRAS DEL REPOSITORIO
samples=(
    "GLDS-608_rna-seq_SMC-F-E1_S1"
    "GLDS-608_rna-seq_SMC-F-E2_S2"
    "GLDS-608_rna-seq_SMC-F-E3_S3"
    "GLDS-608_rna-seq_SMC-G-E1_S5"
    "GLDS-608_rna-seq_SMC-G-E2_S6"
    "GLDS-608_rna-seq_SMC-G-E3_S7"
)

# Se establece un contador para ver el número de muestra que se está analizando
total_samples=${#samples[@]}
counter=1

# 5. COMIENZO DEL BUCLE QUE ANALIZARÁ CADA MUESTRA
for sample in "${samples[@]}";
do
	echo "---------------------PROCESANDO MUESTRA $counter DE $total_samples: $sample-----------------------------------"
    echo ""

    # Se definen las rutas de los archivos como variables debido a su uso repetitivo
    R1="$project_directory/data/raw/${sample}_R1_raw.fastq.gz"
    R2="$project_directory/data/raw/${sample}_R2_raw.fastq.gz"
    T1="$project_directory/data/processed/2_trimming/${sample}_trimmed_R1.fastq.gz"
    T2="$project_directory/data/processed/2_trimming/${sample}_trimmed_R2.fastq.gz"

    # Se definen las rutas de los archivos .log que muestran la ejecución de los programas como variables
    log_wget="$project_directory/data/raw/${sample}_descarga.log"
    log_fastqc_pre="$project_directory/data/processed/1_quality_control_pre/${sample}_fastqc.log"
    log_fastp="$project_directory/data/processed/2_trimming/${sample}_fastp.log"
    log_fastqc_post="$project_directory/data/processed/1_quality_control_post/${sample}_fastqc.log"
    log_salmon="$project_directory/data/processed/3_pseudoalignment/${sample}_salmon.log"


	# 5.1 Descarga de los archivos FASTQ crudos
	echo "--> [1/7] Descarga de los archivos FASTQ crudos R1 y R2"
	wget -nv -c "https://osdr.nasa.gov/geode-py/ws/studies/OSD-635/download?file=${sample}_R1_raw.fastq.gz&version=1" -O "$R1" > "$log_wget" 2>&1
	wget -nv -c "https://osdr.nasa.gov/geode-py/ws/studies/OSD-635/download?file=${sample}_R2_raw.fastq.gz&version=1" -O "$R2" >> "$log_wget" 2>&1

	# 5.2 Control de calidad mediante FASTQC
	echo "--> [2/7] Control de calidad mediante FastQC pre-trimming"
	fastqc -o "$project_directory/data/processed/1_quality_control_pre" "$R1" "$R2" > "$log_fastqc_pre" 2>&1

    # 5.3 Preprocesado de las lecturas crudas mediante fastp
    echo "--> [3/7] Limpieza de las lecturas de baja calidad y recorte de adaptadores mediante FastP"
    fastp \
        -i "$R1" \
        -I "$R2" \
        -o "$T1" \
        -O "$T2" \
        --detect_adapter_for_pe \
        --trim_poly_g \
        --trim_poly_x \
        --cut_right --cut_window_size 4 --cut_mean_quality 20 \
        --average_qual 20 \
        -l 30 -w 2 \
        -h "$project_directory/data/processed/2_trimming/${sample}_fastp.html" \
        -j "$project_directory/data/processed/2_trimming/${sample}_fastp.json" > "$log_fastp" 2>&1
    
    # 5.4 Eliminación de archivos crudos originales
    echo "--> [4/7] Eliminación de los archivos FASTQ crudos"
	rm "$R1" "$R2"

    # 5.5 Control de calidad post-trimming
    echo "--> [5/7] Control de calidad mediante FastQC post-trimming"
	fastqc -o "$project_directory/data/processed/1_quality_control_post" "$T1" "$T2" > "$log_fastqc_post" 2>&1

    # 5.6 Cuantificación de transcritos mediante Salmon
    echo "--> [6/7] Cuantificación de los transcritos mediante Salmon"
    salmon quant -i "$transcriptome" -l ISR \
        -1 "$T1" \
        -2 "$T2" \
        -p 2 --validateMappings --gcBias --seqBias \
        -o "$project_directory/data/processed/3_pseudoalignment/salmon_results_${sample}" > "$log_salmon" 2>&1

    # 5.7 Eliminación de archivos procesados
    echo "--> [7/7] Eliminación de los archivos FASTQ procesados"
	rm "$T1" "$T2"

    # 5.8 Se hace una copia de los resultados de salmon al directorio results
    cp "$project_directory/data/processed/3_pseudoalignment/salmon_results_${sample}/quant.sf" \
       "$project_directory/results/${sample}_quant.sf"
    
    echo ""
    echo "---------------------ANÁLISIS DE LA MUESTRA $counter DE $total_samples COMPLETADO: $sample-------------------------"
    ((counter++))

    echo ""
done 

echo ""

# 6. CREACIÓN DE REPORTE MULTIQC FINAL
log_multiqc="$project_directory/data/processed/4_multiqc/multiqc.log"
echo "--> Creación de reporte MultiQC final"
multiqc "$project_directory/data/processed" -n multiqc_report_OSD635 -o "$project_directory/data/processed/4_multiqc" > "$log_multiqc" 2>&1

# 6.1 Se hace una copia del reporte MULTIQC al directorio results
cp "$project_directory/data/processed/4_multiqc/multiqc_report_OSD635.html" "$project_directory/results/"

echo "===================================================================================================================="
echo "                                 ANÁLISIS DEL DATASET OSD-635 COMPLETADO "
echo "===================================================================================================================="
