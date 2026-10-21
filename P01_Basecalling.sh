#!/usr/bin/env bash
#SBATCH -J P01_Basecalling                          # Job name
#SBATCH -o P01_Basecalling.%a.%A.out               # File to which stdout will be written
#SBATCH -e P01_Basecalling.%a.%A.err               # File to which stderr will be written
#SBATCH -p gpu_all                                # Partition
#SBATCH -N 1                                      # Number of node - 1 Milan node has 32-core, 500 Go RAM and GPU Ampere A100
#SBATCH --gres gpu:8                              # Number GPU
#SBATCH --ntasks-per-node 1                       # Number of cores/cpus
#SBATCH --cpus-per-task 28                        # Number of cpus per task
#SBATCH -t 00-20:00                               # Runtime in DD-HH:MM
#SBATCH --mem 100G                                # Memory for all cores in Gbytes
#SBATCH --mail-type=ALL                           # BEGIN,END,ALL
#SBATCH --mail-user=nicolas.nesi@unicaen.fr       # Email address

# ---------------------------------
# Version v1.0
# ---------------------------------
# authors: Nicolas Nesi
# University Caen Normandy, DYNAMICURE INSERM UMR 1311
# Date: 19/10/2026
# ---------------------------------

# ---------------------------------
# Variables and Paths
# ---------------------------------
# safety for failing steps
set -euo pipefail

# conda swtich environments
conda_switch() {
    set +u
    conda deactivate 2>/dev/null || true
    conda activate "$1"
    set -u
}

# Positional parameters
BASECALLER=""
GROUP="ViroCaen"
PROJECT="HumanGenetic"
RUN=""
MODEL=""
BARCODE="SQK-NBD114-24" 
QSCOREBASECALL="15"

usage() {
   local exit_code="${1:-1}"
   {
    echo "Usage: sbatch $0 -b BASECALLER -r RUN [-m MODEL] [-br BARCODE] [-qs QSCORE]"
	echo "   -b, --basecaller  Mandatory. Either Guppy or Dorado"
	echo "   -r, --run    Mandatory. Name of run directory"
	echo "   -m, --model Optional. Name of the model for basecalling"
	echo "   -br, --barcode Optional. Name of the barcode kit [default SQK-NBD114-24]"
	echo "   -qs, --qscore Optional. Qscore value [default 15]"
	echo "   -h, --help   Show this help"
    } >&2
	exit "$exit_code"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
	-b|--basecaller) BASECALLER="$2"; shift 2 ;;
    -r|--run) RUN="$2"; shift 2 ;;
	-m|--model) MODEL="$2"; shift 2 ;;
	-qs|--qscore) QSCOREBASECALL="$2"; shift 2 ;;
	-br|--barcode) BARCODE="$2"; shift 2 ;;
	-h|--help)  usage 0 ;;
    *) 
         echo "Invalid option: $1" >&2
         usage 1
         ;;
    esac
done

# Validation du choix de basecaller
if [[ "${BASECALLER}" != "guppy" && "${BASECALLER}" != "dorado" ]]; then
    echo "Erreur: BASECALLER doit être 'guppy' ou 'dorado', reçu: '${BASECALLER}'" >&2
    exit 1
fi

if [[ "$BASECALLER" == "guppy" ]]; then
    if [[ -z "$MODEL" ]]; then
        MODEL="dna_r9.4.1_450bps_hac.cfg"
        echo "MODEL not specified, using Guppy default: $MODEL"
    elif [[ "$MODEL" != *.cfg ]]; then
        echo -e "\e[31mError: BASECALLER=guppy requires a Guppy config file (.cfg), got: '$MODEL'\e[0m"
        exit 1
    fi
elif [[ "$BASECALLER" == "dorado" ]]; then
    if [[ -z "$MODEL" ]]; then
        MODEL="dna_r10.4.1_e8.2_400bps_sup@v5.2.0"
        echo "MODEL not specified, using Dorado default: $MODEL"
    elif [[ "$MODEL" != *@* ]]; then
        echo -e "\e[31mError: BASECALLER=dorado requires a Dorado model name (containing '@'), got: '$MODEL'\e[0m"
        exit 1
    fi
fi

[[ -z "$RUN" ]] && { echo "Option -r obligatoire" >&2; usage 1; }

# Paths
FOLDER="/dlocal/home/2019013/Data/$GROUP/$PROJECT"
SCRIPTS="$FOLDER/Scripts"
PATH2="$FOLDER/RawReads_Nanopore/$RUN"

if [[ "${BASECALLER}" == "guppy" ]]; then
    BASECALL_DIR="${FOLDER}/Guppy_Basecalling/${RUN}"
else
    BASECALL_DIR="${FOLDER}/Dorado_Basecalling/${RUN}"
fi

# Variables
RUNDATE="$(date '+%Y-%m-%d')"

# Threads
export THREADS="${SLURM_CPUS_PER_TASK}"
echo "number of threads used ${THREADS}"
# ---------------------------------

# ---------------------------------
# Sample table: Excel -> Dorado sample sheet or Guppy barcode table
# ---------------------------------
echo "Formatting Exel file for Basecalling"

SAMPLES_SAMPLES_XLSX="${SCRIPTS}/Samples_${RUN}.xlsx"
SAMPLES_DORADO="${SCRIPTS}/Samples_Dorado_${RUN}.csv"
SAMPLES_GUPPY="${SCRIPTS}/Samples_Guppy_${RUN}.csv"

[[ -f "$SAMPLES_XLSX" ]] || { echo "ERROR: ${SAMPLES_XLSX} not found" >&2; exit 1; }

module load py_env/miniconda3/25.7.0
set +u
eval "$(conda shell.bash hook)"
set -u
conda_switch NN_env
# include:
# csvtk=0.31.0

# Keep the 5 required columns (by name) and drop rows without alias
csvtk xlsx2csv "$SAMPLES_XLSX" \
    | csvtk cut -f flow_cell_id,kit,experiment_id,barcode,alias \
    | awk -F',' 'NR==1 || $5 !~ /^[[:space:]]*$/' > "$CLEAN"

set +u; conda deactivate; set -u

# --- Checks ---
[[ "$(tail -n +2 "$CLEAN" | wc -l)" -gt 0 ]] \
    || { echo "ERROR: no sample with an alias in ${SAMPLES_XLSX}" >&2; exit 1; }

DUP="$(tail -n +2 "$CLEAN" | cut -d',' -f4 | sort | uniq -d)"
[[ -z "$DUP" ]] || { echo "ERROR: duplicated barcode(s): $DUP" >&2; exit 1; }

DUP="$(tail -n +2 "$CLEAN" | cut -d',' -f5 | sort | uniq -d)"
[[ -z "$DUP" ]] || { echo "ERROR: duplicated alias(es): $DUP" >&2; exit 1; }

BAD="$(tail -n +2 "$CLEAN" | cut -d',' -f5 | grep -vE '^[A-Za-z0-9._-]+$' || true)"
[[ -z "$BAD" ]] || { echo "ERROR: invalid alias (use letters, digits, . _ -): $BAD" >&2; exit 1; }

KITS="$(tail -n +2 "$CLEAN" | cut -d',' -f2 | sort -u)"
[[ "$KITS" == "$BARCODE" ]] \
    || { echo "ERROR: kit in table (${KITS}) differs from -br (${BARCODE})" >&2; exit 1; }

# --- Format for the chosen basecaller ---
if [[ "$BASECALLER" == "dorado" ]]; then
    cp "$CLEAN" "${SAMPLE_DORADO}"
else
    tail -n +2 "$CLEAN" | cut -d',' -f4,5 > "${SCRIPTS}/${SAMPLES_GUPPY}"
    tail -n +2 "$CLEAN" | cut -d',' -f5   > "${SCRIPTS}/List_samples_${RUN}.txt"
fi

rm -f "$CLEAN"
# ---------------------------------

# ---------------------------------
# Copy input data to LOCAL_WORK_DIR
# ---------------------------------
if [[ "${BASECALLER}" == "guppy" ]]; then
    cp -R "${FOLDER}"/Nanopore_outputRun/"${RUN}"/no_sample_id/*/fast5 "${LOCAL_WORK_DIR}"
else
    cp -R "${FOLDER}"/Nanopore_outputRun/"${RUN}"/no_sample_id/*/pod5_skip "${LOCAL_WORK_DIR}"
    mkdir -p "${LOCAL_WORK_DIR}"/demultiplex
fi

if [[ "${BASECALLER}" == "guppy" ]]; then
    if [[ ! -d "${LOCAL_WORK_DIR}/fast5" ]]; then
        echo "ERROR: fast5 directory not found in ${LOCAL_WORK_DIR}" >&2
        exit 1
    fi
else
    if [[ ! -d "${LOCAL_WORK_DIR}/pod5_skip" ]]; then
        echo "ERROR: pod5_skip directory not found in ${LOCAL_WORK_DIR}" >&2
        exit 1
    fi
fi

cd "${LOCAL_WORK_DIR}"/

echo "Analysis run the: ${RUNDATE}"
echo "Working directory: ${PWD}"
echo "Submitted used:"
echo "Basecaller ${BASECALLER} Group ${GROUP} Project ${PROJECT} Run ${RUN}"
echo "Config/Model ${MODEL} NanoporeKit ${BARCODE} Qscore ${QSCOREBASECALL}"
# ---------------------------------

# ---------------------------------
# Make subdirectories
# ---------------------------------
for dir in \
    "${BASECALL_DIR}" \
    "${FOLDER}/RawReads_Nanopore/${RUN}" \
    "${FOLDER}/CleanedReads_Nanopore/${RUN}" \
    "${FOLDER}/minimap2/${RUN}" \
    "${FOLDER}/QC_Stats/${RUN}/PycoQC" \
    "${FOLDER}/QC_Stats/${RUN}/BUSCO" \
    "${FOLDER}/QC_Stats/${RUN}/Seqkit" \
    "${FOLDER}/QC_Stats/${RUN}/RawReads_Nanopore" \
    "${FOLDER}/QC_Stats/${RUN}/CleanedReads_Nanopore" \
    "${FOLDER}/Log/${RUN}_run${RUNDATE}" \
    "${FOLDER}/Scripts/Sbatch_log"
do
    if [ ! -d "$dir" ]; then
        mkdir -p "$dir"
        echo "Created directory: $dir"
    else
        echo "Directory already exists: $dir"
    fi
done
# ---------------------------------

# ---------------------------------
# Basecalling
# ---------------------------------
if [[ "${BASECALLER}" == "guppy" ]]; then
    # =========================================
    # GUPPY
    # =========================================
    module load bio/guppy/6.5.7_gpu

    echo "Running Basecalling, Trimming (barcodes, adapters) and Demultiplexing using Guppy"

    guppy_basecaller \
        -i "${LOCAL_WORK_DIR}"/fast5 \
        -s "${LOCAL_WORK_DIR}"/demultiplex \
        -c "${MODEL}" \
        -x "cuda:all" \
        --num_callers 28 \
        --gpu_runners_per_device 8 \
        --chunks_per_runner 512 \
        --chunk_size 3000 \
        --trim_strategy dna \
        --trim_adapters \
        --detect_mid_strand_barcodes \
        --enable_trim_barcodes \
        --barcode_kits "${BARCODE}" \
        --min_qscore "${QSCOREBASECALL}"

    # --- Concatenate pass fastq ---
    if [[ -d demultiplex/pass/unclassified ]]; then
        mv demultiplex/pass/unclassified demultiplex/pass/barcode00
    fi

    for barcode_dir in demultiplex/pass/barcode* ; do
        if [[ -d "${barcode_dir}" ]]; then
            cat "${barcode_dir}"/*.fastq > "${barcode_dir}".fastq
            rm -r "${barcode_dir}"
            mv "${barcode_dir}".fastq demultiplex/
        fi
    done

    # --- Move demultiplex outputs ---
	if [ -d "${BASECALL_DIR}/demultiplex" ]; then
        rm -rf "${BASECALL_DIR}/demultiplex"
	fi
    mv demultiplex/ "${BASECALL_DIR}"
    cd "${BASECALL_DIR}"/demultiplex/
    
    cp "${SAMPLES_GUPPY}" "${BASECALL_DIR}"/demultiplex/

    cd "${BASECALL_DIR}"/demultiplex/

    awk -F',' '{ print $1".fastq", $2".fastq" }' "Samples_Guppy_${RUN}.csv" > OldNewNames.txt

    while read -r old new; do
        mv "$old" "$new"
    done < OldNewNames.txt

    # --- Move Fastq to RawReads ---
    mv "${BASECALL_DIR}"/demultiplex/*.fastq "${PATH2}"/
    rm -f "${PATH2}"/barcode*.fastq
    pigz -f "${PATH2}"/*.fastq

elif [[ "${BASECALLER}" == "dorado" ]]; then

    # =========================================
    # DORADO
    # =========================================

    module load bio/dorado/1.2.0 # back to 0.9.6 if needed

    echo "Running Basecalling, Demultiplexing and Adapter/barcode trimming using Dorado"

    dorado basecaller \
        /soft/2019013/Logiciels/Dorado-models/"${MODEL}" \
        "${LOCAL_WORK_DIR}"/pod5_skip/ \
        --device "cuda:all" \
        --min-qscore "${QSCOREBASECALL}" \
        --chunksize 10000 \
        --no-trim \
        --verbose > "${RUN}".basecalls.bam

    dorado summary "${RUN}".basecalls.bam > demultiplex/sequencing_summary.txt

    dorado demux \
        --kit-name "${BARCODE}" \
        --threads "${THREADS}" \
        --sample-sheet "${SAMPLES_DORADO}" \
        --emit-summary \
        --output-dir demultiplex/ "${RUN}".basecalls.bam

    for f in demultiplex/*.bam; do
        base=$(basename "$f" .bam)
        dorado trim \
            --no-trim-primers \
            --sequencing-kit "${BARCODE}" \
            --threads "${THREADS}" \
            "$f" > demultiplex/"${base}"-trimmed.bam
    done

    # --- Environments ---
    module load py_env/miniconda3/25.7.0
	set +u
    eval "$(conda shell.bash hook)"
    conda activate clean-ONT_env
	set -u
    # include:
    # samtools=1.18

    # --- Convert uBAM to Fastq ---
    echo "Convert uBAM to Fastq"
    for f in demultiplex/*-trimmed.bam; do
        base=$(basename "$f" -trimmed.bam)
        samtools bam2fq -@ "${THREADS}" "$f" > demultiplex/"${base}".fastq
    done

    # --- Move and rename fastq ---
    echo "Move and rename fastq"

    if [ -d "${BASECALL_DIR}/demultiplex" ]; then
        rm -rf "${BASECALL_DIR}/demultiplex"
	fi
    mv demultiplex/ "${BASECALL_DIR}"
    mv "${RUN}".basecalls.bam "${BASECALL_DIR}"

    cd "${BASECALL_DIR}"/demultiplex/

    for f in *.fastq; do
        new_name=$(echo "$f" | sed "s/^.*${BARCODE}_//")
        [[ "$f" == "$new_name" ]] || mv "$f" "$new_name"
    done
	
	[[ -f unclassified.fastq ]] && mv unclassified.fastq NoBarcode.fastq || true

    mv "${BASECALL_DIR}"/demultiplex/*.fastq "${PATH2}"/
    pigz -f "${PATH2}"/*.fastq
	
	for f in "${PATH2}"/*.fastq.gz; do
      basename "$f" .fastq.gz
    done | grep -v '^NoBarcode$' > "${SCRIPTS}/List_samples_${RUN}.txt"
fi
# ---------------------------------

# ---------------------------------
# Quality control with PycoQC
# ---------------------------------
conda_switch QC-ONT_env
# include:
# PycoQC=2.5.2

if [[ "${BASECALLER}" == "guppy" ]]; then

    pycoQC \
        --summary_file "${BASECALL_DIR}/demultiplex/sequencing_summary.txt" \
        --report_title "${PROJECT}_${RUN}_${MODEL}_qscore${QSCOREBASECALL}_report" \
        --html_outfile "${FOLDER}/QC_Stats/${RUN}/PycoQC/${PROJECT}_${RUN}_run${RUNDATE}.html"

elif [[ "${BASECALLER}" == "dorado" ]]; then

    cd "${BASECALL_DIR}"/demultiplex/

    # format summary files for pycoQC
    { head -n 1 sequencing_summary.txt; tail -n +2 sequencing_summary.txt | sort -k2,2; } > sequencing_summary_sorted.txt
    { head -n 1 barcoding_summary.txt; tail -n +2 barcoding_summary.txt | sort -k2,2; } > barcoding_summary_sorted.txt
    header=$(head -n 1 sequencing_summary_sorted.txt)
    awk 'NR==FNR {barcodes[FNR]=$NF; next} FNR > 1 {$NF=barcodes[FNR]; print}' barcoding_summary_sorted.txt sequencing_summary_sorted.txt > dorado_temp
    echo "$header" | cat - dorado_temp > dorado_temp2
    awk '{$1=$1}1' OFS='\t' dorado_temp2 > dorado_summary.tsv
    rm dorado_temp*

    pycoQC \
        --summary_file "${BASECALL_DIR}/demultiplex/dorado_summary.tsv" \
        --report_title "${PROJECT}_${RUN}_${MODEL}_qscore${QSCOREBASECALL}_report" \
        --html_outfile "${FOLDER}/QC_Stats/${RUN}/PycoQC/${PROJECT}_${RUN}_run${RUNDATE}.html"
fi
# ---------------------------------

# ---------------------------------
# Cleaning Working Directory
# ---------------------------------
rm -rf "${LOCAL_WORK_DIR:?}"/*
# ---------------------------------

# ---------------------------------
mv "${SLURM_SUBMIT_DIR}/P01_Basecalling.${SLURM_JOB_ID}."* "${FOLDER}/Scripts/Sbatch_log/" 2>/dev/null || true
# ---------------------------------