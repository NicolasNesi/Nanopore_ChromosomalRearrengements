#!/usr/bin/env bash
#SBATCH -J P02b_Clair3                          # Job name 
#SBATCH -o P02b_Clair3.%a.%A.out                # File to which stdout will be written
#SBATCH -e P02b_Clair3.%a.%A.err                # File to which stderr will be written
#SBATCH -p tcourt                               # Partition
#SBATCH -n 1                                    # Number of cores/cpus
#SBATCH -c 24                                   # Number of cpus ask per task
#SBATCH -t 01-00:00                             # Runtime in DD-HH:MM
#SBATCH --mem 100G                              # Memory for all cores in Gbytes (or --mem 752000 in Mbytes)
#SBATCH --mail-type=ALL                         # BEGIN,END,ALL,NONE
#SBATCH --mail-user=nicolas.nesi@unicaen.fr     # Email address

# ---------------------------------
# Version v1.0
# ---------------------------------
# authors: Nicolas Nesi
# University Caen Normandy, DYNAMICURE INSERM UMR 1311
# Date: 20/10/2026
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
RUN=""
MODEL="r1041_e82_400bps_sup_v520"
GROUP="ViroCaen"
PROJECT="HumanGenetic"
REFERENCE="GRCh38_no_alt.fa"
REFNAME="${REFERENCE%.fa}"

usage() {
   local exit_code="${1:-1}"
   {
    echo "Usage: sbatch $0 -r RUN [-m MODEL]"
	echo "   -r, --run    Mandatory. Name of run directory"
	echo "   -m, --model Optional. Name of the model for clair3 [default r1041_e82_400bps_sup_v520]"
	echo "   -h, --help   Show this help"
    } >&2
	exit "$exit_code"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
    -r|--run) RUN="$2"; shift 2 ;;
	-m|--model) MODEL="$2"; shift 2 ;;
	-h|--help)  usage 0 ;;
    *) 
         echo "Invalid option: $1" >&2
         usage 1
         ;;
    esac
done

[[ -z "$RUN" ]] && { echo "Option -r obligatoire" >&2; usage 1; }

# Paths
FOLDER="/dlocal/home/2019013/Data/$GROUP/$PROJECT"
PATH1="${FOLDER}/Mapping/${RUN}"
PATH2="/soft/2019013/Logiciels/Clair3"

# Database
MAPPING="/dlocal/home/2019013/Databases/ReferenceGenomes_Mapping"

# Array job
SAMPLE="$(sed -n -e "${SLURM_ARRAY_TASK_ID}p" ${FOLDER}/Scripts/List_samples_$RUN.txt)"


# Variables
RUNDATE="$(date '+%Y-%m-%d')"

# Threads
export THREADS="${SLURM_CPUS_PER_TASK}"
echo "number of threads used ${THREADS}"
# ---------------------------------

# ---------------------------------
# Copy input data to LOCAL_WORK_DIR
# ---------------------------------
cp "${PATH1}"/"${SAMPLE}_mapped${REFNAME}"_sorted.ba* "${LOCAL_WORK_DIR}"
cp "${MAPPING}/${REFERENCE}" "${MAPPING}/${REFERENCE}.fai" "${LOCAL_WORK_DIR}"
cd "${LOCAL_WORK_DIR}"/

echo "Analysis run the: ${RUNDATE}"
echo "Working directory: ${PWD}"
echo "Variables used:"
echo "group:$GROUP project:$PROJECT run:$RUN model:$MODEL"
# ---------------------------------

# ---------------------------------
# Environments
# ---------------------------------
module load py_env/miniconda3/25.7.0
set +u
eval "$(conda shell.bash hook)"
conda activate sniffles_env
set -u
# include:
# mosdepth=0.3.14
# ont-spectre=0.3.3
# ---------------------------------

# ---------------------------------
# Calling CNV
# ---------------------------------
echo "Calling CNV for duplication deletion events"

samtools coverage -q 7 "${SAMPLE}_mapped${REFNAME}_sorted.bam" > "CoverageTable_${SAMPLE}_to_${REFNAME}.tsv"

MOSDIR="${SAMPLE}_mosdepth"
mkdir -p "${SAMPLE}_mosdepth"

mosdepth \
--by 1000 \
--no-per-base \
--fast-mode \
--mapq 20 \
--threads 4 \
"${MOSDIR}/mosdepth" \
"${SAMPLE}_mapped${REFNAME}_sorted.bam"

mkdir -p clair3_output

OUTPUT_DIR="${LOCAL_WORK_DIR}/clair3_output"

singularity exec \
--bind "${LOCAL_WORK_DIR}" \
"${PATH2}"/clair3_latest.sif \
/opt/bin/run_clair3.sh \
--bam_fn="${LOCAL_WORK_DIR}/${SAMPLE}_mapped${REFNAME}_sorted.bam" \
--ref_fn="${LOCAL_WORK_DIR}/${REFERENCE}" \
--threads="${THREADS}" \
--platform="ont" \
--model_path=/opt/models/"${MODEL}" \
--output=clair3_output/ \
--sample_name="${SAMPLE}" \
--remove_intermediate_dir

mkdir -p "${SAMPLE}_spectre"

spectre CNVCaller \
--bin-size 1000 \
--coverage "${SAMPLE}_mosdepth"/ \
--snv "${LOCAL_WORK_DIR}/clair3_output/merge_output.vcf.gz" \
--reference "${MAPPING}"/"${REFERENCE}" \
--sample_id "${SAMPLE}" \
--output-dir "${SAMPLE}_spectre"/ \
--threads "${THREADS}" 2>&1 | tee "${SAMPLE}_spectre.log"
# ---------------------------------

# ---------------------------------
# Move outputs files
# ---------------------------------
# VCF files
mkdir -p "${FOLDER}/SNV/${RUN}"

cp "${OUTPUT_DIR}/merge_output.vcf.gz"     "${FOLDER}/SNV/${RUN}/${SAMPLE}.vcf.gz"
cp "${OUTPUT_DIR}/merge_output.vcf.gz.tbi" "${FOLDER}/SNV/${RUN}/${SAMPLE}.vcf.gz.tbi"

mkdir -p "${FOLDER}/CNV/${RUN}/${SAMPLE}"

cp -r "${SAMPLE}_spectre/." "${FOLDER}/CNV/${RUN}/${SAMPLE}/"
cp "${SAMPLE}_spectre.log" "${FOLDER}/CNV/${RUN}/${SAMPLE}/"

# Coverage
mkdir -p "${FOLDER}/Coverage/${RUN}"

mv "CoverageTable_${SAMPLE}_to_${REFERENCE}.tsv" "${FOLDER}/Coverage/${RUN}"
cp -r "${MOSDIR}" "${FOLDER}/Coverage/${RUN}"
# ---------------------------------

# ---------------------------------
# Cleaning Working Directory
# ---------------------------------
rm -rf "${LOCAL_WORK_DIR:?}"/*
# ---------------------------------

# ---------------------------------
mv "${FOLDER}/Scripts/P02b_Clair3.${SLURM_ARRAY_TASK_ID}.${SLURM_ARRAY_JOB_ID}."* \
   "${FOLDER}/Scripts/Sbatch_log/" 2>/dev/null || true
# ---------------------------------
