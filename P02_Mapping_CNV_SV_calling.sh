#!/usr/bin/env bash
#SBATCH -J P02_MappingCNVSV                            # Job name 
#SBATCH -o P02_MappingCNVSV.%a.%A.out                  # File to which stdout will be written
#SBATCH -e P02_MappingCNVSV.%a.%A.err                  # File to which stderr will be written
#SBATCH -p tcourt                                      # Partition
#SBATCH -n 1                                           # Number of cores/cpus
#SBATCH -c 24                                          # Number of cpus ask per task
#SBATCH -t 01-00:00                                    # Runtime in DD-HH:MM
#SBATCH --mem 50G                                      # Memory for all cores in Gbytes (or --mem 752000 in Mbytes)
#SBATCH --mail-type=ALL                                # BEGIN,END,ALL,NONE
#SBATCH --mail-user=mymail                             # Email address

# ---------------------------------
# Version v1.0
# ---------------------------------
# author: Nicolas Nesi
# University Caen Normandy, DYNAMICURE INSERM UMR 1311
# Date: 19/10/2026
# ---------------------------------

# ---------------------------------
# Variables and Paths
# ---------------------------------
# safety for failing steps
set -euo pipefail

# conda switch environments
conda_switch() {
    set +u
    conda deactivate 2>/dev/null || true
    conda activate "$1"
    set -u
}

# Positional parameters
GROUP="ViroCaen"
PROJECT="HumanGenetic"
RUN=""
QSCORE="15"
REFERENCE="GRCh38_no_alt.fa"
CUTOFFVCF="5"

usage() {
   local exit_code="${1:-1}"
   {
    echo "Usage: sbatch --array=1-x $0 -r RUN -cv CUTOFFVCF [-q QSCORE]"
	echo "   -r, --run    Mandatory. Name of run directory"
	echo "   -q, --qscore Optional. Qscore threshold for reads quality filter [default 15]"
	echo "   -cv, --cutvcf Optional. DV minimal (reads variants) [default 5]"
	echo "   -h, --help   Show this help"
    } >&2
	exit "$exit_code"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
    -r|--run) RUN="$2"; shift 2 ;;
	-q|--qscore) QSCORE="$2"; shift 2 ;;
	-cv|--vcf) CUTOFFVCF="$2"; shift 2 ;;
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
PATH1="$FOLDER/RawReads_Nanopore/$RUN"
PATH2="$FOLDER/CleanedReads_Nanopore/$RUN"
PATH3="$FOLDER/QC_Stats/$RUN"

# Database
MAPPING="/dlocal/home/2019013/Databases/ReferenceGenomes_Mapping"

# Array job
SAMPLE="$(sed -n -e "${SLURM_ARRAY_TASK_ID}p" ${FOLDER}/Scripts/List_samples_$RUN.txt)"

# Variables
RUNDATE="$(date '+%Y-%m-%d')"

# Threads
export THREADS="${SLURM_CPUS_PER_TASK}"
echo "number of threads used:$THREADS"
# ---------------------------------

# ---------------------------------
# Copy input data to LOCAL_WORK_DIR
# ---------------------------------
cp "${PATH1}"/"${SAMPLE}".fastq.gz "${LOCAL_WORK_DIR}"
cd "${LOCAL_WORK_DIR}"/
gunzip *.fastq.gz

echo "Analysis run the: ${RUNDATE}"
echo "Working directory: ${PWD}"
echo "Variables used:"
echo "group:$GROUP project:$PROJECT run:$RUN Qscore:$QSCORE"
echo "reference for mapping:$REFERENCE cut-off VCF:$CUTOFFVCF"
# ---------------------------------

# ---------------------------------
# Environments
# ---------------------------------
module load py_env/miniconda3/25.7.0
set +u
eval "$(conda shell.bash hook)"
conda activate clean-ONT_env
set -u
# include:
# fastp=1.4.0 
# nanoq=0.10.0 
# ---------------------------------

# ---------------------------------
# Adapter trimming
# ---------------------------------
echo "Starting Nanopore adaptor trimming for sample: $SAMPLE"

porechop_abi -abi \
-i "${SAMPLE}".fastq \
-t "${THREADS}" \
-v 0 \
-o "${SAMPLE}_TrimAdapter.fastq.gz" 2>&1 | tee -a "${SAMPLE}"_porechop-abi_run"${RUNDATE}".log
# ---------------------------------

# ---------------------------------
# Size and quality trimming
# ---------------------------------
echo "Starting quality and size trimming on $SAMPLE"

fastp \
--in1 "${SAMPLE}_TrimAdapter.fastq.gz" \
--disable_adapter_trimming \
--length_required 500 \
--average_qual "${QSCORE}" \
--low_complexity_filter \
--thread "${THREADS}" \
--html "${SAMPLE}"_ReportFastp.html \
--json "${SAMPLE}"_ReportFastp.json \
--report_title "$SAMPLE fastq report" \
--out1 "${SAMPLE}_TrimAdapter_CleanedReads.fastq" 2>&1 | tee "${SAMPLE}_fastp_run${RUNDATE}.log"
# ---------------------------------

# ---------------------------------
# Statistic reads
# ---------------------------------
echo "Starting statistic of Nanopore raw fastq for sample: $SAMPLE"

nanoq \
-i "${SAMPLE}.fastq" \
--stats \
--header \
-vvv > "${SAMPLE}"_RawReads_Report.txt

gawk -i inplace 'BEGINFILE {split(FILENAME,f,".")} NR==1{print f[1];print;next} {print}' "${SAMPLE}_RawReads_Report.txt"
sed -i '2,5d' "${SAMPLE}_RawReads_Report.txt"
sed -i 's/_RawReads_Report//' "${SAMPLE}_RawReads_Report.txt"
sed -i 's/ \+ /\t/g' "${SAMPLE}_RawReads_Report.txt"

echo "Starting statistic of cleaned fastq with repeat for $SAMPLE"

nanoq \
-i "${SAMPLE}_TrimAdapter_CleanedReads.fastq" \
--stats \
--header \
-vvv > "${SAMPLE}_TrimAdapter_CleanedReads_Report.txt"

gawk -i inplace 'BEGINFILE {split(FILENAME,f,".")} NR==1{print f[1];print;next} {print}' "${SAMPLE}_TrimAdapter_CleanedReads_Report.txt"
sed -i '2,5d' "${SAMPLE}_TrimAdapter_CleanedReads_Report.txt"
sed -i "s/^$SAMPLE.*/$SAMPLE/" "${SAMPLE}_TrimAdapter_CleanedReads_Report.txt"
sed -i 's/ \+ /\t/g' "${SAMPLE}_TrimAdapter_CleanedReads_Report.txt"
# ---------------------------------

# ---------------------------------
# Environments
# ---------------------------------
conda_switch sniffles_env
# include:
# bcftools=1.24
# cutesv=2.1.4
# minimap2=v2.31
# samtools=1.24
# sniffles=2.8.1
# survivor=1.0.7
# ---------------------------------

# ---------------------------------
# Mapping ONT reads with Minimap2
# ---------------------------------
echo "Start mapping cleaned reads of sample $SAMPLE on the reference $REFERENCE"

REFNAME="${REFERENCE%.fa}"

# Mapping ONT reads
minimap2 -ax map-ont \
"${MAPPING}"/"${REFERENCE}.mmi" \
--MD \
-Y \
-t "${THREADS}" \
"${SAMPLE}_TrimAdapter_CleanedReads.fastq" \
2> "${SAMPLE}_minimap2.log" \
 | samtools sort --threads "${THREADS}" -o "${SAMPLE}_mapped${REFNAME}_sorted.bam" -
samtools index "${SAMPLE}_mapped${REFNAME}_sorted.bam"

# Statistics mapping
samtools flagstat -O tsv "${SAMPLE}_mapped${REFNAME}_sorted.bam" > "MappedReads_${SAMPLE}_to_${REFNAME}.tsv"
samtools coverage -m -A -q 7 "${SAMPLE}_mapped${REFNAME}_sorted.bam" > "CoverageHist_${SAMPLE}_to_${REFNAME}.txt"
# ---------------------------------

# ---------------------------------
# Calling for Single Variants
# ---------------------------------
echo "Single Variant calling with Sniffles 2"

sniffles \
--input "${SAMPLE}_mapped${REFNAME}_sorted.bam" \
--reference "${MAPPING}"/"${REFERENCE}" \
--tandem-repeats "${MAPPING}"/human_GRCh38_no_alt_analysis_set.trf.bed \
--vcf "${SAMPLE}_sniffles.vcf" \
--output-rnames \
--threads "${THREADS}"

echo "Structural Variant calling with cuteSV"

mkdir -p tmp_cutesv

cuteSV \
"${SAMPLE}_mapped${REFNAME}_sorted.bam" \
"${MAPPING}"/"${REFERENCE}" \
"${SAMPLE}_cutesv.vcf" \
tmp_cutesv \
--genotype \
--threads "${THREADS}" \
--min_support 3 \
--min_mapq 20 \
--report_readid \
--max_cluster_bias_INS 100 \
--diff_ratio_merging_INS 0.3

echo "Merge Sniffles and cuteSV VCFs"

ls "${SAMPLE}_sniffles.vcf" "${SAMPLE}_cutesv.vcf" > "${SAMPLE}_vcf_list.txt"

SURVIVOR merge "${SAMPLE}_vcf_list.txt" 1000 2 1 1 0 50 "${SAMPLE}_vcf_merged.vcf"
# ---------------------------------	

# ---------------------------------
# Filter VCF
# ---------------------------------
echo "Filter VCF file for translocation"

extract_bnd() {
    local vcf="$1" out="$2"
    bcftools view -i 'SVTYPE="BND"' "$vcf" \
    | bcftools query -f '%CHROM\t%POS\t%ID\t%ALT\t%FILTER\t[%GT]\t[%DR]\t[%DV]\n' \
    | awk -F'\t' 'BEGIN{OFS="\t";
        print "chrom1","pos1","id","chrom2","pos2","alt","filter","GT","DR","DV","VAF"}
      {
        c2=""; p2="";
        tmp=$4; gsub(/[][]/," ",tmp); n=split(tmp,a," ");
        for(i=1;i<=n;i++) if(a[i] ~ /:/){split(a[i],p,":"); c2=p[1]; p2=p[2]}
        if(c2!="" && $1!=c2){
          vaf=($7+$8>0)? sprintf("%.2f",$8/($7+$8)) : "NA";
          print $1,$2,$3,c2,p2,$4,$5,$6,$7,$8,vaf
        }
      }' > "$out"
}

extract_bnd "${SAMPLE}_sniffles.vcf" "${SAMPLE}_sniffles_bnd_interchrom.tsv"
extract_bnd "${SAMPLE}_cutesv.vcf"   "${SAMPLE}_cutesv_bnd_interchrom.tsv"

concordance_bnd() {
    local tol="$1" cute="$2" snif="$3" out="$4"
    awk -F'\t' -v OFS='\t' -v tol="$tol" '
    function canon(c1,p1,c2,p2) {
        if (c1 > c2) { C1=c2; P1=p2; C2=c1; P2=p1 }
        else         { C1=c1; P1=p1; C2=c2; P2=p2 }
    }
    function abs(x) { return (x < 0) ? -x : x }
    FNR==NR {                                  # 1er fichier : cuteSV
        if (FNR > 1) {
            canon($1,$2,$4,$5)
            key = C1 SUBSEP C2
            i = ++n[key]
            cp1[key,i]=P1; cp2[key,i]=P2
            cdv[key,i]=$10; cvaf[key,i]=$11; cflt[key,i]=$7
        }
        next
    }
    FNR==1 {                                   # 2e fichier : Sniffles2
        print "chrom1","pos1","chrom2","pos2",
              "sniffles_DV","sniffles_VAF","sniffles_filter",
              "cutesv_pos1","cutesv_pos2","cutesv_DV","cutesv_VAF","cutesv_filter",
              "delta_pos1","delta_pos2"
        next
    }
    {
        canon($1,$2,$4,$5)
        key = C1 SUBSEP C2
        for (i = 1; i <= n[key]; i++) {
            if (abs(cp1[key,i]-P1) <= tol && abs(cp2[key,i]-P2) <= tol) {
                junction = key SUBSEP P1 SUBSEP P2
                if (!(junction in done)) {
                    done[junction] = 1
                    print C1,P1,C2,P2,$10,$11,$7,
                          cp1[key,i],cp2[key,i],cdv[key,i],cvaf[key,i],cflt[key,i],
                          cp1[key,i]-P1,cp2[key,i]-P2
                }
                break
            }
        }
    }' "$cute" "$snif" > "$out"
}

TOL=1000   # distance max (bp) entre les points de cassure des deux callers

concordance_bnd "$TOL" \
"${SAMPLE}_cutesv_bnd_interchrom.tsv" \
"${SAMPLE}_sniffles_bnd_interchrom.tsv" \
"${SAMPLE}_bnd_concordant.tsv"

echo "Translocations found by the two callers: $(( $(wc -l < "${SAMPLE}_bnd_concordant.tsv") - 1 ))"

# Translocations : DV >= seuil chez les deux callers (colonnes 5 et 10)
awk -F'\t' -v c="$CUTOFFVCF" 'NR==1 || ($5>=c && $10>=c)' \
"${SAMPLE}_bnd_concordant.tsv" > "${SAMPLE}_bnd_concordant_DV${CUTOFFVCF}.tsv"

echo "Filtering for Deletion Duplication"

extract_deldup() {
    local vcf="$1" out="$2" minlen="$3"
    bcftools view -i 'SVTYPE="DEL" || SVTYPE="DUP"' "$vcf" \
    | bcftools query -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\t%INFO/SVLEN\t%ID\t%FILTER\t[%GT]\t[%DR]\t[%DV]\n' \
    | awk -F'\t' -v OFS='\t' -v minlen="$minlen" '
        BEGIN{print "chrom","start","end","type","svlen","id","filter","GT","DR","DV","VAF"}
        $1 ~ /^chr([0-9]+|X|Y)$/ {
            len = ($5=="." || $5=="") ? $3-$2 : $5
            if (len < 0) len = -len
            if (len >= minlen) {
                vaf = ($9+$10>0) ? sprintf("%.2f",$10/($9+$10)) : "NA"
                print $1,$2,$3,$4,len,$6,$7,$8,$9,$10,vaf
            }
        }' > "$out"
}

extract_deldup "${SAMPLE}_sniffles.vcf" "${SAMPLE}_sniffles_deldup.tsv" 50
extract_deldup "${SAMPLE}_cutesv.vcf" "${SAMPLE}_cutesv_deldup.tsv" 50

# DEL/DUP : DV >= seuil et FILTER=PASS (colonnes 10 et 7)
for caller in sniffles cutesv; do
    awk -F'\t' -v c="$CUTOFFVCF" 'NR==1 || ($10>=c && $7=="PASS")' \
        "${SAMPLE}_${caller}_deldup.tsv" > "${SAMPLE}_${caller}_deldup_DV${CUTOFFVCF}.tsv"
done
# ---------------------------------

# ---------------------------------
# Move outputs files
# ---------------------------------
# Cleaned Fastq

mkdir -p "${PATH2}"

pigz -f "${SAMPLE}_TrimAdapter_CleanedReads.fastq"
mv "${SAMPLE}_TrimAdapter_CleanedReads.fastq.gz" "${PATH2}"

# Statistic
mkdir -p "${PATH3}/RawReads"

mv "${SAMPLE}_RawReads_Report.txt" "${PATH3}"/RawReads/

mkdir -p "${PATH3}/CleanedReads"

mv "${SAMPLE}_TrimAdapter_CleanedReads_Report.txt" "${PATH3}/CleanedReads"

mkdir -p "${FOLDER}/Log/${RUN}_run${RUNDATE}"

# bam Minimap2
mkdir -p "${FOLDER}/Mapping/${RUN}"

mv "${SAMPLE}_mapped${REFERENCE}_sorted.bam" "${FOLDER}/Mapping/${RUN}/"
mv "${SAMPLE}_mapped${REFERENCE}_sorted.bam.bai" "${FOLDER}/Mapping/${RUN}/"
mv "Coverage_${SAMPLE}_to_${REFERENCE}.txt" "${FOLDER}/Mapping/${RUN}/"

cp "${MAPPING}"/"${REFERENCE}" "${FOLDER}/Mapping/${RUN}/"

# Coverage
mkdir -p "${FOLDER}/Coverage/${RUN}"

mv "MappedReads_${SAMPLE}_to_${REFERENCE}.tsv" "${FOLDER}/Mapping/${RUN}/"
mv "CoverageHist_${SAMPLE}_to_${REFERENCE}.txt" "${FOLDER}/Mapping/${RUN}/"

# VCF files
mkdir -p "${FOLDER}/VCF/${RUN}"

mv "${SAMPLE}_sniffles.vcf" "${FOLDER}/VCF/${RUN}"
mv "${SAMPLE}_cutesv.vcf" "${FOLDER}/VCF/${RUN}"
mv "${SAMPLE}_vcf_list.txt" "${FOLDER}/VCF/${RUN}"
mv "${SAMPLE}_vcf_merged.vcf" "${FOLDER}/VCF/${RUN}"
mv "${SAMPLE}_cutesv_bnd_interchrom.tsv" "${FOLDER}/VCF/${RUN}"
mv "${SAMPLE}_sniffles_bnd_interchrom.tsv" "${FOLDER}/VCF/${RUN}"
mv "${SAMPLE}_bnd_concordant.tsv" "${FOLDER}/VCF/${RUN}"
mv "${SAMPLE}_sniffles_deldup.tsv" "${FOLDER}/VCF/${RUN}"
mv "${SAMPLE}_cutesv_deldup.tsv" "${FOLDER}/VCF/${RUN}"

# Log files
mv *.log "${FOLDER}/Log/${RUN}_run${RUNDATE}"
# ---------------------------------

# ---------------------------------
# Cleaning Working Directory
# ---------------------------------
rm -rf "${LOCAL_WORK_DIR:?}"/*
# ---------------------------------

# ---------------------------------
mv "${FOLDER}"/Scripts/P02_MappingCNVSV."${SLURM_ARRAY_TASK_ID}".* "${FOLDER}"/Scripts/Sbatch_log
# ---------------------------------
