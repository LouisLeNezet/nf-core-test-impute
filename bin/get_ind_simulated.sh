#!/bin/sh

#SBATCH --job-name=phasing_panel
#SBATCH --chdir=/groups/dog/llenezet/test-datasets/
#SBATCH --ntasks=1
#SBATCH --mem=40G
#SBATCH --cpus-per-task=2
#SBATCH --constraint=avx2
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-type=fail         # send email if job fails
#SBATCH --mail-user=louislenezet@gmail.com

source /local/miniconda3/etc/profile.d/conda.sh
conda activate env_tools

REF_GENOME_LOC=./data/ref_gen/
REF_GENOME_NAME=${REF_GENOME_LOC}reg_gen.chr38.fa

IND_LOC=/groups/dog/data/canFam3/NGS/Epilepsy/DNA/WGS/BAM/BERN-EPIL-2015_T_12559RE_CACO.bam
IND_NAME=./data/ind/12559
IND_S=${IND_NAME}.chr38
IND_S_1X=${IND_S}.1x
REGION=38
MAP_SNP=./data/illumina/chr38.snp.map
PANEL_NAME=/groups/dog/llenezet/test-datasets/data/panel/DVDBC.chr38
VCF=${PANEL_NAME}.sites.vcf.gz
TSV=${PANEL_NAME}.sites.tsv.gz

# Downsampling the individual data to 1X
MEAN_DEPTH=$(samtools coverage ${IND_S}.bam -r ${REGION} | \
    awk -F'\t' '(NR==2){ print $7}')
FRAC_DEPTH=$(echo "scale=5; 1/$MEAN_DEPTH" | bc)
samtools view -T ${REF_GENOME_NAME} -s 1${FRAC_DEPTH} -bo ${IND_S_1X}.bam ${IND_S}.bam ${REGION}
samtools index ${IND_S_1X}.bam

# Compute genotype likelihood based on the panel
bcftools mpileup -f ${REF_GENOME_NAME} -I -E -a 'FORMAT/DP' -T ${VCF} -r ${REGION} ${IND_S_1X}.bam -Ou |
bcftools call -Aim -C alleles -T ${TSV} -Oz -o ${IND_S_1X}.vcf.gz
bcftools index -f ${IND_S_1X}.vcf.gz

# Get individual SNP
plink2 --chr-set 38 --bcf ${IND_S}.bcf \
    --set-missing-var-ids @:# \
    --max-alleles 2 \
    --output-chr '26' \
    --snps-only --allow-extra-chr \
    --extract ${MAP_SNP} \
    --chr 1-38,X --recode vcf bgz\
    --out ${IND_S}.snp
bcftools index -f ${IND_S}.snp.vcf.gz


bcftools view -e 'GT="./."||GT="."' ${IND_S}.snp.vcf.gz -Oz -o ${IND_S}.snp.filtered.vcf.gz
bcftools index -f ${IND_S}.snp.filtered.vcf.gz

# Phase the SNP data
SHAPEIT5_phase_common -I ${IND_S}.snp.filtered.vcf.gz -H ${PANEL_NAME}.phased.bcf -O ${IND_S}.phased.vcf -R ${REGION}