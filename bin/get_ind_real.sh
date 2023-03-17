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

#Variables
REF_GENOME_LOC=./data/ref_gen/
REF_GENOME_NAME=${REF_GENOME_LOC}reg_gen.chr38.fa

IND_LOC=/groups/dog/data/canFam3/NGS/Epilepsy/DNA/WGS/BAM/BERN-EPIL-2015_T_12559RE_CACO.bam
IND_NAME=./data/ind/12559
IND_S=${IND_NAME}.chr38
REGION=38

# Filter out the region of interest and format to BAM
samtools view -T ${REF_GENOME_NAME} -bo ${IND_S}.bam ${IND_LOC} ${REGION}
samtools index ${IND_S}.bam

# Get the genotype likelihood based on the panel for the validation file
PANEL_NAME=/groups/dog/llenezet/test-datasets/data/panel/DVDBC.chr38
VCF=${PANEL_NAME}.sites.vcf.gz
TSV=${PANEL_NAME}.sites.tsv.gz

bcftools mpileup -f ${REF_GENOME_NAME} -I -E -a 'FORMAT/DP' -T ${VCF} -r ${REGION} ${IND_S}.bam -Ou |
bcftools call -Aim -C alleles -T ${TSV} -Ob -o ${IND_S}.bcf
bcftools index -f ${IND_S}.bcf