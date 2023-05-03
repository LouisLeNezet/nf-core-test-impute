#!/bin/bash

#SBATCH --job-name=simul
#SBATCH --chdir=/groups/dog/llenezet/imputation/script/test_quality/wf_test
#SBATCH --ntasks=1
#SBATCH --mem=50G
#SBATCH --cpus-per-task=8
#SBATCH --constraint=avx2
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-type=fail         # send email if job fails
#SBATCH --mail-user=louislenezet@gmail.com

source /local/miniconda3/etc/profile.d/conda.sh
conda activate env_nf


nextflow \
    run main.nf \
    -c nextflow.config \
    --input /groups/dog/llenezet/imputation/script/test_quality/wf_test/assets/samplesheet.csv \
    --region /groups/dog/llenezet/imputation/script/test_quality/wf_test/assets/regionsheet.csv \
    --depth /groups/dog/llenezet/imputation/script/test_quality/wf_test/assets/depthsheet.csv \
    --panel /groups/dog/llenezet/imputation/script/test_quality/wf_test/assets/panelsheet.csv \
    --outdir /scratch/llenezet/nf/data/simulation \
    -work-dir /scratch/llenezet/nf/work \
    --max-cpus 8 \
    --max-memory '50.GB' \
    -profile singularity \
    -resume
