#!/bin/bash

#SBATCH --job-name=impute_test_quality
#SBATCH --chdir=/groups/dog/llenezet/imputation/script/test_quality/nf-core-testquality
#SBATCH --ntasks=1
#SBATCH --mem=40G
#SBATCH --cpus-per-task=8
#SBATCH --constraint=avx2
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-type=fail         # send email if job fails
#SBATCH --mail-user=louislenezet@gmail.com

source /local/miniconda3/etc/profile.d/conda.sh
conda activate env_nf

nohup nextflow \
    run main.nf \
    -c nextflow.config \
    --outdir ./data \
    -work-dir ./work \
    -profile singularity \
    -resume &> nohupBQSR.out