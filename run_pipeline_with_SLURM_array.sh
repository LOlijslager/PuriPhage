#!/bin/bash
#SBATCH --job-name=phage_pipeline
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=3:00:00
#SBATCH --array=0-23%6
#SBATCH --output=phage_pipeline_%A_%a.out
#SBATCH --error=phage_pipeline_%A_%a.err

source ~/.bashrc
conda activate purity_pipeline

INPUT_ROOT="/path/to/input_data"

SAMPLE_DIRS=("${INPUT_ROOT}"/*/)

INPUT_DIR="${SAMPLE_DIRS[$SLURM_ARRAY_TASK_ID]}"

echo "Processing: ${INPUT_DIR}"

python PuriPhage.py \
    --input "${INPUT_DIR}" \
    --sample-metadata sample_data.tsv \
    --viralflye-hmm-db Pfam-A.hmm.gz \
    --threads 8 \
    --mode full \