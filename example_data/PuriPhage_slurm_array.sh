#!/bin/bash
#SBATCH --job-name=PuriPhage
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=3:00:00
#SBATCH --array=0-23%6
#SBATCH --output=PuriPhage_%A_%a.out
#SBATCH --error=PuriPhage_%A_%a.err

source ~/.bashrc
conda activate PuriPhage

INPUT_ROOT="/example_data"

SAMPLE_DIRS=("${INPUT_ROOT}"/*/)

INPUT_DIR="${SAMPLE_DIRS[$SLURM_ARRAY_TASK_ID]}"

echo "Processing: ${INPUT_DIR}"

PuriPhage --input "${INPUT_DIR}" --sample-metadata example_sample_data.tsv
