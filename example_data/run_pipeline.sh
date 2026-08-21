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

INPUT_ROOT="/hpc/dla_mm/lolijslager/phage_assemblies/purity_pipeline/260723_Nanophage260723"

SAMPLE_DIRS=("${INPUT_ROOT}"/input_data/*/)

INPUT_DIR="${SAMPLE_DIRS[$SLURM_ARRAY_TASK_ID]}"
OUTPUT_DIR="${INPUT_ROOT}/output"

echo "Processing: ${INPUT_DIR}"
echo "Saving to: ${OUTPUT_DIR}"

python PuriPhage.py \
    --input "${INPUT_DIR}" \
    --output "${OUTPUT_DIR}" \
    --viralflye-hmm-db /hpc/dla_mm/lolijslager/phage_assemblies/viralFlye/viralFlyeCode/Pfam-A.hmm.gz \
    --threads 8 \
    --mode full \