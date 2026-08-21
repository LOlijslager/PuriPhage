#!/bin/bash
#SBATCH --job-name=phage_pipeline
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=72:00:00
#SBATCH --output=phage_pipeline_%j.out
#SBATCH --error=phage_pipeline_%j.err

# Activate environment
source ~/.bashrc
conda activate purity_pipeline

python PuriPhage.py \
    --input /input_data/ \
    --sample-metadata sample_data.tsv \
    --viralflye-hmm-db Pfam-A.hmm.gz \
    --threads 8