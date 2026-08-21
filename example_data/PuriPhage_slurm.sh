#!/bin/bash
#SBATCH --job-name=phage_pipeline
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=72:00:00
#SBATCH --output=PuriPhage_%j.out
#SBATCH --error=PuriPhage_%j.err

# Activate environment
source ~/.bashrc
conda activate PuriPhage

PuriPhage --input example_data/barcode17 --sample-metadata example_sample_data.tsv
