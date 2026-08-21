# PuriPhage

PuriPhage is an automated pipeline to analyse phage DNA sequenced using long-read sequencing approaches. It creates assemblies, selects phage-based assemblies, and determines the purity of the identified phage. Although created to analyse ONT data, the purity part will also run on Illumina data.

## Citations

If you use the assembly part of PuriPhage, please cite:

Kolmogorov, M., Bickhart, D. M., Behsaz, B., Gurevich, A., Rayko, M., Shin, S. B., ... & Pevzner, P. A. (2020). metaFlye: scalable long-read metagenome assembly using repeat graphs. Nature methods, 17(11), 1103-1110.

Antipov, D., Rayko, M., Kolmogorov, M., & Pevzner, P. A. (2022). viralFlye: assembling viruses and identifying their hosts from long-read metagenomics data. Genome biology, 23(1), 57.

## Installation

Install PuriPhage using the following commands:

```
git clone https://github.com/LOlijslager/PuriPhage.git

# Install Puriphage
cd PuriPhage
conda env create -f PuriPhage.yml
conda activate PuriPhage
pip install -e .
```

To run the assembly part of the pipeline, you will need the viralFlye reference database, which can be downloaded via:
```
wget http://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam34.0/Pfam-A.hmm.gz
```

## Verifying installation and example usage.
```
cd example_data
conda activate PuriPhage
PuriPhage --input example_data --sample-metadata example_sample_data.tsv
conda deactivate PuriPhage
```
In example_data, there will now be a folder named "output", with subdirectories: "assemblies", "logs", and "reports". If everything is working as it should, de Reports folder will match the output in "example_reports".

## Running PuriPhage
```
# Before each use:
conda activate PuriPhage

# Basic command
PuriPhage --input input_reads

# Basic command extended
PuriPhage --input input_reads --sample-metadata sample_data.tsv --viralflye-hmm-db path/to/Pfam-A.hmm.gz --references references  --output results --mode full

# Make only the assembly
PuriPhage --input input_reads --mode assembly

# Check only the purity of the sample
PuriPhage --input input_reads --mode purity

#After each use:
conda deactivate PuriPhage
```
Some notes: 
- Usage of the sample_data.tsv is optional. If it is excluded, the code will determine which DNA element from the reference database best describes the data and uses that for comparison, and all reads matching to any Host will be considered Host DNA. If your sample is included in sample_data.tsv (matching the name of a fasta file in the reference database), it will use that one for comparison instead.
- To conserve processing time, the creation of an assembly will be skipped if a completed assembly is already found in the assembly dir from a previous run.
- In case you're working with SLURM, there are some example files in the example_dir as well.

## Interpreting output

If the complete pipeline is run, the output folder will consist of three subfolders: assemblies, logs, and reports.

### Assemblies

Contains, for each sample, an output folder and a fasta file of all identified phage assemblies. The assemblies are structures as [samplename]\_[assemblytype]\_[assemblyidentifier].fasta. In this, the assembly type indicates what type of assembly viralFlye determined the assembly to be. 

### Logs
Contains, for each sample, a log file for the run with which it was created.

### reports
Contains four output file: 
- all_assembly_hits.tsv:
Contains some basic information on sequences in the reference database to which the created assembly maps. Purpose is to determine if a contamination occurred with another phage you have in your lab if data is not as expected. 
- assembly_summary.tsv:
Contains details of all phage assemblies created during this run. If no assembly is mentioned for a particular sample, no assembly could be made that had sufficient coverage to be considered a full phage assembly, as determined by viralFlye. 
- prophage_analysis.tsv:
Contains the raw data used to determine whether a prophage is induced (defined as prophage coverage is >5 more expressed than is expected based on the non-prophage host coverage.) Result of this is included in purity_summary.tsv.
- purity_summary.tsv:
Contains all information relevant to determine whether the sample is pure.
- run_metadata.txt:
Contains the metadata, with all parameters used for this run.

## Options
```
$ PuriPhage --help
usage: PuriPhage.py [-h] --input INPUT [--output OUTPUT] [--references REFERENCES] [--sample-metadata SAMPLE_METADATA] [--blast-exe BLAST_EXE] [--makeblastdb-exe MAKEBLASTDB_EXE]
                    [--viralflye-hmm-db VIRALFLYE_HMM_DB] [--mode {full,purity,assembly}] [--min-read-length MIN_READ_LENGTH] [--min-qscore MIN_QSCORE] [--min-percent-identity MIN_PERCENT_IDENTITY]
                    [--terminal-trim-bp TERMINAL_TRIM_BP] [--disable-barcode-filtering] [--enable-downsampling] [--target-bases TARGET_BASES] [--viralflye-completeness VIRALFLYE_COMPLETENESS]
                    [--threads THREADS] [--export_unmappable_reads] [--troubleshooting] [--force]

Phage purification pipeline

options:
  -h, --help            show this help message and exit
  --input INPUT         Input FASTQ file (.fastq or .fastq.gz), a directory containing FASTQ files for a single sample, or a directory containing one subdirectory per sample.
  --output OUTPUT       Output directory (Default: output).
  --references REFERENCES
                        Reference sequence directory containing 'Phage' 'Host' and 'Prophage' subdirectories with FASTA files (default: input_sequences )
  --sample-metadata SAMPLE_METADATA
                        Optional TSV file containing sample information. Required columns: Sample, Phage, Host. Phage and Host names must match corresponding filenames in the reference database.
  --blast-exe BLAST_EXE
                        Path to blastn executable.
  --makeblastdb-exe MAKEBLASTDB_EXE
                        Path to makeblastdb executable.
  --viralflye-hmm-db VIRALFLYE_HMM_DB
                        Path to viralFlye HMM database
  --mode {full,purity,assembly}
                        Pipeline mode (default: full). In purity mode, assembly and assembly comparison are skipped and only read mapping against the reference database is performed. In Assembly mode,
                        mapping against the reference database is skipped.
  --min-read-length MIN_READ_LENGTH
                        Default: 1000
  --min-qscore MIN_QSCORE
                        Default: 20
  --min-percent-identity MIN_PERCENT_IDENTITY
                        Minimum percent identity required for a hit to be accepted (default: 98). This is defined as by how much of the read is explained by a reference sequence.
  --terminal-trim-bp TERMINAL_TRIM_BP
                        Number of bases removed from both ends of each read prior to barcode detection and read mapping. Read ends are often lower quality than the rest of the sequence. Default: 100.
  --disable-barcode-filtering
                        Disable barcode/chimera filtering. Internal barcode sequences are normally treated as evidence of chimeric reads caused by PCR artefacts or sequencing errors. Disable this option
                        if a barcode sequence is genuinely expected within the phage genome.
  --enable-downsampling
                        Downsample reads used. This can reduce assembly time and may improve assembly success for very high-coverage datasets or prevent memory issues.
  --target-bases TARGET_BASES
                        Target number of bases retained after downsampling for assembly (default: 30,000,000).
  --viralflye-completeness VIRALFLYE_COMPLETENESS
                        See the viralFlye --completeness parameter (default: 0.01).
  --threads THREADS     Default: 1
  --export_unmappable_reads
                        If reads can't be mapped, they will be written to a new fastq file.
  --troubleshooting     Enable detailed logging and display external commands used during pipeline execution.
  --force               Overwrite existing summary files. Will not delete assemblies. Not compatible with HPC array jobs.
```


