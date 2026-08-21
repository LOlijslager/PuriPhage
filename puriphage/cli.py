
"""
PuriPhage.py

Configure parsed arguments and send them to run_pipeline
"""


from dataclasses import dataclass
from datetime import datetime
import argparse

from run_pipeline import run_pipeline



@dataclass
class PipelineConfig:

    #
    # Input
    #

    input_path: str
    output_dir: str
    reference_sequences: str
    sample_metadata: str

    #
    # External tools
    #

    blast_exe: str
    makeblastdb_exe: str

    #
    # Assembly
    #

    viralflye_hmm_db: str
    viralflye_completeness: float = 0.01

    #
    # Read filtering
    #

    min_read_length: int = 1000
    min_qscore: int = 20
    min_percent_identity: float = 98

    terminal_trim_bp: int = 100

    enable_barcode_filtering: bool = True

    enable_downsampling: bool = False
    export_unmappable_reads: bool = False
    target_bases: int = 30_000_000

    #
    # Runtime
    #

    mode: str = "full"
    threads: int = 1
    troubleshooting: bool = False
    delete_old_files: bool = False

def build_config(args):

    return PipelineConfig(
        input_path=args.input,
        output_dir=args.output,
        reference_sequences=args.references,
        sample_metadata=args.sample_metadata,

        blast_exe=args.blast_exe,
        makeblastdb_exe=args.makeblastdb_exe,

        viralflye_hmm_db=args.viralflye_hmm_db,
        viralflye_completeness=args.viralflye_completeness,

        min_read_length=args.min_read_length,
        min_qscore=args.min_qscore,
        min_percent_identity=args.min_percent_identity,

        terminal_trim_bp=args.terminal_trim_bp,

        enable_barcode_filtering=(
            not args.disable_barcode_filtering
        ),

        enable_downsampling=args.enable_downsampling,
        export_unmappable_reads=args.export_unmappable_reads,
        target_bases=args.target_bases,

        mode=args.mode,
        threads=args.threads,
        troubleshooting=args.troubleshooting,
        delete_old_files=args.force,
    )

def parse_arguments():
    """
    Parse command-line arguments.
    """

    parser = argparse.ArgumentParser(
        description="Phage purification pipeline"
    )

    #
    # Required arguments
    #

    parser.add_argument(
        "--input",
        required=True,
        help=(
        "Input FASTQ file (.fastq or .fastq.gz), a directory "
        "containing FASTQ files for a single sample, or a "
        "directory containing one subdirectory per sample."
        )
    )

    parser.add_argument(
        "--output",
        default="output",
        help="Output directory (Default: output)."
    )

    parser.add_argument(
        "--references",
        default="input_sequences",
        help=(
        "Reference sequence directory containing 'Phage' "
        " 'Host' and 'Prophage' subdirectories with " 
        " FASTA files (default: input_sequences )"
        )
    )

    parser.add_argument(
    "--sample-metadata",
    default=None,
    help=(
        "Optional TSV file containing sample information. "
        "Required columns: Sample, Phage, Host. "
        "Phage and Host names must match corresponding filenames"
        " in the reference database."
        )
    )

    parser.add_argument(
        "--blast-exe",
        default="blastn",
        help="Path to blastn executable."
    )

    parser.add_argument(
        "--makeblastdb-exe",
        default="makeblastdb",
        help="Path to makeblastdb executable."
    )

    parser.add_argument(
        "--viralflye-hmm-db",
        help="Path to viralFlye HMM database"
    )

    #
    # Pipeline mode
    #

    parser.add_argument(
        "--mode",
        choices=["full", "purity", "assembly"],
        default="full",
        help=(
        "Pipeline mode (default: full). In purity mode, "
        "assembly and assembly comparison are skipped and "
        "only read mapping against the reference database "
        "is performed. In Assembly mode, mapping "
        "against the reference database is skipped."
        )
        )

    #
    # Read filtering
    #

    parser.add_argument(
        "--min-read-length",
        type=int,
        default=1000,
        help = "Default: 1000"
    )

    parser.add_argument(
        "--min-qscore",
        type=int,
        default=20,
        help = "Default: 20"
    )

    parser.add_argument(
        "--min-percent-identity",
        type=float,
        default=98,
        help=(
        "Minimum percent identity required for a hit "
        "to be accepted (default: 98). This is defined as "
        "by how much of the read is explained by a "
        "reference sequence."
        )
    )

    parser.add_argument(
        "--terminal-trim-bp",
        type=int,
        default=100,
        help=(
        "Number of bases removed from both ends of each "
        "read prior to barcode detection and read mapping. "
        "Read ends are often lower quality than the rest "
        "of the sequence. Default: 100."
        )
    )

    #
    # Barcode filtering
    #

    parser.add_argument(
        "--disable-barcode-filtering",
        action="store_true",
        help=(
        "Disable barcode/chimera filtering. Internal barcode "
        "sequences are normally treated as evidence of "
        "chimeric reads caused by PCR artefacts or sequencing "
        "errors. Disable this option if a barcode sequence "
        "is genuinely expected within the phage genome."
        )
    )

    #
    # Assembly downsampling
    #

    parser.add_argument(
        "--enable-downsampling",
        action="store_true",
        help=(
        "Downsample reads used. This can reduce assembly "
        "time and may improve assembly success for very "
        "high-coverage datasets or prevent memory issues." 
        )
    )

    parser.add_argument(
        "--target-bases",
        type=int,
        default=30_000_000,
        help=(
        "Target number of bases retained after "
        "downsampling for assembly (default: 30,000,000)."
        )
    )
        
    #
    # viralFlye
    #

    parser.add_argument(
        "--viralflye-completeness",
        type=float,
        default=0.01,
        help="See the viralFlye --completeness parameter (default: 0.01)."
    )

    #
    # Runtime
    #

    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Default: 1"
    )
    
    parser.add_argument(
        "--export_unmappable_reads",
        action="store_true",
        help=(
        "If reads can't be mapped, they will be "
        "written to a new fastq file." 
        )
    )

    parser.add_argument(
        "--troubleshooting",
        action="store_true",
        help=(
        "Enable detailed logging and display external "
        "commands used during pipeline execution."
        )
    )

    parser.add_argument(
        "--force",
        action="store_true",
        help=(
        "Overwrite existing summary files. Will not delete assemblies. "
        "Not compatible with HPC array jobs."
        )
    )

    return parser.parse_args()



def main():

    args = parse_arguments()

    config = build_config(args)

    run_pipeline(config)


if __name__ == "__main__":
    main()
