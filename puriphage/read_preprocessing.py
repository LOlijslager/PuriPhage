"""
read_preprocessing.py

Preprocess Oxford Nanopore reads for:

- Flye assembly
- BLAST read mapping

Processing steps:

1. Length filtering
2. Barcode/chimera filtering
3. Assembly quality filtering
4. Optional assembly downsampling
5. Mapping quality filtering
6. FASTQ and FASTA generation
"""

from pathlib import Path
import gzip
import statistics
import os

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import logging

logger = logging.getLogger("phage_pipeline")    


BARCODE_SEQUENCES = {
    "barcode01": Seq("AAGAAAGTTGTCGGTGTCTTTGTG"),
    "barcode02": Seq("TCGATTCCGTTTGTAGTCGTCTGT"),
    "barcode03": Seq("GAGTCTTGTGTCCCAGTTACCAGG"),
    "barcode04": Seq("TTCGGATTCTATCGTGTTTCCCTA"),
    "barcode05": Seq("CTTGTCCAGGGTTTGTGTAACCTT"),
    "barcode06": Seq("TTCTCGCAAAGGCAGAAAGTAGTC"),
    "barcode07": Seq("GTGTTACCGTGGGAATGAATCCTT"),
    "barcode08": Seq("TTCAGGGAACAAACCAAGTTACGT"),
    "barcode09": Seq("AACTAGGCACAGCGAGTCTTGGTT"),
    "barcode10": Seq("AAGCGTTGAAACCTTTGTCCTCTC"),
    "barcode11": Seq("GTTTCATCTATCGGAGGGAATGGA"),
    "barcode12": Seq("GTTGAGTTACAAAGCACCGATCAG"),
    "barcode13": Seq("AGAACGACTTCCATACTCGTGTGA"),
    "barcode14": Seq("AACGAGTCTCTTGGGACCCATAGA"),
    "barcode15": Seq("AGGTCTACCTCGCTAACACCACTG"),
    "barcode16": Seq("CGTCAACTGACAGTGGTTCGTACT"),
    "barcode17": Seq("ACCCTCCAGGAAAGTACCTCTGAT"),
    "barcode18": Seq("CCAAACCCAACAACCTAGATAGGC"),
    "barcode19": Seq("GTTCCTCGTGCAGTGTCAAGAGAT"),
    "barcode20": Seq("TTGCGTCCTGTTACGAGAACTCAT"),
    "barcode21": Seq("GAGCCTCTCATTGTCCGTTCTCTA"),
    "barcode22": Seq("ACCACTGCCATGTATCAAAGTACG"),
    "barcode23": Seq("CTTACTACCCAGTGAACCTCCTCG"),
    "barcode24": Seq("GCATAGTTCTGCATGATGGGTTAG"),
}

ALL_BARCODES = tuple(str(seq) for seq in BARCODE_SEQUENCES.values())


def find_fastq_files(input_path):
    """
    Accept:

    - directory of fastq.gz files
    - directory of fastq files
    - single fastq.gz file
    - single fastq file
    """

    input_path = Path(input_path)

    if input_path.is_dir():

        fastq_files = sorted(input_path.glob("*.fastq.gz"))

        if not fastq_files:
            fastq_files = sorted(input_path.glob("*.fastq"))

        return fastq_files

    if input_path.is_file():

        if input_path.name.endswith(".fastq.gz"):
            return [input_path]

        if input_path.suffix == ".fastq":
            return [input_path]

    raise ValueError(f"Unsupported input: {input_path}")


def open_fastq_file(filepath):
    """Open compressed or uncompressed FASTQ files."""

    if str(filepath).endswith(".gz"):
        return gzip.open(filepath, "rt")

    return open(filepath, "r")


def preprocess_reads(
    input_path,
    output_directory,
    sample_name,
    min_read_length,
    min_qscore,
    min_percent_identity,
    enable_barcode_filtering=True,
    enable_downsampling=False,
    target_bases=30_000_000,
    terminal_trim_bp=100,
):
    """
    Generate:

    - assembly FASTQ
    - mapping FASTA files

    Returns:
        fasta_directory
        assembly_fastq
        total_reads
        mapping_reads
    """

    total_reads = 0

    assembly_reads = 0
    mapping_reads = 0

    total_assembly_bases = 0
    total_mapping_bases = 0

    filtered_by_length = 0
    filtered_by_quality = 0
    filtered_as_chimera = 0
    filtered_as_palindrome = 0

    output_directory = Path(output_directory)
    output_directory.mkdir(exist_ok=True)

    fasta_directory = output_directory / sample_name
    fasta_directory.mkdir(exist_ok=True)

    assembly_fastq = output_directory / f"{sample_name}.fastq"

    fastq_files = find_fastq_files(input_path)
    
    assembly_handle = None
    try:
        for fastq_file in fastq_files:

            fasta_file = (
                fasta_directory
                / f"{fastq_file.stem.replace('.fastq', '')}.fasta"
            )
            
            fasta_handle = None

            try:

                with open_fastq_file(fastq_file) as input_handle:

                    for record in SeqIO.parse(
                        input_handle,
                        "fastq"
                    ):
                        #
                        # Stop once both outputs have reached
                        # their target sizes.
                        #                        

                        assembly_limit_reached = (
                            enable_downsampling
                            and total_assembly_bases >= target_bases
                        )

                        mapping_limit_reached = (
                            enable_downsampling
                            and total_mapping_bases >= target_bases
                        )

                        if (
                            assembly_limit_reached
                            and mapping_limit_reached
                        ):
                            break

                        total_reads += 1

                        sequence = record.seq
                        quality_scores = record.letter_annotations["phred_quality"]

                        #
                        # Length filter
                        #
                        # Require enough sequence to remain
                        # after mapping-end trimming.
                        #

                        min_required_length = (
                            min_read_length
                            + (2 * terminal_trim_bp)
                        )

                        if len(sequence) <= min_required_length:
                            filtered_by_length += 1
                            continue

                        #
                        # Chimera detection
                        #
                        # Ignore barcode sequence near read ends,
                        # as genuine barcode remnants may remain.
                        # Internal barcodes are treated as likely
                        # PCR chimeras.
                        #

                        inspection_sequence = str(
                            sequence[
                                terminal_trim_bp:
                                -terminal_trim_bp
                            ]
                        )

                        if enable_barcode_filtering:

                            contains_barcode = any(
                                barcode in inspection_sequence
                                for barcode in ALL_BARCODES
                            )

                            if contains_barcode:
                                filtered_as_chimera += 1
                                continue
                            
                        #
                        # Palindrome detection
                        #
                        # Filter out (likely PCR) artifacts in
                        # the form of perfectly palindromic reads.

                        if detect_palindrome(inspection_sequence): 
                            filtered_as_palindrome += 1
                            continue

                        #
                        # =================================
                        # Assembly branch
                        # =================================
                        #
                        # Flye benefits from additional reads,
                        # therefore assembly filtering is more
                        # permissive than mapping filtering.
                        #

                        median_qscore = statistics.median(
                            quality_scores
                        )

                        if (
                            median_qscore >= min_qscore
                            and not assembly_limit_reached
                        ):

                            if assembly_handle is None:
                                assembly_handle = open(
                                    assembly_fastq,
                                    "w"
                                )

                            SeqIO.write(
                                record,
                                assembly_handle,
                                "fastq"
                            )

                            assembly_reads += 1
                            total_assembly_bases += len(
                                sequence
                            )

                        #
                        # =================================
                        # Mapping branch
                        # =================================
                        #
                        # Trim low-quality read ends before
                        # evaluating mapping quality.
                        #

                        if terminal_trim_bp > 0:
                            trimmed_sequence = sequence[
                                terminal_trim_bp:
                                -terminal_trim_bp
                            ]

                            trimmed_qualities = (
                                quality_scores[
                                    terminal_trim_bp:
                                    -terminal_trim_bp
                                ]
                            )
                        else:
                            
                            trimmed_sequence = sequence
                            trimmed_qualities = quality_scores


                        high_quality_bases = sum(
                            q >= min_qscore
                            for q in trimmed_qualities
                        )

                        percent_high_quality = (
                            100
                            * high_quality_bases
                            / len(trimmed_qualities)
                        )
                        

                        if (
                            percent_high_quality >= min_percent_identity
                            and not mapping_limit_reached
                        ):
                            if fasta_handle is None:
                                fasta_handle = open(
                                    fasta_file,
                                    "w"
                                )

                            fasta_record = SeqRecord(
                                trimmed_sequence,
                                id=record.id,
                                description=""
                            )

                            SeqIO.write(
                                fasta_record,
                                fasta_handle,
                                "fasta"
                            )

                            mapping_reads += 1
                            total_mapping_bases += len(trimmed_sequence)
                            
                        elif (
                            percent_high_quality < min_percent_identity
                            and not mapping_limit_reached
                        ):
                            
                            filtered_by_quality += 1
            except ValueError:
                logger.error("Corrupted FASTQ file detected.")
                
            finally:

                if fasta_handle is not None:
                    fasta_handle.close()
                
            #
            # If both targets have been reached,
            # stop processing additional FASTQ files.
            #

            assembly_limit_reached = (
                enable_downsampling
                and total_assembly_bases >= target_bases
            )

            mapping_limit_reached = (
                enable_downsampling
                and total_mapping_bases >= target_bases
            )

            if (
                assembly_limit_reached
                and mapping_limit_reached
            ):
                break
    
    finally:

        if assembly_handle is not None:
            assembly_handle.close()

    logger.info(
        "Read preprocessing summary"
    )

    logger.info(
        "Total reads: %s",
        f"{total_reads:,}",
    )

    logger.info(
        "Assembly reads: %s",
        f"{assembly_reads:,}",
    )

    logger.info(
        "Mapping reads: %s",
        f"{mapping_reads:,}",
    )

    logger.info(
        "Assembly bases: %s",
        f"{total_assembly_bases:,}",
    )

    logger.info(
        "Mapping bases: %s",
        f"{total_mapping_bases:,}",
    )

    logger.info(
        "Rejected (length): %s",
        f"{filtered_by_length:,}",
    )

    logger.info(
        "Rejected (chimera): %s",
        f"{filtered_as_chimera:,}",
    )

    logger.info(
        "Rejected (palindrome): %s",
        f"{filtered_as_palindrome:,}",
    )

    logger.info(
        "Rejected (quality): %s",
        f"{filtered_by_quality:,}",
    )


    return (
        str(fasta_directory),
        str(assembly_fastq),
        total_reads,
        assembly_reads,
        total_assembly_bases,
        mapping_reads,
        total_mapping_bases
    )


from Bio.Seq import Seq
from Bio import Align


# Create the aligner once, outside the function
aligner = Align.PairwiseAligner(
    mode="global",
    match_score=1,
    mismatch_score=-1,
    open_gap_score=-2,
    extend_gap_score=-0.5,
)


def detect_palindrome(seq,
                      k=15,
                      min_kmer_frac=0.2,
                      min_identity=0.85):
    """
    Detect foldback palindrome reads of the form:

        A -> B -> C -> C' -> B' -> A'

    Parameters
    ----------
    record : SeqRecord
        Biopython FASTQ/FASTA record.

    k : int
        Kmer size for rapid screening.

    min_kmer_frac : float
        Fraction of kmers that must overlap before alignment is attempted.

    min_identity : float
        Alignment identity threshold.

    Returns
    -------
    bool
        True if likely palindrome, False otherwise.
    """

    mid = len(seq) // 2

    left = seq[:mid]
    right_rc = str(Seq(seq[mid:]).reverse_complement())

    # ---------
    # FAST KMER SCREEN
    # ---------

    def kmers(s):
        return {
            s[i:i+k]
            for i in range(len(s) - k + 1)
        }

    k1 = kmers(left)
    k2 = kmers(right_rc)

    if not k1 or not k2:
        return False

    kmer_similarity = len(k1 & k2) / min(len(k1), len(k2))

    if kmer_similarity < min_kmer_frac:
        return False

    # ---------
    # ALIGNMENT
    # ---------

    score = aligner.score(left, right_rc)

    max_score = min(len(left), len(right_rc))

    identity = score / max_score

    return identity >= min_identity

def filter_file_by_read_id(
    outdir,
    sample_name,
    raw_data_path,
    unmappable_reads,
    file_format
):
    """
    Export unmappable reads to a FASTQ or FASTA file.
    """

    fastq_files = find_fastq_files(raw_data_path)

    output_file = os.path.join(
        outdir,
        f"{sample_name}_unmappable_reads.{file_format}"
    )

    with open(output_file, "w") as outfile:

        for fastq_file in fastq_files:
            with open_fastq_file(fastq_file) as input_handle:
                for record in SeqIO.parse(
                    input_handle,
                    "fastq"
                ):

                    if record.id in unmappable_reads:
                        SeqIO.write(
                            record,
                            outfile,
                            file_format,
                        )
