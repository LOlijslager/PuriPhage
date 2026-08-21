"""
make_viralFlye_assembly.py

Wrapper functions for:

- Flye assembly
- viralFlye phage detection
- Extraction of candidate phage assemblies

The resulting assemblies are renamed and exported
in a consistent format for downstream analysis.
"""

import os
import subprocess

from Bio import SeqIO


import logging

logger = logging.getLogger("phage_pipeline")

def make_viralFlye_assemblies(
    sample_name,
    assembly_dir,
    assembly_fastq,
    threads,
    viralflye_hmm_db,
    viralflye_completeness,
    troubleshooting
    ):

    """
    Make assemblies from Nanopore data. metaFlye is
    an adapted setting from Flye that deals better
    both short sequences (which phages are) and with
    unevenly distributed data, as is produced
    by the PCR-based ONT necessary to sequence most
    phages.

    After making assemblies with metaFlye, allow
    viralFlye to make some phage-specific corrections.
    Although I'm not sure how important this is for
    most phage assemblies, viralFlye is also pretty
    good at recognising which assemblies are phages.
    It is therefore primarily used as a filtering step.
    The assemblies and the coverage calculated by flye
    are returned.
    """     

    
    #
    # Flye assembly
    #

    flye_dir = run_metaflye(
        sample_name,
        assembly_dir,
        assembly_fastq,
        threads,
        troubleshooting,
    )

    #
    # Get assembly coverage
    #
    flye_summary_file = os.path.join(flye_dir,"assembly_info.txt")
    assembly_meta_dict = parse_flye_metadata(flye_summary_file)

    #
    # viralFlye filtering and corrections
    #

    assemblies = get_assemblies_of_interest(
        sample_name,
        flye_dir,
        assembly_meta_dict,
        assembly_fastq,
        assembly_dir,
        viralflye_hmm_db,
        threads,
        viralflye_completeness,
        troubleshooting
    )


    
    return assemblies


def run_metaflye(
    sample_name,
    assembly_dir,
    input_fastq,
    threads,
    troubleshooting=False,
):
    """
    Run Flye assembly.

    Existing assemblies are reused if
    assembly.fasta already exists.
    """
    logger.info(
        "Running metaFlye assembly for %s",
        sample_name,
    )

    output_dir = os.path.join(assembly_dir, sample_name)
    os.makedirs(output_dir, exist_ok=True)

    assembly_file = os.path.join(output_dir, "assembly.fasta")

    if not os.path.isfile(assembly_file):

        command = (
            f"flye "
            f"--nano-hq {input_fastq} "
            f"--out-dir {output_dir} "
            f"--meta "
            f"--threads {threads}"
        )

        logger.debug(command)

        subprocess.run(
            command,
            shell=True,
            check=True
        )
    else:
        logger.info(
            "Using existing Flye assembly for %s",
            sample_name,
        )


    return output_dir


def get_assemblies_of_interest(
    sample_name,
    flye_output_dir,
    assembly_meta_dict,
    input_reads,
    phage_assembly_dir,
    viralflye_hmm_db,
    threads,
    viralflye_completeness=0.01,
    troubleshooting=False,
):
    """
    Run viralFlye and extract candidate phage assemblies.

    A low completeness threshold is used by default
    because viralFlye may otherwise miss genuine
    phage sequence when assemblies are incomplete.
    """
    
    logger.info(
        "Running viralFlye for %s",
        sample_name,
    )

    viralflye_completed = False

    for filename in os.listdir(flye_output_dir):

        if not filename.endswith("viralFlye.fasta"):
            continue

        filepath = os.path.join(
            flye_output_dir,
            filename,
        )

        if any(SeqIO.parse(filepath, "fasta")):
            viralflye_completed = True
            break

    if not viralflye_completed:

        command = (
            f"viralFlye.py "
            f"--dir {flye_output_dir} "
            f"--hmm {viralflye_hmm_db} "
            f"--reads {input_reads} "
            f"--outdir {flye_output_dir} "
            f"--threads {threads} "
            f"--completeness {viralflye_completeness}"
        )

        

        subprocess.run(
            command,
            shell=True,
            check=True
        )
    else:
        logger.info(
            "Using existing viralFlye assembly for %s",
            sample_name,
        )

    assemblies_of_interest = []

    for filename in os.listdir(flye_output_dir):

        if not filename.endswith("viralFlye.fasta"):
            continue

        filepath = os.path.join(flye_output_dir, filename)

        #
        # Determine assembly type from viralFlye output
        #

        if "circular" in filename:
            assembly_type = "circ"

        elif "linear" in filename:
            assembly_type = "lin"

        elif "components" in filename:
            assembly_type = "undetermined"

        else:
            assembly_type = "unknown"

        #
        # Export each assembly individually
        #

        for i, record in enumerate(
            SeqIO.parse(filepath, "fasta"),
            start=1
        ):

            assembly_name = (
                f"{sample_name}_"
                f"{assembly_type}"
                f"assembly_{i}"
            )

            assembly_coverage = assembly_meta_dict[record.id][0]
            assembly_length = assembly_meta_dict[record.id][1]

            record.id = assembly_name
            record.description = ""

            output_file = os.path.join(
                phage_assembly_dir,
                f"{assembly_name}.fasta"
            )

            SeqIO.write(
                record,
                output_file,
                "fasta"
            )

            assemblies_of_interest.append(
                {
                    "assembly_name": assembly_name,
                    "seq_len": assembly_length,
                    "file_location": output_file,
                    "coverage": assembly_coverage
                }
            )

    return assemblies_of_interest


def parse_flye_metadata(flye_summary_file):
    """
    Parse assembly_info.txt to obtain the coverages
    of the produced assemblies.

    Outputs a dictionary with contig names as keys
    and assembly coverage as values.
    """
    
    coverage_dict = {}
    
    with open(flye_summary_file, "r") as f:
        # Read and parse header
        header = next(f).strip().lstrip("#").split()

        seq_name_idx = header.index("seq_name")
        coverage_idx = header.index("cov.")
        length_idx = header.index("length")

        for line in f:
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            fields = line.split()

            seq_name = fields[seq_name_idx]
            coverage = int(fields[coverage_idx])
            length = int(fields[length_idx])

            coverage_dict[seq_name] = [coverage,length]

    return coverage_dict
