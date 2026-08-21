"""
make_report.py

Functions for creating pipeline reports.

Reports generated:

- run_metadata.txt
- purity_summary.tsv
- assembly_summary.tsv
- all_assembly_hits.tsv
"""

from datetime import datetime
import os
import fcntl

def create_report_paths(report_dir):
    """
    Standardise report filenames across the pipeline.
    """

    return {
        "metadata":
            os.path.join(
                report_dir,
                "run_metadata.txt",
            ),

        "purity_summary":
            os.path.join(
                report_dir,
                "purity_summary.tsv",
            ),

        "assembly_summary":
            os.path.join(
                report_dir,
                "assembly_summary.tsv",
            ),

        "all_assembly_hits":
            os.path.join(
                report_dir,
                "all_assembly_hits.tsv",
            ),

        "prophage_analysis":
            os.path.join(
                report_dir,
                "prophage_analysis.tsv",
            ),
    }

def write_run_metadata(config, output_file):
    """
    Write reproducibility information for a run.
    """

    with open(output_file, "w") as file:

        fcntl.flock(file, fcntl.LOCK_EX) #Lock file so array jobs don't simultaneously write

        file.write(
            "Phage purification pipeline\n"
        )

        file.write(
            f"Run date\t{datetime.now()}\n\n"
        )

        file.write(
            "Configuration\n"
            "-------------\n"
        )

        for key, value in config.__dict__.items():
            file.write(f"{key}\t{value}\n")

        fcntl.flock(file, fcntl.LOCK_UN)

def initialise_assembly_summary(output_file, delete_old_files):
    """
    Create assembly summary table.
    """

    if delete_old_files or not os.path.isfile(output_file):

        with open(output_file, "w") as file:

            fcntl.flock(file, fcntl.LOCK_EX) #Lock file so array jobs don't simultaneously write

            file.write(
                "sample_name\t" #TODO: Add assembly length, comparison phage length, reads, bp, coverage
                "number_of_reads_used\t"
                "number_of_bp_used\t"
                "assembly_name\t"
                "assembly_length\t"
                "assembly_coverage\t"
                "comparison_phage\t"
                "comparison_phage_length\t"
                "sequence_identity\t"
                "residues_mismatched\n"
            )

            fcntl.flock(file, fcntl.LOCK_UN)


def append_assembly_summary(
    output_file,
    sample_name,
    total_assembly_reads,
    total_assembly_bp,
    assembly_length,
    comparison_phage_length,
    assembly_coverage,
    assembly_name,
    best_result,
):
    """
    Add one assembly result.
    """

    with open(output_file, "a") as file:

        fcntl.flock(file, fcntl.LOCK_EX) #Lock file so array jobs don't simultaneously write

        file.write(
            f"{sample_name}\t"
            f"{total_assembly_reads}\t"
            f"{total_assembly_bp}\t"
            f"{assembly_name}\t"
            f"{assembly_length}\t"
            f"{assembly_coverage}\t"
            f"{best_result['hit_id']}\t" #comparison phage name
            f"{comparison_phage_length}\t"
            f"{best_result['percent_identity']:.2f}\t" #sequence identity between assembly and comparison phage
            f"{best_result['residues_mismatched']}\n" #number of mismatched residues between assembly and comparison phage
        )

        fcntl.flock(file, fcntl.LOCK_UN)


def initialise_hit_table(output_file, delete_old_files):
    """
    Create table containing all assembly hits.
    """
    if delete_old_files or not os.path.isfile(output_file):

        with open(output_file, "w") as file:

            fcntl.flock(file, fcntl.LOCK_EX) #Lock file so array jobs don't simultaneously write

            file.write(
                "assembly_name\t"
                "hit_id\t"
                "sequence_identity\t"
                "residues_mismatched\n"
            )

            fcntl.flock(file, fcntl.LOCK_UN)


def initialise_purity_summary(output_file, delete_old_files):
    """
    Create purity summary table if it doesn't already exist.
    """

    if delete_old_files or not os.path.isfile(output_file):

        with open(output_file, "w") as file:

            fcntl.flock(file, fcntl.LOCK_EX) #Lock file so array jobs don't simultaneously write

            file.write(
                "sample\t"
                "total_reads\t"
                "mapped_reads\t"
                "mapped_bp\t"
                "compared_phage\t"
                "phage_reads\t"
                "phage_bp\t"
                "%phage_reads\t"
                "%host_reads\t"
                "%unmappable_reads\t"
                "%phage_bp\t"
                "%host_bp\t"
                "%unmappable_bp\t"
                "coverage\t"
                "number_of_mobilised_prophages\n"
            )

            fcntl.flock(file, fcntl.LOCK_UN)


def append_purity_summary(
    output_file,
    sample_name,
    total_reads,
    mapped_reads,
    mapped_bp,
    phage_stats,
):
    """
    Add one sample to the summary.
    """
    percentage_reads_matching_phage = (phage_stats['phage_read_total'] / mapped_reads)*100
    percentage_reads_matching_host = (phage_stats['host_read_total'] / mapped_reads)*100
    percentage_reads_matching_NA = (phage_stats['unmappable_read_total'] / mapped_reads)*100

    percentage_bp_matching_phage = (phage_stats['phage_bp_total'] / mapped_bp)*100
    percentage_bp_matching_host = (phage_stats['host_bp_total'] / mapped_bp)*100
    percentage_bp_matching_NA = (phage_stats['unmappable_bp_total'] / mapped_bp)*100

    with open(output_file, "a") as file:

        fcntl.flock(file, fcntl.LOCK_EX) #Lock file so array jobs don't simultaneously write

        file.write(
            f"{sample_name}\t"
            f"{total_reads}\t"
            f"{mapped_reads}\t"
            f"{mapped_bp}\t"
            f"{phage_stats['suspected_phage']}\t" #expected phage, if supplied
            f"{phage_stats['phage_read_total']}\t" #reads mapping to expected phage
            f"{phage_stats['phage_bp_total']}\t" #bp matching phage
            f"{percentage_reads_matching_phage:.1f}\t" #% reads matching phage
            f"{percentage_reads_matching_host:.1f}\t" #% reads matching host
            f"{percentage_reads_matching_NA:.1f}\t" #% reads that are unmappable
            f"{percentage_bp_matching_phage:.1f}\t" #% sequenced bp phage. Can be used as stand in for masss.
            f"{percentage_bp_matching_host:.1f}\t" #% sequenced bp host
            f"{percentage_bp_matching_NA:.1f}\t" #% sequenced bp that are unmappable
            f"{phage_stats['coverage']:.2f}\t" #coverage obtained for phage.
            f"{phage_stats['mobilised_prophages']}\n" #ratio of prophage reads to non-prophage host reads
        )

        fcntl.flock(file, fcntl.LOCK_UN)

def initialise_prophage_summary(output_file, delete_old_files):
    """
    Create purity summary table if it doesn't already exist.
    """

    if delete_old_files or not os.path.isfile(output_file):

        with open(output_file, "w") as file:

            fcntl.flock(file, fcntl.LOCK_EX) #Lock file so array jobs don't simultaneously write

            file.write(
                "Sample\t"
                "Prophage\t"
                "Prophage_coverage\t"
                "NonProphageHost_coverage\t"
                "Total_host_coverage\t"
                "Ratio_Prophage_to_NonProphageHost\t"
                "Prophage_mobilised\n"
            )

            fcntl.flock(file, fcntl.LOCK_UN)


def append_prophage_summary(
    output_file,
    sample_name,
    prophage_name,
    prophage_coverage,
    coverage_non_prophage_host,
    host_coverage,
    prophage_ratio,
    prophage_mobilised
):
    """
    Add one prophage to the summary.
    """
    with open(output_file, "a") as file:

        fcntl.flock(file, fcntl.LOCK_EX) #Lock file so array jobs don't simultaneously write

        file.write(
                f"{sample_name}\t"
                f"{prophage_name}\t"
                f"{prophage_coverage:.2f}\t"
                f"{coverage_non_prophage_host:.2f}\t"
                f"{host_coverage:.2f}\t"
                f"{prophage_ratio:.2f}\t"
                f"{prophage_mobilised}\n"
            )

        fcntl.flock(file, fcntl.LOCK_UN)
