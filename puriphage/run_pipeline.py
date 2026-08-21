"""
run_pipeline.py

Pipeline for:

1. Read preprocessing
2. Flye assembly
3. viralFlye filtering
4. BLAST identification
5. Purity assessment
6. Reporting
"""
from pathlib import Path
from tempfile import TemporaryDirectory
import os
import csv

from .read_preprocessing import preprocess_reads
from .logging_utils import setup_logging

from .make_viralFlye_assembly import (
    make_viralFlye_assemblies
)

from .run_blast import (
    create_reference_database,
    create_assembly_database,
    run_blast_best_hit,
)

from .analyse_results import (
    count_hits,
    process_expected_or_best_phage,
    filter_file_by_read_id
)

from .make_report import (
    create_report_paths,
    write_run_metadata,
    
    initialise_purity_summary,
    initialise_assembly_summary,
    initialise_hit_table,
    initialise_prophage_summary,
    
    append_purity_summary,
    append_assembly_summary,
)

def run_pipeline(config):

    validate_config(config)
    
    directories = create_directories(
        config.output_dir,
        config.export_unmappable_reads
    )

    if config.sample_metadata:
        #
        # Configure run details/folders/etc
        #
        
        sample_metadata = read_sample_metadata(
            config.sample_metadata,
        )
    else:
        sample_metadata = None

    report_files = create_report_paths(
    str(directories["reports"])
    )

    write_run_metadata(
        config,
        report_files["metadata"],
    )

    if config.mode in ("purity", "full"):
        initialise_purity_summary(
            report_files["purity_summary"],
            config.delete_old_files
        )
        initialise_prophage_summary(
            report_files["prophage_analysis"],
            config.delete_old_files
        )
        
    if config.mode in ("assembly", "full"):
        initialise_assembly_summary(
            report_files["assembly_summary"],
            config.delete_old_files
        )

        initialise_hit_table(
            report_files["all_assembly_hits"],
            config.delete_old_files
        )

    samples = discover_samples(
        config.input_path
    )

    for sample in samples:
        #
        # Start analysing sample. Find sample details first.
        #

        sample_name = sample["sample_name"]

        logger = setup_logging(
            str(directories["logs"]),
            sample_name,
            config.troubleshooting,
        )

        logger.info(
            "Processing sample: %s",
            sample_name,
        )

        if config.sample_metadata != None:
            metadata_supplied = True
            metadata = sample_metadata.get(
                sample_name,
                {}
            )
            logger.info(
                    "Sample %s found in sample metadata.",
                    sample_name,
                )

            expected_phage = metadata.get("phage")
            if expected_phage != "":
                logger.info(
                    "Comparing to phage: %s",
                    expected_phage,
                )
            else:
                expected_phage = None
            if metadata.get("host") != "": 
                supplied_host = metadata.get("host")
                logger.info(
                    "Using host: %s",
                    supplied_host,
                )
            else:
                supplied_host = None
                
        else:
            metadata_supplied = False
            expected_phage = None
            supplied_host = None
        
        with TemporaryDirectory(dir=config.output_dir) as temp_dir:

            #
            # Start configuring for BLAST
            # Save intermediate data in a temporary directory
            #

            temp_dir = Path(temp_dir)
            
            temporary_directories = {
                    "mapping": temp_dir / "mapping",
                    "blast_db": temp_dir / "blast_db",
                }
            
            for directory in temporary_directories.values():
                    directory.mkdir(
                        parents=True,
                        exist_ok=True,
                    )

            temp_mapping_dir = temporary_directories["mapping"]
            temp_blast_dir = temporary_directories["blast_db"]

            logger.info(
                "Using temporary working directory: %s",
                temp_dir,
            )

            #
            # Create reference database
            #

            (
                reference_db,
                reference_fasta,
                sequence_lengths,
                phage_list,
                host_list,
                prophage_list
            ) = create_reference_database(
                config.makeblastdb_exe,
                config.reference_sequences,
                str(temp_blast_dir),
                config.troubleshooting,
            )

            if expected_phage != None:
                if expected_phage not in phage_list:
                    raise ValueError(
                        f"Invalid Phage '{expected_phage}' in --sample-metadata "
                        f"for sample {sample_name}. "
                        "Phage names must match corresponding filenames"
                        " in the reference database."
                        )
            if supplied_host != None:
                if supplied_host not in host_list:
                    raise ValueError(
                        f"Invalid Host '{supplied_host}' in --sample-metadata. "
                        f"for sample {sample_name}. "
                        "Host names must match corresponding filenames"
                        " in the reference database."
                        )
                else:
                    host_list = [supplied_host]
                    new_prophage_list = []
                    prophage_found = False
                    for prophage in prophage_list:
                        if supplied_host+"_prophage" in prophage:
                            new_prophage_list.append(prophage)
                            prophage_found = True
                    prophage_list = new_prophage_list
                    if prophage_found:
                        logger.info(
                            "Prophages %s found for host %s",
                            (", ").join(prophage_list),
                            supplied_host,
                        )
                        
                    else:
                        logger.info(
                            "No prophages found for host %s",
                            supplied_host
                        )
                        
            #
            # Read preprocessing
            #

            (
                fasta_directory,
                assembly_fastq,
                total_reads,
                total_assembly_reads,
                total_assembly_bp,
                mapping_reads,
                mapping_bp
            ) = preprocess_reads(
                input_path=sample["input_path"],
                output_directory=str(temp_mapping_dir),
                sample_name=sample_name,
                min_read_length = config.min_read_length,
                min_qscore = config.min_qscore,
                min_percent_identity = config.min_percent_identity,
                enable_barcode_filtering = config.enable_barcode_filtering,
                enable_downsampling = config.enable_downsampling,
                target_bases = config.target_bases,
                terminal_trim_bp = config.terminal_trim_bp,
            )

            #
            #initialise the phage_assembly variable.
            #Default is expected phage, overwritten if a valid assembly is created.
            #
            
            phage_assembly = expected_phage 

            #
            # Create phage assembly
            #

            if config.mode in ("assembly", "full"):

                assemblies = make_viralFlye_assemblies(
                    sample_name,
                    str(directories["assemblies"]),
                    assembly_fastq,
                    config.threads,
                    config.viralflye_hmm_db,
                    config.viralflye_completeness,
                    config.troubleshooting
                )


                if not assemblies:

                    logger.warning(
                        "No candidate phage assemblies found for %s",
                        sample_name,
                    )

                    continue

                #
                # Make BLAST database with assemblies as wel as references
                #

                assembly_files = [
                    assembly["file_location"]
                    for assembly in assemblies
                ]
                
                assembly_db, sequence_lengths = (
                    create_assembly_database(
                        config.makeblastdb_exe,
                        reference_fasta,
                        assembly_files,
                        str(temp_blast_dir),
                        sequence_lengths,
                        config.troubleshooting,
                    )
                )

                #
                # Evaluate assemblies
                #

                for assembly in assemblies:

                    blast_results = (
                        run_blast_best_hit(
                            config.blast_exe,
                            assembly_db,
                            assembly["file_location"],
                            config.min_percent_identity,
                            True,
                            host_list,
                            config.threads,
                            config.troubleshooting,
                        )
                    )

                    matching_results = process_expected_or_best_phage(
                        expected_phage,
                        sequence_lengths,
                        assembly["assembly_name"],
                        blast_results,
                        os.path.join(
                            directories["reports"],
                            "all_assembly_hits.tsv",
                        ),
                        "\t",
                    )

                    logger.info(
                        "%s assembly mapping result: %s (%.1f%%)",
                        assembly["assembly_name"],
                        matching_results["hit_id"],
                        matching_results["percent_identity"],
                    )

                    assembly_length = assembly["seq_len"]
                    comparison_phage_length = sequence_lengths["Phage"].get(
                        matching_results["hit_id"].split("/")[0],
                        "N.A."
                    )
                    assembly_coverage = assembly["coverage"]

                    append_assembly_summary(
                        report_files["assembly_summary"],
                        sample_name,
                        total_assembly_reads,
                        total_assembly_bp,
                        assembly_length,
                        comparison_phage_length,
                        assembly_coverage,
                        assembly["assembly_name"],
                        matching_results
                        )
                #
                # Overwrite phage assembly to determine purity of with the new assembly
                #
                if assemblies != []:
                    phage_assembly = assemblies[0]["assembly_name"]
                    reference_db = assembly_db
                    
            #
            # Determine phage purity, either of created assembly or (of none is provided) of phage assembly supplied
            #

            if config.mode in ("purity", "full"):
                #
                # BLAST all FASTA files
                #

                merged_results = {}
                logger.info(
                    "Running BLAST on %s",
                    sample_name,
                )

                for fasta_file in Path(fasta_directory).glob("*.fasta"):
                    
                    blast_results = run_blast_best_hit(
                        config.blast_exe,
                        reference_db,
                        str(fasta_file),
                        config.min_percent_identity,
                        False,
                        host_list,
                        config.threads,
                        config.troubleshooting,
                    )

                    merged_results.update(blast_results)

                counter, bp_number_dict, unmappable_reads, suspected_phage_numbers = count_hits( #TODO: does unmappable reads do anything?
                    merged_results,
                    phage_assembly,
                    host_list,
                    prophage_list,
                    sequence_lengths,
                    config.min_percent_identity,
                    report_files["prophage_analysis"],
                    sample_name
                )

                logger.info(
                    "Purity analysis completed for %s",
                    sample_name,
                )

                if config.export_unmappable_reads:
                    filter_file_by_read_id(
                        directories["unmappable_reads"],
                        assembly_fastq,
                        unmappable_reads
                        )
                    
                    logger.info(
                    "Unmappable reads written to %s",
                    str(directories["unmappable_reads"])
                    )

                

                append_purity_summary(
                    report_files["purity_summary"],
                    sample_name,
                    total_reads,
                    mapping_reads,
                    mapping_bp,
                    suspected_phage_numbers,
                )
                    
def validate_config(config):

    #
    # Full mode requires viralFlye.
    #

    if config.mode == "full":

        if not config.viralflye_hmm_db:

            raise ValueError(
                "--viralflye-hmm-db is required "
                "when running in full mode."
            )

        if not os.path.isfile(config.viralflye_hmm_db):

            raise FileNotFoundError(
                f"HMM database not found: "
                f"{config.viralflye_hmm_db}"
            )

def create_directories(output_dir,export_unmappable_reads):

    output_dir = Path(output_dir)

    directories = {
        "assemblies": output_dir / "assemblies",
        "reports": output_dir / "reports",
        "logs": output_dir / "logs",
    }

    if export_unmappable_reads:
        directories["unmappable_reads"] = output_dir / "unmappable_reads"

    for directory in directories.values():
        directory.mkdir(
            parents=True,
            exist_ok=True,
        )

    return directories

def read_sample_metadata(tsv_file):
    """
    Read sample metadata TSV.

    Returns:
        {
            "sample_name": {
                "phage": "phage_name",
                "host": "host_name",
            }
        }
    """

    sample_metadata = {}

    with open(tsv_file, newline="") as handle:
        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        required_columns = {
            "Sample",
            "Phage",
            "Host",
        }

        missing = required_columns - set(
            reader.fieldnames or []
        )

        if missing:
            raise ValueError(
                f"Missing columns: {', '.join(sorted(missing))}"
            )

        for row in reader:
            sample_metadata[row["Sample"]] = {
                "phage": row["Phage"],
                "host": row["Host"],
            }

    return sample_metadata

def discover_samples(input_path):
    """
    Support:

    raw_data/
        barcode01/
        barcode02/

    sample_directory/
        *.fastq.gz

    sample.fastq.gz

    sample.fastq
    """

    input_path = Path(input_path)

    #
    # Single FASTQ file
    #

    if input_path.is_file():

        sample_name = input_path.stem.replace(".fastq", "")

        return [
            {
                "sample_name": sample_name,
                "input_path": str(input_path),
            }
        ]

    #
    # Single-sample directory
    #

    fastq_files = list(input_path.glob("*.fastq.gz"))
    fastq_files.extend(input_path.glob("*.fastq"))

    if fastq_files:

        return [
            {
                "sample_name": input_path.name,
                "input_path": str(input_path),
            }
        ]

    #
    # Multi-sample directory
    #

    samples = []

    for entry in sorted(input_path.iterdir()):

        if entry.is_dir():

            samples.append(
                {
                    "sample_name": entry.name,
                    "input_path": str(entry),
                }
            )

    return samples
