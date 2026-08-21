"""
run_blast.py

Functions for:

- Creating BLAST databases
- Running BLAST searches
- Calculating percent identity

The pipeline intentionally evaluates how much of a
query sequence can be explained by a reference rather
than using BLAST's native percent identity metric.
"""

import os
import subprocess
import tempfile

from collections import defaultdict
from Bio import SeqIO

import logging

logger = logging.getLogger("phage_pipeline")


def create_reference_database(
    makeblastdb_exe,
    reference_sequence_dir,
    blast_db_dir,
    troubleshooting=False,
):
    """
    Create a BLAST database containing:

    - reference phages
    - host genomes

    Phage genomes are duplicated to allow matching
    across circular genome boundaries.
    """

    logger.info(
        "Creating reference BLAST database"
    )

    os.makedirs(blast_db_dir, exist_ok=True)

    reference_fasta = os.path.join(
        blast_db_dir,
        "reference_sequences.fasta"
    )

    sequence_lengths = {
    "Phage": {},
    "Prophage": {},
    "Host": {},
    }

    phage_list = set()
    host_list = set()
    prophage_list = set()

    with open(reference_fasta, "w") as output_handle:

        for directory in os.listdir(reference_sequence_dir):

            directory_path = os.path.join(
                reference_sequence_dir,
                directory
            )

            for filename in os.listdir(directory_path):

                if not (
                    filename.endswith(".fasta")
                    or filename.endswith(".fna")
                ):
                    continue

                filepath = os.path.join(
                    directory_path,
                    filename
                )
                index = 0
                for record in SeqIO.parse(filepath, "fasta"):

                    record.id = os.path.splitext(filename)[0]

                    if directory == "Phage" or directory == "Prophage":

                        if directory == "Phage":

                            sequence_lengths["Phage"][record.id] = len(record.seq)
                            phage_list.add(record.id)

                        elif directory == "Prophage":

                            index += 1
                            record.id = os.path.splitext(filename)[0]+str(index)
                            prophage_list.add(record.id)
                            sequence_lengths["Prophage"][record.id] = len(record.seq)
                            
                        #
                        # Duplicate circular genomes so
                        # BLAST can align reads spanning
                        # the end/start boundary.
                        #
                        record.seq = record.seq + record.seq

                    elif directory == "Host":

                        host_list.add(record.id)
                        sequence_lengths["Host"][record.id] = len(record.seq)

                    SeqIO.write(
                        record,
                        output_handle,
                        "fasta"
                    )

    blast_db = os.path.join(
        blast_db_dir,
        "reference_db"
    )

    command = (
        f'"{makeblastdb_exe}" '
        f'-in "{reference_fasta}" '
        f'-dbtype nucl '
        f'-out "{blast_db}"'
    )

    logger.debug(command)

    subprocess.run(
        command,
        shell=True,
        check=True
    )

    return (
        blast_db,
        reference_fasta,
        sequence_lengths,
        phage_list,
        host_list,
        prophage_list
    )


def create_assembly_database(
    makeblastdb_exe,
    reference_fasta,
    assembly_files,
    blast_db_dir,
    sequence_lengths,
    troubleshooting=False,
):
    """
    Create a temporary BLAST database containing:

    reference sequences
    + assemblies of interest
    """

    logger.info(
        "Creating assembly BLAST database"
    )

    combined_fasta = os.path.join(
        blast_db_dir,
        "reference_plus_assemblies.fasta"
    )
    
    with open(combined_fasta, "w") as outfile:

        with open(reference_fasta) as infile:
            outfile.write(infile.read())

        for assembly_file in assembly_files:

            for record in SeqIO.parse(
                assembly_file,
                "fasta"
            ):

                sequence_lengths["Phage"][record.id] = len(record.seq)

                #
                # Assemblies are assumed to be
                # circular candidate phages.
                #
                record.seq = record.seq + record.seq

                SeqIO.write(
                    record,
                    outfile,
                    "fasta"
                )

    blast_db = os.path.join(
        blast_db_dir,
        "assembly_db"
    )

    command = (
        f'"{makeblastdb_exe}" '
        f'-in "{combined_fasta}" '
        f'-dbtype nucl '
        f'-out "{blast_db}"'
    )

    logger.debug(command)

    subprocess.run(
        command,
        shell=True,
        check=True
    )

    return blast_db, sequence_lengths


def run_blast_best_hit(
    blast_exe,
    blast_db,
    input_fasta_file,
    min_percent_identity,
    full_assembly,
    host_list,
    threads,
    troubleshooting=False,
):
    """
    Run BLAST and calculate percent identity.

    Percent identity is based on the fraction
    of the query/reference explained by all
    non-overlapping HSPs.
    """

    results = {}

    for record in SeqIO.parse(input_fasta_file, "fasta"):
        results[record.id] = [len(record.seq)]

    #
    # Structure:
    #
    # query_id
    #   └─ hit_id
    #        ├─ query_length
    #        ├─ hit_length
    #        └─ covered_positions
    #

    blast_hits = defaultdict(
        lambda: defaultdict(
            lambda: {
                "query_length": 0,
                "hit_length": 0,
                "covered_positions": set(),
            }
        )
    )

    with tempfile.NamedTemporaryFile(
        suffix=".tsv",
        delete=False
    ) as tmp:

        blast_output = tmp.name

    command = (
        f'"{blast_exe}" '
        f'-db "{blast_db}" '
        f'-query "{input_fasta_file}" '
        f'-out "{blast_output}" '
        f'-outfmt "6 qseqid sseqid qlen slen qstart qend bitscore evalue" '
        f'-num_threads {threads}'
    )

    logger.debug(command)

    subprocess.run(
        command,
        shell=True,
        check=True
    )

    with open(blast_output) as infile:

        for line in infile:

            (
                query_id,
                hit_id,
                query_length,
                hit_length,
                qstart,
                qend,
                bitscore,
                evalue,
            ) = line.rstrip().split("\t")

            bitscore = float(bitscore)

            #
            # Ignore weak alignments.
            #
            if bitscore <= 100:
                continue

            hit = blast_hits[query_id][hit_id]

            hit["query_length"] = int(query_length)

            #
            # Reference phages, prophages and assemblies
            # were duplicated when creating
            # the BLAST database.
            #
            hit["hit_length"] = int(hit_length) / 2

            for position in range(
                int(qstart),
                int(qend) + 1
            ):
                hit["covered_positions"].add(position)

    os.remove(blast_output)

    #
    # Convert coverage information into
    # percent identity values.
    #

    for query_id, hits in blast_hits.items():

        for hit_id, hit in hits.items():

            residues_matched = len(hit["covered_positions"])

            query_gap = abs(
                residues_matched
                - hit["query_length"]
            )

            if hit_id in host_list:

                #
                # Host matches
                #
                # Hosts are much larger than reads.
                # Only determine how much of the
                # read can be explained.
                #

                residues_mismatched = query_gap

            elif full_assembly:
                #
                # Assembly comparisons
                #
                # Both query and hit are expected
                # to represent complete phages.
                #

                hit_gap = abs(
                    residues_matched
                    - hit["hit_length"]
                )

                residues_mismatched = max(
                    query_gap,
                    hit_gap
                )


            else:
                #
                # Read mapping
                #
                # Reads may only cover part of the
                # reference sequence.
                #

                residues_mismatched = query_gap

            if full_assembly:

                expected_length = max(
                    hit["query_length"],
                    hit["hit_length"]
                )

            else:

                expected_length = hit["query_length"]

            percent_identity = (
                100
                * (expected_length - residues_mismatched)
                / expected_length
            )

            if (
                not full_assembly
                and percent_identity
                < min_percent_identity
            ):
                continue

            results[query_id].append(
                {
                    "hit_id": hit_id,
                    "contig_id": query_id,
                    "percent_identity":
                        percent_identity,
                    "residues_matched":
                        residues_matched,
                    "residues_mismatched":
                        residues_mismatched,
                }
            )

    return results
