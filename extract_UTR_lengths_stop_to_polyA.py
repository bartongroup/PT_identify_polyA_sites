#!/usr/bin/env python3

import pysam
import re
import os
import logging
from Bio import SeqIO
import pandas as pd
import argparse


def get_args():
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(description="Extract poly(A) sites and calculate distance from stop codon.",
                                     add_help=False)
    optional = parser.add_argument_group('optional arguments')

    optional.add_argument("--bam", dest='bam',
                          action="store",
                          nargs='+',
                          required=True,
                          type=str,
                          help="List of BAM files to be parsed.")

    optional.add_argument("--output", dest='output',
                          action="store",
                          type=str,
                          required=True,
                          help="Directory to store output files.")

    optional.add_argument("--fasta", dest='fasta',
                          action="store",
                          type=str,
                          required=True,
                          help="FASTA file of the reference genome.")

    optional.add_argument("--reference_transcript", dest='reference_transcript',
                          action="store",
                          type=str,
                          required=True,
                          help="FASTA file of the reference transcripts with UTR.")

    optional.add_argument("--polyA_length", dest='polyA_length',
                          action="store",
                          type=int,
                          default=7,
                          help="Minimum length of poly(A) tail to consider.")

    optional.add_argument("--log", dest='log',
                          action="store",
                          type=str,
                          default="extract_polA_site.log",
                          help="Log file to store the logging information.")

    parser.add_argument("-h", "--help",
                        action="help",
                        default=argparse.SUPPRESS,
                        help="Show this help message and exit.")

    return parser.parse_args()


def index_reference_transcripts(reference_transcript_file):
    """
    Index reference transcripts from a FASTA file.

    Parameters:
    reference_transcript_file (str): Path to the FASTA file containing reference transcript sequences.

    Returns:
    dict: A dictionary of reference transcripts indexed by transcript ID.
    """
    return SeqIO.to_dict(SeqIO.parse(reference_transcript_file, "fasta"))


def find_stop_codon_position(ref_seq, reference_transcript):
    """
    Find the stop codon position in the reference sequence.

    Parameters:
    ref_seq (Seq): Reference sequence as a Biopython Seq object.
    reference_transcript (SeqRecord): Reference transcript as a Biopython SeqRecord object.

    Returns:
    int: Position of the stop codon in the reference sequence.

    Raises:
    ValueError: If the reference transcript is not found in the reference sequence.
    """
    ref_seq_str = str(ref_seq)
    reference_transcript_str = str(reference_transcript.seq)
    match_start = ref_seq_str.find(reference_transcript_str)
    
    if match_start == -1:
        raise ValueError(f"Reference transcript not found in sequence for {reference_transcript.id}")

    stop_codon_position = match_start + len(reference_transcript_str)
    return stop_codon_position


def find_polyA_from_right(seq, polyA_length):
    """
    Search for poly(A) sequences from the right-hand side of the read.

    This function searches for a poly(A) tail starting from the right-hand (3') side of the read.
    If a poly(A) sequence is found, it returns the start position of the poly(A) tail in the read.

    Parameters:
    seq (str): The read sequence.
    polyA_length (int): Minimum length of poly(A) tail to consider.

    Returns:
    int: The start position of the poly(A) tail in the read sequence.
    """
    # Search for poly(A) stretch at the end of the sequence
    right_match = re.search(r'(A{' + str(polyA_length) + ',})$', seq)
    
    if right_match:
        # Return the start position of the poly(A) within the read
        return right_match.start()
    else:
        return None



def extract_polyA_sites_with_distance(bam_file, fasta_file, reference_transcripts, 
                                      polyA_length, output_dir):
    """
    Extract poly(A) sites from a BAM file and calculate distance from stop codon.

    Parameters:
    bam_file (str): Path to the BAM file to be processed.
    fasta_file (str): Path to the FASTA file containing reference genome sequences.
    reference_transcripts (dict): Dictionary of reference transcripts indexed by transcript ID.
    polyA_length (int): Minimum length of poly(A) tail to consider.
    output_dir (str): Directory to output the results.

    Returns:
    None
    """
    logging.info(f"Processing file: {bam_file}")
    bam = pysam.AlignmentFile(bam_file, "rb")
    fasta = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

    output_file = os.path.join(output_dir, f"{os.path.basename(bam_file)}_polyA_distances.tsv")
    with open(output_file, 'w') as out_f:
        out_f.write("transcript_ID\tread_name\tdistance_to_polyA\n")

        for read in bam.fetch():
            if not read.is_unmapped:
                seq = read.query_sequence
                transcript_id = read.reference_name

                if transcript_id not in reference_transcripts:
                    continue

                ref_seq = fasta[transcript_id].seq
                reference_transcript = reference_transcripts[transcript_id]
                
                try:
                    stop_codon_position = find_stop_codon_position(ref_seq, reference_transcript)
                except ValueError as e:
                    logging.debug(e)
                    continue

                # Only consider poly(A) sites after the stop codon
                if read.reference_start < stop_codon_position:
                    match = re.search(r'(A{' + str(polyA_length) + ',})', seq)
                    if match:
                        start_pos = match.start()
                        coordinate = read.reference_start + start_pos

                        # Search again after stop codon if match within CDS
                        if coordinate < stop_codon_position:
                            match = re.search(r'(A{' + str(polyA_length) + ',})',
                                              seq[stop_codon_position - read.reference_start:])
                            if match:
                                start_pos = match.start() + (stop_codon_position - read.reference_start)
                                coordinate = read.reference_start + start_pos

                        if match and coordinate >= stop_codon_position:
                            read_name = read.query_name
                            distance_to_stop = coordinate - stop_codon_position
                            out_f.write(f"{transcript_id}\t{read_name}\t{distance_to_stop}\n")

                        if match and coordinate >= stop_codon_position:
                            read_name = read.query_name
                            polyA_start = read.reference_start + start_pos
                            polyA_length_detected = len(match.group(0))
                            distance_to_stop = coordinate - stop_codon_position
                            old_polyA_start = polyA_start
                            # Step 2: Call the new function to check for poly(A) from the right side
                            right_start_pos = find_polyA_from_right(seq, polyA_length)
                        
                            # If a valid right-side polyA site is found and positions differ, update the start position
                            if right_start_pos is not None and right_start_pos != start_pos:
                                logging.info(f"Right-side poly(A) site found at different position for read {read_name}. Updating poly(A) start position.")
                                polyA_start = read.reference_start + right_start_pos
                                coordinate = read.reference_start + right_start_pos  # Update the coordinate to the new polyA start
                                distance_to_stop = coordinate - stop_codon_position  # Recalculate distance to stop

                                # Update the sequence preceding the polyA tail from the read based on the new right-hand polyA start
                                pre_polyA_seq_from_read = seq[max(0, right_start_pos - distance_to_stop):right_start_pos]

                                # Update the sequence preceding the polyA tail from reference genome
                                pre_polyA_seq_from_ref = str(fasta[transcript_id].seq[stop_codon_position:coordinate]) if distance_to_stop > 0 else ""

                            out_f.write(f"{transcript_id}\t{read_name}\t{distance_to_stop}\n")

                        # Sequence from reference genome
                        pre_polyA_seq_from_ref = str(fasta[transcript_id].seq[stop_codon_position:coordinate]) if distance_to_stop > 0 else ""

    bam.close()
    logging.info(f"Processed {bam_file} and saved results to {output_file}")


def main():
    """
    Main function to handle the overall script execution.
    """
    args = get_args()

    # Setup logging
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s',
                        filename=args.log, filemode='w')
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    logging.getLogger('').addHandler(console)

    # Index the reference transcripts
    reference_transcripts = index_reference_transcripts(args.reference_transcript)

    output_dir = args.output  # Directory to store outputs

    # Process each BAM file
    for bam_file in args.bam:
        extract_polyA_sites_with_distance(bam_file, args.fasta, reference_transcripts,
                                          args.polyA_length, output_dir)


if __name__ == "__main__":
    main()
