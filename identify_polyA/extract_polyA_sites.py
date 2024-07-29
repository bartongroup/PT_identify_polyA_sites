import re
import pysam
import logging
from Bio.Seq import Seq

def extract_polyA_sites(bam_file, fasta_file, stop_codons, group, min_polya_length=10):
    """
    Extract poly(A) sites from the BAM file and calculate distances from stop codons.

    Args:
        bam_file (str): Path to the BAM file.
        fasta_file (str): Path to the reference genome FASTA file.
        stop_codons (dict): Dictionary of stop codon positions and strands.
        group (str): Sample group (e.g., WT or MUT).
        min_polya_length (int): Minimum length of the poly(A) tail to be considered.

    Returns:
        list: A list of poly(A) site information.
    """
    logging.info(f"Processing file: {bam_file} as {group}")
    bam = pysam.AlignmentFile(bam_file, "rb")
    fasta = pysam.FastaFile(fasta_file)
    
    polyA_sites = []
    read_count = 0
    matched_read_count = 0
    polyA_pattern = f"(A{{{min_polya_length},}})$"

    for read in bam.fetch():
        read_count += 1
        if not read.is_unmapped:
            seq = read.query_sequence
            if not isinstance(seq, str):
                seq = str(seq)
            match = re.search(polyA_pattern, seq)
            
            if match:
                print("match = ", match)
                matched_read_count += 1
                read_name = read.query_name
                polyA_start_in_read = match.start()
                polyA_length = len(match.group(0))
                transcript_id = read.reference_name
                chrom = read.reference_name

                # Calculate the genomic coordinate considering the CIGAR string
                coordinate = calculate_genomic_coordinate(read, polyA_start_in_read)

                stop_codon_info = stop_codons.get(transcript_id, None)
                if stop_codon_info:
                    stop_codon_pos, strand = stop_codon_info
                    if strand == '+':
                        distance_from_stop = coordinate - stop_codon_pos
                        region_start = min(stop_codon_pos, coordinate)
                        region_end = max(stop_codon_pos, coordinate)
                        reverse_complement = False
                    else:
                        distance_from_stop = stop_codon_pos - coordinate
                        region_start = min(coordinate, stop_codon_pos)
                        region_end = max(coordinate, stop_codon_pos)
                        reverse_complement = True

                    # Fetch the sequence from the stop codon to the poly(A) site
                    try:
                        pre_polyA_seq = fasta.fetch(chrom, region_start, region_end)
                        if reverse_complement:
                            pre_polyA_seq = str(Seq(pre_polyA_seq).reverse_complement())
                    except KeyError:
                        logging.error(f"Chromosome '{chrom}' not found in FASTA file.")
                        continue

                    polyA_sites.append([read_name, transcript_id, coordinate, polyA_start_in_read, polyA_length, pre_polyA_seq, distance_from_stop, reverse_complement])
    
    logging.info(f"Processed {read_count} reads and found {matched_read_count} reads with poly(A) tails.")
    logging.info(f"Extracted {len(polyA_sites)} poly(A) sites from {bam_file}")
    return polyA_sites




def calculate_genomic_coordinate(read, polyA_start_in_read):
    """
    Calculate the genomic coordinate of the poly(A) start site considering the CIGAR string.

    Args:
        read (pysam.AlignedSegment): The aligned read.
        polyA_start_in_read (int): The start position of the poly(A) tail within the read.

    Returns:
        int: The genomic coordinate of the poly(A) start site.
    """
    genomic_coordinate = read.reference_start
    read_pos = 0

    for operation, length in read.cigartuples:
        print(f"CIGAR Operation: {operation}, Length: {length}")
        
        if operation == 0:  # M (alignment match)
            if read_pos + length > polyA_start_in_read:
                genomic_coordinate += (polyA_start_in_read - read_pos)
                print(f"Final genomic coordinate: {genomic_coordinate}")
                return genomic_coordinate
            genomic_coordinate += length
            read_pos += length
        elif operation == 2:  # D (deletion)
            genomic_coordinate += length
        elif operation == 3:  # N (skipped region)
            genomic_coordinate += length
        elif operation == 1:  # I (insertion)
            read_pos += length
        elif operation == 4:  # S (soft clipping)
            read_pos += length
        elif operation == 5:  # H (hard clipping)
            # Hard clipping does not consume read bases, so no change
            pass
    
    # If the polyA start is not found within the CIGAR operations
    print(f"Final genomic coordinate after loop: {genomic_coordinate}")
    return genomic_coordinate
