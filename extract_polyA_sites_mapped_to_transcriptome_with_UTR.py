import pysam
import re
import pandas as pd
import argparse
import os
from scipy.stats import wasserstein_distance, mannwhitneyu
import numpy as np
import logging
from statsmodels.stats.multitest import multipletests
from Bio import SeqIO


def get_args():
    parser = argparse.ArgumentParser(description="Extract poly(A) sites from nanopore direct RNAseq data",
                                     add_help=False)
    file_directory = os.path.realpath(__file__).split("extract_polyA_sites.py")[0]
    optional = parser.add_argument_group('optional arguments')

    optional.add_argument("--bam", dest='bam',
                          action="store",
                          nargs='+',
                          required=True,
                          type=str,
                          help="List of BAM files to be parsed e.g. --bam file1.bam file2.bam file3.bam")

    optional.add_argument("--output", dest='output',
                          action="store",
                          type=str,
                          required=True,
                          help="Output TSV file to store the poly(A) sites")

    optional.add_argument("--fasta", dest='fasta',
                          action="store",
                          type=str,
                          required=True,
                          help="FASTA file of the reference genome")

    optional.add_argument("--reference_transcript", dest='reference_transcript',
                          action="store",
                          type=str,
                          required=True,
                          help="FASTA file of the reference transcripts with UTR")

    optional.add_argument("--groups", dest='groups',
                          action="store",
                          nargs='+',
                          required=True,
                          type=str,
                          help="List of group names corresponding to BAM files e.g. --groups WT WT WT MUT MUT MUT")

    optional.add_argument("--polyA_length", dest='polyA_length',
                          action="store",
                          type=int,
                          default=7,
                          help="Minimum length of poly(A) tail to consider")

    optional.add_argument("--fdr", dest='fdr',
                          action="store",
                          type=float,
                          default=0.05,
                          help="False discovery rate threshold for multiple testing correction")

    optional.add_argument("--log", dest='log',
                          action="store",
                          type=str,
                          default="extract_polA_site.log",
                          help="Log file to store the logging information")

    parser.add_argument("--log-level", dest='log_level', 
                        action="store", 
                        type=str, 
                        default="INFO",
                        help="Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)")


    optional.add_argument("-h", "--help", 
                          action="help", 
                          default=argparse.SUPPRESS,
                          help="Show this help message and exit")

    return parser.parse_args()


def index_reference_transcripts(reference_transcript_file):
    """
    Index reference transcripts from a FASTA file.

    This function reads a FASTA file containing reference transcript sequences and creates 
    a dictionary where each key is a transcript ID and each value is the corresponding sequence record.

    Parameters:
    reference_transcript_file (str): Path to the FASTA file containing reference transcript sequences.

    Returns:
    dict: A dictionary of reference transcripts indexed by transcript ID.
    """
    reference_transcripts = SeqIO.to_dict(SeqIO.parse(reference_transcript_file, "fasta"))
    return reference_transcripts


def find_stop_codon_position(ref_seq, reference_transcript):
    """
    Find the stop codon position in the reference sequence.

    This function identifies the position of the stop codon in a reference sequence 
    (sequence ATG to Stop codon -  i.e. no UTR on 3 or 5 prime end)
    by finding the end position of the reference transcript_wth-UTR-sequence 
    within the reference sequence.

    Note the read should be mapped to the transcript-with-UTR

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
    # for testing
    # if reference_transcript_str.startswith("ATGCGCGCGCGCCGCCGCCGCCGCCG"):
    #     print(ref_seq_str, reference_transcript_str)
    match_start = ref_seq_str.find(reference_transcript_str)
    if match_start == -1:
        raise ValueError(f"Reference transcript not found in the sequence for {reference_transcript.id}")

    # Stop codon position is the end of the matched sequence
    stop_codon_position = match_start + len(reference_transcript_str)
    return stop_codon_position


def extract_polyA_sites(bam_file, fasta_file, reference_transcripts, polyA_length, group):
    """
    Extract poly(A) sites from a BAM file.

    This function processes a BAM file to identify and extract poly(A) sites that are mapped
    to transcriptome sequences. It filters reads to only consider poly(A) sites located after
    the stop codon of the reference transcripts.

    Parameters:
    bam_file (str): Path to the BAM file to be processed.
    fasta_file (str): Path to the FASTA file containing reference genome sequences.
    reference_transcripts (dict): Dictionary of reference transcripts indexed by transcript ID.
    polyA_length (int): Minimum length of poly(A) tail to consider.
    group (str): Group name (e.g., 'WT' or 'MUT') associated with the BAM file.

    Returns:
    list: A list of lists where each inner list contains the following information about a poly(A) site:
        - read_name (str): Name of the read.
        - transcript_id (str): Transcript ID.
        - coordinate (int): Genomic coordinate of the poly(A) site.
        - polyA_start (int): Start position of the poly(A) site within the read.
        - polyA_length (int): Length of the poly(A) tail.
        - pre_polyA_seq_from_read (str): Sequence preceding the poly(A) tail from the read.
        - pre_polyA_seq_from_ref (str): Sequence preceding the poly(A) tail from the reference genome.
        - distance_to_stop (int): Distance from the poly(A) site to the stop codon.
    """
    logging.info(f"Processing file: {bam_file} as {group}")
    bam = pysam.AlignmentFile(bam_file, "rb")
    fasta = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

    polyA_sites = []

    for read in bam.fetch():
        if not read.is_unmapped:
            seq = read.query_sequence
            transcript_id = read.reference_name

            if transcript_id not in reference_transcripts:
                logging.warning(f"Transcript ID {transcript_id} not found in reference transcripts.")
                continue

            ref_seq = fasta[transcript_id].seq
            reference_transcript = reference_transcripts[transcript_id]
            try:
                stop_codon_position = find_stop_codon_position(ref_seq, reference_transcript)
            except ValueError as e:
                logging.warning(e)
                continue

            # Only consider poly(A) sites after the stop codon
            if read.reference_start < stop_codon_position:
                if seq is None:
                    logging.debug(f"Read {read.query_name} has no sequence.")
                    continue
                if not isinstance(seq, str):
                    logging.debug(f"Read {read.query_name} sequence is not a string. Type: {type(seq)}")
                    continue

                match = re.search(r'(A{' + str(polyA_length) + ',})', seq)
                if match:
                    start_pos = match.start()
                    coordinate = read.reference_start + start_pos

                    # If initial match is within CDS, search again after stop codon
                    if coordinate < stop_codon_position:
                        logging.debug(f"Initial polyA site within CDS for read {read.query_name}. Searching again after stop codon.")
                        match = re.search(r'(A{' + str(polyA_length) + ',})', seq[stop_codon_position - 
                                                                                  read.reference_start:])
                        if match:
                            start_pos = match.start() + (stop_codon_position - read.reference_start)
                            coordinate = read.reference_start + start_pos

                    if match and coordinate >= stop_codon_position:
                        read_name = read.query_name
                        polyA_start = read.reference_start + start_pos
                        polyA_length_detected = len(match.group(0))
                        chrom = read.reference_name
                        distance_to_stop = coordinate - stop_codon_position

                        # Sequence from RNAseq read
                        pre_polyA_seq_from_read = seq[max(0, match.start() - distance_to_stop):match.start()]

                        # Sequence from reference genome
                        pre_polyA_seq_from_ref = str(fasta[transcript_id].seq[stop_codon_position:coordinate]) if distance_to_stop > 0 else ""

                        polyA_sites.append([read_name, transcript_id, coordinate, polyA_start, 
                                            polyA_length_detected, pre_polyA_seq_from_read,
                                            pre_polyA_seq_from_ref, distance_to_stop])
                        logging.debug(f"PolyA site found: {match.group()} at position {start_pos} in read {read_name}")
                    elif match:
                        logging.debug(f"PolyA site found within CDS: {match.group()} at position {start_pos} in read {read.query_name}")
                else:
                    logging.debug(f"No PolyA site found in read {read.query_name}")

    bam.close()
    logging.info(f"Extracted {len(polyA_sites)} poly(A) sites from {bam_file}")
    return polyA_sites


def perform_statistical_analysis(polyA_data, fdr_threshold):
    """
    Perform statistical analysis on poly(A) site data.

    This function processes poly(A) site data to compare the distances from the poly(A) sites to the stop codon
    between wild type (WT) and mutant (MUT) groups. It uses the Mann-Whitney U test for the comparison and applies 
    False Discovery Rate (FDR) correction to account for multiple testing.

    Parameters:
    polyA_data (dict): Dictionary containing poly(A) site data for 'WT' and 'MUT' groups.
                       Each value is a list of lists where each inner list contains information about a poly(A) site.
    fdr_threshold (float): The significance level for FDR correction.

    Returns:
    tuple: A tuple containing two elements:
           - significant_results (list): A list of significant results after FDR correction. Each entry contains:
             - transcript_id (str): Transcript ID.
             - u_statistic (float): Mann-Whitney U statistic.
             - p_value (float): Original p-value from the Mann-Whitney U test.
             - pval_corr (float): Corrected p-value after FDR correction.
           - significant_transcript_ids (set): A set of transcript IDs that passed the FDR threshold.
    """
    results = []

    all_transcripts = set([site[1] for group in polyA_data.values() for site in group])
    logging.info(f"Total transcripts found: {len(all_transcripts)}")

    for transcript_id in all_transcripts:
        # site[7] is the distance to stop codon. This is what we are comparing ...
        wt_sites = [site[7] for site in polyA_data['WT'] if site[1] == transcript_id]
        mut_sites = [site[7] for site in polyA_data['MUT'] if site[1] == transcript_id]

        if len(wt_sites) > 1 and len(mut_sites) > 1:
            u_statistic, p_value = mannwhitneyu(wt_sites, mut_sites, alternative='two-sided')
            results.append((transcript_id, u_statistic, p_value))
            logging.info(f"Transcript {transcript_id}: WT count = {len(wt_sites)}, MUT count = {len(mut_sites)}, p-value = {p_value}")
        else:
            logging.debug(f"Transcript {transcript_id}: insufficient data for WT or MUT (WT count = {len(wt_sites)}, MUT count = {len(mut_sites)})")

    logging.info(f"Performed statistical tests on {len(results)} transcripts")

    # Multiple testing correction
    # extract p value from results
    p_values = [result[2] for result in results]
    # multipletests function from the statsmodels library to perform multiple testing correction. return 4 values
    # alpha=fdr_threshold = the threshold
    # method='fdr_bh: Specifies the method for FDR correction. refers to the Benjamini-Hochberg procedure.
    reject, pvals_corrected, _, _ = multipletests(p_values, alpha=fdr_threshold, method='fdr_bh')

    significant_results = [(result[0], result[1], result[2], pval_corr) for result, pval_corr, reject_flag in zip(results, pvals_corrected, reject) if reject_flag]
    explantion_here = """This line creates a list of significant results after FDR correction.
    It uses a list comprehension to iterate over the original results, pvals_corrected, and reject arrays simultaneously.
    For each tuple, it checks if the test was significant after correction (if reject_flag).
    If the test is significant, it includes a tuple in the significant_results list containing:
    result[0]: The first element of the original result tuple (e.g., transcript ID).
    result[1]: The second element of the original result tuple (e.g., test statistic).
    result[2]: The original p-value.
    pval_corr: The corrected p-value."""
    significant_transcript_ids = {result[0] for result in significant_results}
    logging.info(f"Significant transcripts after FDR correction: {len(significant_results)}")
    logging.info(f"Significant transcript IDs: {significant_transcript_ids}")
    return significant_results, significant_transcript_ids


def main():
    args = get_args()

    # Setup logging to file and console
    log_level = getattr(logging, args.log_level.upper(), logging.INFO)
    logging.basicConfig(level=log_level, 
                        format='%(asctime)s - %(levelname)s - %(message)s', 
                        filename=args.log, filemode='w')
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    console.setFormatter(formatter)
    logging.getLogger('').addHandler(console)

    polyA_data = {'WT': [], 'MUT': []}

    # Index the reference transcripts
    reference_transcripts = index_reference_transcripts(args.reference_transcript)

    for bam_file, group in zip(args.bam, args.groups):
        polyA_sites = extract_polyA_sites(bam_file, args.fasta, reference_transcripts, 
                                          args.polyA_length, group)
        polyA_data[group].extend(polyA_sites)

    # Convert results to DataFrame
    columns = ['Read_Name', 'TranscriptID', 'Genomic_Coordinate', 'PolyA_Start', 
               'PolyA_Length', 'Pre_PolyA_Sequence_From_Read', 
               'Pre_PolyA_Sequence_From_Ref', 'Distance_to_Stop']
    polyA_df_wt = pd.DataFrame(polyA_data['WT'], columns=columns)
    polyA_df_mut = pd.DataFrame(polyA_data['MUT'], columns=columns)

    # Save to TSV
    polyA_df_wt.to_csv(f"WT_{args.output}", sep='\t', index=False)
    polyA_df_mut.to_csv(f"MUT_{args.output}", sep='\t', index=False)
    logging.info(f"Poly(A) sites have been extracted and saved to {args.output}")

    # Perform global statistical comparison of poly(A) site locations
    logging.info("Starting global statistical analysis of poly(A) site locations")
    wt_sites = polyA_df_wt['Distance_to_Stop'].values
    mut_sites = polyA_df_mut['Distance_to_Stop'].values

    w_distance = wasserstein_distance(wt_sites, mut_sites)
    logging.info(f"Wasserstein distance between WT and MUT poly(A) sites: {w_distance}")

    # Additional global statistical tests
    u_statistic, p_value = mannwhitneyu(wt_sites, mut_sites, alternative='two-sided')
    logging.info(f"Mann-Whitney U test p-value: {p_value}")

    # Summary statistics
    summary_stats = {
        'WT_Count': len(wt_sites),
        'MUT_Count': len(mut_sites),
        'WT_Mean': np.mean(wt_sites),
        'MUT_Mean': np.mean(mut_sites),
        'WT_Median': np.median(wt_sites),
        'MUT_Median': np.median(wt_sites),
        'WT_Std': np.std(wt_sites),
        'MUT_Std': np.std(wt_sites)
    }

    logging.info(f"Summary Statistics: {summary_stats}")

    # Perform per-transcript statistical comparison of poly(A) site locations
    significant_results, significant_transcript_ids = perform_statistical_analysis(polyA_data, args.fdr)

    # Save significant results
    significant_df = pd.DataFrame(significant_results, columns=['TranscriptID', 'U_Statistic', 'p_value', 'Adjusted_p_value'])
    significant_df.to_csv(f"significant_{args.output}", sep='\t', index=False)
    logging.info(f"Significant transcripts with different poly(A) site locations have been saved to significant_{args.output}")

    logging.info(f"Total significant transcripts: {len(significant_results)}")

    # Global statistical analysis for significant transcripts
    logging.info("Starting global statistical analysis for significant transcripts")
    significant_wt_sites = polyA_df_wt[polyA_df_wt['TranscriptID'].isin(significant_transcript_ids)]['Distance_to_Stop'].values
    significant_mut_sites = polyA_df_mut[polyA_df_mut['TranscriptID'].isin(significant_transcript_ids)]['Distance_to_Stop'].values

    significant_w_distance = wasserstein_distance(significant_wt_sites, significant_mut_sites)
    logging.info(f"Wasserstein distance for significant transcripts between WT and MUT poly(A) sites: {significant_w_distance}")

    # Additional global statistical tests for significant transcripts
    significant_u_statistic, significant_p_value = mannwhitneyu(significant_wt_sites, significant_mut_sites, alternative='two-sided')
    logging.info(f"Mann-Whitney U test p-value for significant transcripts: {significant_p_value}")

    # Summary statistics for significant transcripts
    significant_summary_stats = {
        'WT_Count': len(significant_wt_sites),
        'MUT_Count': len(significant_mut_sites),
        'WT_Mean': np.mean(significant_wt_sites),
        'MUT_Mean': np.mean(significant_mut_sites),
        'WT_Median': np.median(significant_wt_sites),
        'MUT_Median': np.median(significant_mut_sites),
        'WT_Std': np.std(significant_wt_sites),
        'MUT_Std': np.std(significant_mut_sites)
    }

    logging.info(f"Summary Statistics for significant transcripts only: {significant_summary_stats}")


if __name__ == "__main__":
    main()
