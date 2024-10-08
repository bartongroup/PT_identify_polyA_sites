#!/usr/bin/env python3

import pysam
import re
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import os
from scipy.stats import wasserstein_distance, mannwhitneyu
from scipy.stats import iqr, mode
import numpy as np
import logging
from statsmodels.stats.multitest import multipletests
from Bio import SeqIO

# Author P. Thorpe DAG UoD Dundee 2024

print(" ...   libs loaded ...")

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
                logging.debug(f"Transcript ID {transcript_id} not found in reference transcripts.")
                continue

            ref_seq = fasta[transcript_id].seq
            # print(transcript_id)
            reference_transcript = reference_transcripts[transcript_id]
            try:
                stop_codon_position = find_stop_codon_position(ref_seq, reference_transcript)
            except ValueError as e:
                logging.debug(e)
                continue

            # Only consider poly(A) sites after the stop codon
            if read.reference_start < stop_codon_position:
                if seq is None:
                    logging.debug(f"Read {read.query_name} has no sequence.")
                    continue
                if not isinstance(seq, str):
                    logging.debug(f"Read {read.query_name} sequence is not a string. Type: {type(seq)}")
                    continue

                # Step 1: Look for polyA sequence after the stop codon
                match = re.search(r'(A{' + str(polyA_length) + ',})', seq)
                if match:
                    start_pos = match.start()
                    coordinate = read.reference_start + start_pos

                    # Check if the initial match is before the stop codon, re-search if needed
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
                        distance_to_stop = coordinate - stop_codon_position

                        # Sequence preceding polyA tail in the read
                        pre_polyA_seq_from_read = seq[max(0, match.start() - distance_to_stop):match.start()]

                        # Sequence preceding polyA tail from reference genome
                        pre_polyA_seq_from_ref = str(fasta[transcript_id].seq[stop_codon_position:coordinate]) if distance_to_stop > 0 else ""

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
            logging.debug(f"Transcript {transcript_id}: WT count = {len(wt_sites)}, MUT count = {len(mut_sites)}, p-value = {p_value}")
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
    return results, significant_results, significant_transcript_ids


def perform_emd_analysis(polyA_data):
    """
    Perform Earth Mover's Distance (EMD) analysis on poly(A) site data.

    This function processes poly(A) site data to compare the distributions of distances from the poly(A) sites to the stop codon
    between wild type (WT) and mutant (MUT) groups using Earth Mover's Distance (EMD).

    Parameters:
    polyA_data (dict): Dictionary containing poly(A) site data for 'WT' and 'MUT' groups.
                       Each value is a list of lists where each inner list contains information about a poly(A) site.

    Returns:
    list: A list of results for each transcript. Each entry contains:
          - transcript_id (str): Transcript ID.
          - emd (float): Earth Mover's Distance between WT and MUT distributions.
          - wt_mean (float): Mean distance to stop codon for WT.
          - mut_mean (float): Mean distance to stop codon for MUT.
          - wt_median (float): Median distance to stop codon for WT.
          - mut_median (float): Median distance to stop codon for MUT.
          - wt_count (int): Number of WT sites.
          - mut_count (int): Number of MUT sites.
          - wt_min (float): Minimum distance to stop codon for WT.
          - mut_min (float): Minimum distance to stop codon for MUT.
          - wt_max (float): Maximum distance to stop codon for WT.
          - mut_max (float): Maximum distance to stop codon for MUT.
          - wt_mode (float): Mode distance to stop codon for WT.
          - mut_mode (float): Mode distance to stop codon for MUT.
    """
    results = []

    all_transcripts = set([site[1] for group in polyA_data.values() for site in group])
    logging.info(f"Total transcripts found for EMD analysis: {len(all_transcripts)}")

    for transcript_id in all_transcripts:
        # site[7] is the distance to stop codon. This is what we are comparing ...
        wt_sites = [site[7] for site in polyA_data['WT'] if site[1] == transcript_id]
        mut_sites = [site[7] for site in polyA_data['MUT'] if site[1] == transcript_id]

        if len(wt_sites) > 1 and len(mut_sites) > 1:
            emd = wasserstein_distance(wt_sites, mut_sites)
            wt_mean = np.mean(wt_sites)
            mut_mean = np.mean(mut_sites)
            wt_median = np.median(wt_sites)
            mut_median = np.median(mut_sites)
            wt_count = len(wt_sites)
            mut_count = len(mut_sites)
            wt_min = np.min(wt_sites)
            mut_min = np.min(mut_sites)
            wt_max = np.max(wt_sites)
            mut_max = np.max(mut_sites)
            wt_mode_result = mode(wt_sites)
            wt_mode = wt_mode_result.mode[0] if len(wt_mode_result.mode) > 0 else np.nan
            mut_mode_result = mode(mut_sites)
            mut_mode = mut_mode_result.mode[0] if len(mut_mode_result.mode) > 0 else np.nan

            results.append((transcript_id, emd, wt_mean, mut_mean, wt_median, 
                            mut_median, wt_count, mut_count, wt_min, mut_min, 
                            wt_max, mut_max, wt_mode, mut_mode))
            logging.debug(f"Transcript {transcript_id}: EMD = {emd}")
        else:
            logging.debug(f"Transcript {transcript_id}: insufficient data for WT or MUT (WT count = {len(wt_sites)}, MUT count = {len(mut_sites)})")

    logging.info(f"Performed EMD analysis on {len(results)} transcripts")
    return results


def create_violin_plots(wt_sites, mut_sites, significant_wt_sites, 
                        significant_mut_sites, output_prefix):
    """
    Create violin plots for WT and MUT groups, and save them to a PDF file.

    Parameters:
    wt_sites (array-like): Distance to stop codon for WT group for all transcripts.
    mut_sites (array-like): Distance to stop codon for MUT group for all transcripts.
    significant_wt_sites (array-like): Distance to stop codon for WT group for significant transcripts.
    significant_mut_sites (array-like): Distance to stop codon for MUT group for significant transcripts.
    output_prefix (str): Prefix for the output PDF file.
    """
    data_all = pd.DataFrame({
        'Distance_to_Stop': list(wt_sites) + list(mut_sites),
        'Group': ['WT'] * len(wt_sites) + ['MUT'] * len(mut_sites)
    })

    data_significant = pd.DataFrame({
        'Distance_to_Stop': list(significant_wt_sites) + list(significant_mut_sites),
        'Group': ['WT'] * len(significant_wt_sites) + ['MUT'] * len(significant_mut_sites)
    })

    # Plot for all transcripts
    plt.figure(figsize=(10, 12))

    plt.subplot(2, 1, 1)
    sns.violinplot(x='Group', y='Distance_to_Stop', data=data_all)
    plt.title('Violin Plot for All Transcripts')

    plt.subplot(2, 1, 2)
    sns.violinplot(x='Group', y='Distance_to_Stop', data=data_all)
    plt.yscale('log')
    plt.title('Violin Plot for All Transcripts (Log Scale)')

    plt.tight_layout()
    plt.savefig(f"{output_prefix}_all_transcripts_violin_plot.pdf")
    plt.close()

    # Plot for significant transcripts
    plt.figure(figsize=(10, 12))

    plt.subplot(2, 1, 1)
    sns.violinplot(x='Group', y='Distance_to_Stop', data=data_significant)
    plt.title('Violin Plot for Significant Transcripts')

    plt.subplot(2, 1, 2)
    sns.violinplot(x='Group', y='Distance_to_Stop', data=data_significant)
    plt.yscale('log')
    plt.title('Violin Plot for Significant Transcripts (Log Scale)')

    plt.tight_layout()
    plt.savefig(f"{output_prefix}_significant_transcripts_violin_plot.pdf")
    plt.close()


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

    # Summary statistics for all transcripts
    summary_stats = {
        'WT_Count': len(wt_sites),
        'MUT_Count': len(mut_sites),
        'WT_Mean': np.mean(wt_sites),
        'MUT_Mean': np.mean(mut_sites),
        'WT_Median': np.median(wt_sites),
        'MUT_Median': np.median(mut_sites),
        'WT_Std': np.std(wt_sites),
        'MUT_Std': np.std(mut_sites),
        'WT_Min': np.min(wt_sites),
        'MUT_Min': np.min(mut_sites),
        'WT_Max': np.max(wt_sites),
        'MUT_Max': np.max(mut_sites),
        'WT_Mode': mode(wt_sites).mode[0] if len(wt_sites) > 0 else np.nan,
        'MUT_Mode': mode(mut_sites).mode[0] if len(mut_sites) > 0 else np.nan,
        'WT_IQR': iqr(wt_sites),
        'MUT_IQR': iqr(mut_sites)}

    logging.info("Summary Statistics for All Transcripts:")
    logging.info(f"Summary Statistics: {summary_stats}")
    for key, value in summary_stats.items():
        logging.info(f"{key}: {value:.2f}")

    # Perform per-transcript statistical comparison of poly(A) site locations - mann whitney
    results, significant_results, significant_transcript_ids = perform_statistical_analysis(polyA_data, 
                                                                                   args.fdr)
 
    # Perform EMD analysis on poly(A) site locations
    emd_results = perform_emd_analysis(polyA_data)

    # Convert results to DataFrames
    emd_df = pd.DataFrame(emd_results, columns=['TranscriptID', 'EMD', 'WT_Mean', 
                                                'MUT_Mean', 'WT_Median', 'MUT_Median', 
                                                'WT_Count', 'MUT_Count', 'WT_Min', 'MUT_Min', 
                                                'WT_Max', 'MUT_Max', 
                                                'WT_Mode', 'MUT_Mode'])

    stat_df = pd.DataFrame(results, columns=['TranscriptID', 'U_Statistic', 'p_value'])

    # Merge DataFrames on TranscriptID
    combined_df = pd.merge(emd_df, stat_df, on='TranscriptID', how='left')

    # Apply FDR correction to the combined DataFrame
    p_values = combined_df['p_value'].values
    reject, pvals_corrected, _, _ = multipletests(p_values, alpha=args.fdr, method='fdr_bh')
    combined_df['Adjusted_p_value'] = pvals_corrected
    combined_df['Significant'] = reject

    # Save combined results to a single file
    combined_df.to_csv(f"stats_analysis_{args.output}", sep='\t', index=False)
    logging.info(f"Combined EMD and statistical analysis results have been saved to stats_analysis{args.output}")


    # Save significant results
    significant_df = pd.DataFrame(significant_results, columns=['TranscriptID', 'U_Statistic', 
                                                                'p_value', 'Adjusted_p_value'])
    significant_df.to_csv(f"significant_{args.output}", sep='\t', index=False)
    logging.info(f"Significant transcripts with different poly(A) site locations have been saved to significant_{args.output}")

    logging.info(f"Total significant transcripts: {len(significant_results)}")

    logging.info(f"Significant transcripts with different poly(A) site locations have been saved to significant_{args.output}")

    logging.info(f"Total significant transcripts: {len(significant_results)}")

    # Global statistical analysis for significant transcripts
    logging.info("Starting global statistical analysis for significant transcripts")
    significant_wt_sites = polyA_df_wt[polyA_df_wt['TranscriptID'].isin(significant_transcript_ids)]['Distance_to_Stop'].values
    significant_mut_sites = polyA_df_mut[polyA_df_mut['TranscriptID'].isin(significant_transcript_ids)]['Distance_to_Stop'].values

    significant_w_distance = wasserstein_distance(significant_wt_sites, significant_mut_sites)
    logging.info(f"Wasserstein distance for significant transcripts between WT and MUT poly(A) sites: {significant_w_distance}")

    # Additional global statistical tests for significant transcripts
    significant_u_statistic, significant_p_value = mannwhitneyu(significant_wt_sites, 
                                                                significant_mut_sites, 
                                                                alternative='two-sided')
    logging.info(f"Mann-Whitney U test p-value for significant transcripts: {significant_p_value}")

    significant_summary_stats = {
        'WT_Count': len(significant_wt_sites),
        'MUT_Count': len(significant_mut_sites),
        'WT_Mean': np.mean(significant_wt_sites),
        'MUT_Mean': np.mean(significant_mut_sites),
        'WT_Median': np.median(significant_wt_sites),
        'MUT_Median': np.median(significant_mut_sites),
        'WT_Std': np.std(significant_wt_sites),
        'MUT_Std': np.std(significant_mut_sites),
        'WT_Min': np.min(significant_wt_sites),
        'MUT_Min': np.min(significant_mut_sites),
        'WT_Max': np.max(significant_wt_sites),
        'MUT_Max': np.max(significant_mut_sites),
        'WT_Mode': mode(significant_wt_sites).mode[0] if len(significant_wt_sites) > 0 else np.nan,
        'MUT_Mode': mode(significant_mut_sites).mode[0] if len(significant_mut_sites) > 0 else np.nan,
        'WT_IQR': iqr(significant_wt_sites),
        'MUT_IQR': iqr(significant_mut_sites)}

    logging.info(f"Summary Statistics for significant transcripts only: {significant_summary_stats}")
    for key, value in significant_summary_stats.items():
        logging.info(f"{key}: {value:.2f}")

    # call the plot function
    create_violin_plots(wt_sites, mut_sites, 
                        significant_wt_sites, 
                        significant_mut_sites, 
                        args.output)

if __name__ == "__main__":
    main()
