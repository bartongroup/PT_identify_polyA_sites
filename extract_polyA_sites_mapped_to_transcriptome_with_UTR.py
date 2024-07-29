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
                          default=10,
                          help="Minimum length of poly(A) tail to consider")

    optional.add_argument("--fdr", dest='fdr',
                          action="store",
                          type=float,
                          default=0.05,
                          help="False discovery rate threshold for multiple testing correction")

    optional.add_argument("--log", dest='log',
                          action="store",
                          type=str,
                          default="script.log",
                          help="Log file to store the logging information")

    optional.add_argument("-h", "--help", action="help", default=argparse.SUPPRESS,
                          help="Show this help message and exit")

    return parser.parse_args()

def index_reference_transcripts(reference_transcript_file):
    reference_transcripts = SeqIO.to_dict(SeqIO.parse(reference_transcript_file, "fasta"))
    return reference_transcripts

def find_stop_codon_position(ref_seq, reference_transcript):
    ref_seq_str = str(ref_seq)
    reference_transcript_str = str(reference_transcript.seq)
    
    match_start = ref_seq_str.find(reference_transcript_str)
    if match_start == -1:
        raise ValueError(f"Reference transcript not found in the sequence for {reference_transcript.id}")

    # Stop codon position is the end of the matched sequence
    stop_codon_position = match_start + len(reference_transcript_str)
    return stop_codon_position

def extract_polyA_sites(bam_file, fasta_file, reference_transcripts, polyA_length, group):
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
            if read.reference_start >= stop_codon_position:
                match = re.search(r'(A{' + str(polyA_length) + ',})$', seq)
                if match:
                    read_name = read.query_name
                    polyA_start = read.reference_start + match.start()
                    polyA_length = len(match.group(0))
                    chrom = read.reference_name
                    coordinate = read.reference_start + match.start()

                    distance_to_stop = coordinate - stop_codon_position

                    # Sequence from RNAseq read
                    pre_polyA_seq_from_read = seq[match.start() - distance_to_stop:match.start()]

                    # Sequence from reference genome
                    pre_polyA_seq_from_ref = str(fasta[transcript_id].seq[stop_codon_position:coordinate]) if distance_to_stop > 0 else ""

                    polyA_sites.append([read_name, transcript_id, coordinate, polyA_start, polyA_length, pre_polyA_seq_from_read, pre_polyA_seq_from_ref, distance_to_stop])

    logging.info(f"Extracted {len(polyA_sites)} poly(A) sites from {bam_file}")
    return polyA_sites

def perform_statistical_analysis(polyA_data, fdr_threshold):
    results = []

    all_transcripts = set([site[1] for group in polyA_data.values() for site in group])
    logging.info(f"Total transcripts found: {len(all_transcripts)}")

    for transcript_id in all_transcripts:
        wt_sites = [site[7] for site in polyA_data['WT'] if site[1] == transcript_id]
        mut_sites = [site[7] for site in polyA_data['MUT'] if site[1] == transcript_id]

        if len(wt_sites) > 1 and len(mut_sites) > 1:
            u_statistic, p_value = mannwhitneyu(wt_sites, mut_sites, alternative='two-sided')
            results.append((transcript_id, u_statistic, p_value))
            logging.info(f"Transcript {transcript_id}: WT count = {len(wt_sites)}, MUT count = {len(mut_sites)}, p-value = {p_value}")
        else:
            logging.info(f"Transcript {transcript_id}: insufficient data for WT or MUT (WT count = {len(wt_sites)}, MUT count = {len(mut_sites)})")

    logging.info(f"Performed statistical tests on {len(results)} transcripts")

    # Multiple testing correction
    p_values = [result[2] for result in results]
    reject, pvals_corrected, _, _ = multipletests(p_values, alpha=fdr_threshold, method='fdr_bh')

    significant_results = [(result[0], result[1], result[2], pval_corr) for result, pval_corr, reject_flag in zip(results, pvals_corrected, reject) if reject_flag]
    logging.info(f"Significant transcripts after FDR correction: {len(significant_results)}")
    return significant_results

def main():
    args = get_args()

    # Setup logging to file and console
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', 
                        filename=args.log, filemode='w')
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(levellevel)s - %(message)s')
    console.setFormatter(formatter)
    logging.getLogger('').addHandler(console)

    polyA_data = {'WT': [], 'MUT': []}

    # Index the reference transcripts
    reference_transcripts = index_reference_transcripts(args.reference_transcript)

    for bam_file, group in zip(args.bam, args.groups):
        polyA_sites = extract_polyA_sites(bam_file, args.fasta, reference_transcripts, args.polyA_length, group)
        polyA_data[group].extend(polyA_sites)

    # Convert results to DataFrame
    columns = ['Read_Name', 'TranscriptID', 'Genomic_Coordinate', 'PolyA_Start', 'PolyA_Length', 'Pre_PolyA_Sequence_From_Read', 'Pre_PolyA_Sequence_From_Ref', 'Distance_to_Stop']
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
        'MUT_Mean': np.mean(wt_sites),
        'WT_Median': np.median(wt_sites),
        'MUT_Median': np.median(wt_sites),
        'WT_Std': np.std(wt_sites),
        'MUT_Std': np.std(wt_sites)
    }

    logging.info(f"Summary Statistics: {summary_stats}")

    # Perform per-transcript statistical comparison of poly(A) site locations
    significant_results = perform_statistical_analysis(polyA_data, args.fdr)

    # Save significant results
    significant_df = pd.DataFrame(significant_results, columns=['TranscriptID', 'U_Statistic', 'p_value', 'Adjusted_p_value'])
    significant_df.to_csv(f"significant_{args.output}", sep='\t', index=False)
    logging.info(f"Significant transcripts with different poly(A) site locations have been saved to significant_{args.output}")

    logging.info(f"Total significant transcripts: {len(significant_results)}")

if __name__ == "__main__":
    main()
