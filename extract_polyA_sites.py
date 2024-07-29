import pysam
import re
import pandas as pd
import argparse
import os
from scipy.stats import wasserstein_distance, mannwhitneyu
import numpy as np
import logging
from statsmodels.stats.multitest import multipletests
from identify_polyA.parse_gtf import parse_gff_gft
from identify_polyA.extract_stop_codon_positions import extract_stop_codon_positions

def get_args():
    """
    Parse and return the command line arguments.

    Returns:
        argparse.Namespace: Parsed command line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Extract poly(A) sites from nanopore direct RNAseq data",
        add_help=False
    )
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

    optional.add_argument("--gtf", dest='gtf',
                          action="store",
                          type=str,
                          required=True,
                          help="GTF file with genomic annotations")

    optional.add_argument("--groups", dest='groups',
                          action="store",
                          nargs='+',
                          required=True,
                          type=str,
                          help="List of group names corresponding to BAM files e.g. --groups WT WT WT MUT MUT MUT")

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


def extract_polyA_sites(bam_file, fasta_file, stop_codons, group):
    """
    Extract poly(A) sites from the BAM file and calculate distances from stop codons.

    Args:
        bam_file (str): Path to the BAM file.
        fasta_file (str): Path to the reference genome FASTA file.
        stop_codons (dict): Dictionary of stop codon positions and strands.
        group (str): Sample group (e.g., WT or MUT).

    Returns:
        list: A list of poly(A) site information.
    """
    logging.info(f"Processing file: {bam_file} as {group}")
    bam = pysam.AlignmentFile(bam_file, "rb")
    fasta = pysam.FastaFile(fasta_file)
    
    polyA_sites = []
    read_count = 0
    matched_read_count = 0

    for read in bam.fetch():
        read_count += 1
        if not read.is_unmapped:
            seq = read.query_sequence
            logging.debug(f"Read {read.query_name}: Sequence length = {len(seq)}")

            if not isinstance(seq, str):
                seq = str(seq)
            match = re.search(r'(A{10,})$', seq)
            
            if match:
                matched_read_count += 1
                read_name = read.query_name
                polyA_start = read.reference_start + match.start()
                polyA_length = len(match.group(0))
                transcript_id = read.reference_name
                chrom = read.reference_name
                coordinate = read.reference_start + match.start()
                
                if coordinate >= 20:
                    region_start = coordinate - 20
                    region_end = coordinate
                    try:
                        pre_polyA_seq = fasta.fetch(chrom, region_start, region_end)
                    except KeyError:
                        logging.error(f"Chromosome '{chrom}' not found in FASTA file.")
                        continue
                else:
                    try:
                        pre_polyA_seq = fasta.fetch(chrom, 0, coordinate)
                    except KeyError:
                        logging.error(f"Chromosome '{chrom}' not found in FASTA file.")
                        continue

                stop_codon_info = stop_codons.get(transcript_id, None)
                if stop_codon_info:
                    stop_codon_pos, strand = stop_codon_info
                    if strand == '+':
                        distance_from_stop = coordinate - stop_codon_pos
                    else:
                        distance_from_stop = stop_codon_pos - coordinate
                    polyA_sites.append([read_name, transcript_id, coordinate, polyA_start, polyA_length, pre_polyA_seq, distance_from_stop])
    
    logging.info(f"Processed {read_count} reads and found {matched_read_count} reads with poly(A) tails.")
    logging.info(f"Extracted {len(polyA_sites)} poly(A) sites from {bam_file}")
    return polyA_sites


def perform_statistical_analysis(polyA_data, fdr_threshold):
    """
    Perform statistical analysis on poly(A) site data to find significant differences.

    Args:
        polyA_data (dict): Dictionary of poly(A) site data for WT and MUT groups.
        fdr_threshold (float): False discovery rate threshold for multiple testing correction.

    Returns:
        list: A list of significant results after multiple testing correction.
    """
    results = []

    all_transcripts = set([site[1] for group in polyA_data.values() for site in group])
    logging.info(f"Total transcripts found: {len(all_transcripts)}")

    for transcript_id in all_transcripts:
        wt_sites = [site[6] for site in polyA_data['WT'] if site[1] == transcript_id]
        mut_sites = [site[6] for site in polyA_data['MUT'] if site[1] == transcript_id]

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
    """
    Main function to orchestrate the extraction and analysis of poly(A) sites.
    """
    args = get_args()
    
    # Setup logging to file and console
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', 
                        filename=args.log, filemode='w')
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    console.setFormatter(formatter)
    logging.getLogger('').addHandler(console)

    # Extract stop codon positions from GTF file
    stop_codons = extract_stop_codon_positions(args.gtf)

    polyA_data = {'WT': [], 'MUT': []}

    for bam_file, group in zip(args.bam, args.groups):
        polyA_sites = extract_polyA_sites(bam_file, args.fasta, stop_codons, group)
        polyA_data[group].extend(polyA_sites)

    # Convert results to DataFrame
    columns = ['Read_Name', 'TranscriptID', 'Genomic_Coordinate', 'PolyA_Start', 'PolyA_Length', 'Pre_PolyA_Sequence', 'Distance_From_Stop']
    polyA_df_wt = pd.DataFrame(polyA_data['WT'], columns=columns)
    polyA_df_mut = pd.DataFrame(polyA_data['MUT'], columns=columns)
    
    # Save to TSV
    polyA_df_wt.to_csv(f"WT_{args.output}", sep='\t', index=False)
    polyA_df_mut.to_csv(f"MUT_{args.output}", sep='\t', index=False)
    logging.info(f"Poly(A) sites have been extracted and saved to {args.output}")

    # Check if there are any poly(A) sites
    if polyA_df_wt.empty or polyA_df_mut.empty:
        logging.error("No poly(A) sites extracted. Ensure the input BAM files and GTF annotations are correct.")
        return

    # Perform global statistical comparison of poly(A) site locations
    logging.info("Starting global statistical analysis of poly(A) site locations")
    wt_sites = polyA_df_wt['Distance_From_Stop'].values
    mut_sites = polyA_df_mut['Distance_From_Stop'].values

    if len(wt_sites) == 0 or len(mut_sites) == 0:
        logging.error("No poly(A) sites found in one or both groups. Cannot perform statistical analysis.")
        return

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
        'MUT_Median': np.median(mut_sites),
        'WT_Std': np.std(wt_sites),
        'MUT_Std': np.std(mut_sites)
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
