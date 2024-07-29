#!/usr/bin/env python3

import os
from identify_polyA.parse_gtf import parse_gff_gft
import re

def extract_stop_codon_positions(gtf_file):
    """
    Extract stop codon positions from the GTF file.

    Args:
        gtf_file (str): Path to the GTF file.

    Returns:
        dict: A dictionary with transcript IDs as keys and stop codon positions as values.
    """
    features = parse_gff_gft(gtf_file)
    stop_codons = {}

    for feature in features:
        seqname, source, feature_type, start, end, score, strand, frame, attribute = feature
        if feature_type == 'CDS':
            transcript_id = re.search(r'transcript_id "([^"]+)"', attribute).group(1)
            if transcript_id not in stop_codons:
                stop_codons[transcript_id] = (end if strand == '+' else end, strand)
            else:
                if strand == '+':
                    stop_codons[transcript_id] = (max(stop_codons[transcript_id][0], end), strand)
                else:
                    stop_codons[transcript_id] = (max(stop_codons[transcript_id][0], end), strand)

    return stop_codons
