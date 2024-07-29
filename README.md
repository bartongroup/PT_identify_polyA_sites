# identify_polyA_sites
script to identify polyA start sites in nanopore direct RNAseq data mapped to a genome



# PolyA Site Extraction and Analysis

This project provides scripts for extracting poly(A) sites from nanopore direct RNAseq data, performing statistical analysis, and identifying differences between wild-type (WT) and mutant (MUT) groups.

## Requirements

- Python 3.6+
- `pysam` (for handling BAM and FASTA files)
- `pandas` (for data manipulation)
- `scipy` (for statistical tests)
- `statsmodels` (for multiple testing correction)
- `nose2` (for running unit tests)

You can install the required Python packages using `pip`:

```bash
pip install pysam pandas scipy statsmodels nose2
```

## file structure

.
├── identify_polyA
│   ├── extract_stop_codon_positions.py
│   ├── parse_gtf.py
├── scripts
│   ├── extract_polyA_sites.py
├── tests
│   ├── test_extract_stop_codon_positions.py
├── README.md



### Scripts Description
1. identify_polyA/parse_gtf.py
This script contains the function parse_gff_gft for parsing GTF/GFF files.

2. identify_polyA/extract_stop_codon_positions.py
This script contains the function extract_stop_codon_positions which extracts stop codon positions from a GTF file.

3. scripts/extract_polyA_sites.py
This is the main script that:

Extracts poly(A) sites from BAM files.
Calculates distances from stop codons.
Performs statistical analysis to identify differences between WT and MUT groups.
Usage
1. Extracting Stop Codon Positions
To extract stop codon positions, you need to use the extract_stop_codon_positions function from the identify_polyA module.

2. Extracting PolyA Sites and Performing Statistical Analysis
Run the main script extract_polyA_sites.py with the following command-line arguments:

--bam: List of BAM files to be parsed (e.g., --bam file1.bam file2.bam file3.bam)
--output: Output TSV file to store the poly(A) sites
--fasta: FASTA file of the reference genome
--gtf: GTF file with genomic annotations
--groups: List of group names corresponding to BAM files (e.g., --groups WT WT WT MUT MUT MUT)
--fdr: False discovery rate threshold for multiple testing correction (default: 0.05)
--log: Log file to store the logging information (default: script.log)
Example Usage

```bash
python scripts/extract_polyA_sites.py --bam file1.bam file2.bam file3.bam --output polyA_sites.tsv --fasta reference.fasta --gtf annotations.gtf --groups WT WT WT MUT MUT MUT --fdr 0.05 --log script.log
```