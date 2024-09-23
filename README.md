
**PT_identify_polyA_sites** is a Python-based tool designed to identify and analyze poly(A) sites in RNA-seq data. 
This tool processes BAM files to extract poly(A) sites and performs statistical analysis to compare the poly(A) 
site distributions between wild-type (WT) and mutant (MUT) samples.

## Quick start - example usage

Nanopore direct RNAseq reads need to be mapped to the transcriptome with the UTR sequences. This would correspond to this in the example "transcript_with_UTR.fa". Then for the input for this tool --reference_transcript  need sot the the same sequence_ID names but with out the UTR.
These sequnces are used to index where the stop codon is. 

## Example

```bash
python extract_polyA_sites_mapped_to_transcriptome_with_UTR.py --bam transcriptome_tests/test.bam 
--fasta transcriptome_tests/transcript_with_UTR.fa --polyA_length 7 
--reference_transcript transcriptome_tests/transcript.fa
--fdr 0.05 --log analysis.log --log-level INFO --output results.tsv

```

## Features

- Extract poly(A) sites from BAM files mapped to transcriptome sequences.
- Compare the distances of poly(A) start site to the stop codon between WT and MUT groups.
- Perform statistical tests, including Mann-Whitney U test and False Discovery Rate (FDR) correction. And Earth-mover distance
- Generate summary statistics for both overall and significant transcripts.


# PolyA Site Extraction and Analysis

This project provides scripts for extracting poly(A) sites from nanopore direct RNAseq data, performing statistical analysis, and identifying differences between wild-type (WT) and mutant (MUT) groups.

## Requirements

- Python 3.6+
- `pysam` (for handling BAM and FASTA files)
- `pandas` (for data manipulation)
- `scipy` (for statistical tests)
- `statsmodels` (for multiple testing correction)
- `nose2` (for running unit tests)
- `biopython` for the transcriptome approach 
- `matplotlib` for the violin plot 
- `seaborn` for the violin plot 
You can install the required Python packages using `pip`:

```bash
pip install pysam pandas scipy statsmodels nose2 biopython matplotlib seaborn
```

## unit test
in the parent directory run

```bash
nose2
```

# Purpose and Objective:
The script is designed to identify polyadenylation (polyA) sites from RNA sequencing (RNA-seq) data, specifically from nanopore direct RNA sequencing. PolyA tails are important in RNA biology as they indicate the end of an mRNA transcript and play key roles in mRNA stability and degradation. By identifying where the polyA tail starts in reads aligned to reference transcripts, the script helps researchers understand the transcriptome's 3' untranslated regions (UTRs) and the locations of polyA sites.

## How It Works:
The script processes BAM files, which contain RNA-seq reads aligned to reference transcripts, and identifies polyA sites by scanning for stretches of adenine (A) residues. 

It specifically looks for polyA tails that occur after the stop codon in the mRNA sequence, as polyA tails are generally located in the 3' UTR after the coding sequence. The script starts by identifying potential polyA sequences based on a user-defined minimum length of consecutive adenines (e.g., 7 or 10 As in a row). Additionally, to handle cases where the polyA stretch might appear on the right-hand side of the read, the script also looks for polyA sequences from the end of the read and updates the polyA start position if necessary.

## The script extracts various information for each polyA site it detects:

Read Name: The name of the read in which the polyA site was found.
Transcript ID: The transcript to which the read was aligned.
Genomic Coordinate: The genomic coordinate of the polyA start.
PolyA Start: The starting position of the polyA site in the read.
PolyA Length: The length of the polyA tail detected in the read.
Pre-PolyA Sequence from Read: The sequence of nucleotides in the read just before the polyA tail starts.
Pre-PolyA Sequence from Reference: The corresponding sequence from the reference genome or transcript that aligns to the pre-polyA region in the read.
Distance to Stop Codon: The distance from the polyA site to the stop codon in the reference transcript.
Why It’s Done:
## Understanding where polyA tails occur is crucial for several reasons:

Transcript Processing: The polyA tail is involved in the maturation and stability of mRNA. Identifying its position helps in understanding mRNA post-transcriptional modifications.
UTR Annotation: PolyA sites generally occur in the 3' UTR of mRNAs, and their precise location can help researchers refine transcript annotations.
Differences Between Conditions: In comparative studies (e.g., wild type vs mutant), polyA sites may shift or differ in length, which could have biological significance. The script allows researchers to collect polyA site data for further analysis and statistical comparison.
Inputs:
BAM Files (--bam): These files contain the RNA-seq reads aligned to a reference transcriptome.
FASTA File for the Reference Genome (--fasta): This file contains the reference genome sequences in FASTA format.
Reference Transcripts with UTRs (--reference_transcript): This FASTA file contains reference transcript sequences, including their untranslated regions (UTRs), which are crucial for detecting polyA sites after the stop codon.
Groups (--groups): This allows specifying different experimental conditions (e.g., wild type and mutant) for polyA site comparison.
PolyA Length (--polyA_length): The minimum number of consecutive adenines required to identify a polyA tail.
Log Level (--log-level): Defines the logging verbosity for debugging purposes.
Outputs:
The script produces several output files in TSV format, which summarize the detected polyA sites:

TSV File for PolyA Sites: A file listing all identified polyA sites, their genomic coordinates, and the corresponding details (e.g., pre-polyA sequences and distances to stop codons). The file is generated separately for different experimental groups (e.g., WT_output.tsv for wild type).
Additional Analysis Files: If further analyses are run, such as statistical comparisons between groups, additional files can be produced, showing significant differences in polyA site locations.
Why This Is Useful:
This script is especially useful for researchers studying post-transcriptional regulation, mRNA decay, or alternative polyadenylation. By providing precise positions of polyA sites, researchers can investigate how mRNA stability and expression levels might be influenced under different conditions. Moreover, the script helps in understanding UTR dynamics and refining transcript annotations in various organisms or cell types.

This modular and automated approach allows for a high-throughput analysis of polyA sites, reducing manual effort and ensuring consistent, reliable results across large datasets.


## file structure

```bash

.
├── LICENSE
├── MUT_output.test
├── README.md
├── analysis.log
├── data
├── extract_UTR_lengths_stop_to_polyA.py
├── extract_polA_site.log
├── extract_polyA_sites.py
├── extract_polyA_sites_mapped_to_transcriptome_with_UTR.py
├── identify_polyA
│   ├── __pycache__
│   │   ├── extract_polyA_sites.cpython-310.pyc
│   │   ├── extract_stop_codon_positions.cpython-310.pyc
│   │   ├── parse_gff_gft.cpython-310.pyc
│   │   └── parse_gtf.cpython-310.pyc
│   ├── extract_polyA_sites.py
│   ├── extract_stop_codon_positions.py
│   ├── parse_gff_gft.py
│   └── parse_gtf.py
├── script.log
├── test_data
│   ├── reads.fa
│   ├── reads.fq
│   ├── test.bam
│   ├── test.bam.bai
│   ├── test.genome.fasta
│   ├── test.genome.fasta.fai
│   ├── test.gtf
│   └── transcript.fa
├── tests
│   ├── __pycache__
│   │   ├── test_extract_polyA_sites.cpython-310.pyc
│   │   ├── test_extract_stop_codon_positions.cpython-310.pyc
│   │   └── test_parse_gtf.cpython-310.pyc
│   ├── data
│   ├── input
│   ├── test_extract_polyA_sites.py
│   ├── test_extract_stop_codon_positions.py
│   └── test_parse_gtf.py
└── transcriptome_tests
    ├── reads.fa
    ├── test.bam
    ├── test.bam.bai
    ├── transcript.fa
    └── transcript_with_UTR.fa

```
### transcritpome approach (genome version is NOT WORKING YET!)

```
python extract_polyA_sites_mapped_to_transcriptome_with_UTR.py \
    --bam <BAM_FILES> \
    --fasta <FASTA_FILE> \
    --reference_transcript <REFERENCE_TRANSCRIPT_FILE> \
    --polyA_length: Minimum length of poly(A) tail to consider \
    --groups <GROUP_NAMES> \
    --log: Log file to store the logging information. '
    --log-level: Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL). \
    --output <OUTPUT_PREFIX>

```


### Arguments more indepth info 
    --bam: List of BAM files to be parsed (e.g., --bam file1.bam file2.bam file3.bam).
    --fasta: FASTA file of the reference genome.
    --reference_transcript: FASTA file of the reference transcripts with UTR.
    --groups: List of group names corresponding to BAM files (e.g., --groups WT WT WT MUT MUT MUT).
    --output: Output prefix for the generated TSV files.
    --polyA_length: Minimum length of poly(A) tail to consider (default: 10).
    --fdr: False discovery rate threshold for multiple testing correction (default: 0.05).
    --log: Log file to store the logging information (default: script.log).

### Output
The script generates three main output files:

1) WT_<output_prefix>.tsv: Poly(A) site information for the WT group.
2) MUT_<output_prefix>.tsv: Poly(A) site information for the MUT group.
3) significant_<output_prefix>.tsv: Significant transcripts with different poly(A) site locations between WT and MUT groups after FDR correction.
Each TSV file contains the following columns:

- `Read_Name`: Name of the read.
- `TranscriptID`: Identifier of the transcript.
- `Genomic_Coordinate`: Genomic coordinate of the poly(A) site.
- `PolyA_Start`: Start position of the poly(A) tail.
- `PolyA_Length`: Length of the poly(A) tail.
- `Pre_PolyA_Sequence_From_Read`: Sequence preceding the poly(A) site from the RNAseq read.
- `Pre_PolyA_Sequence_From_Ref`: Sequence preceding the poly(A) site from the reference genome.
- `Distance_to_Stop`: Distance from the poly(A) site to the stop codon.

### Statistical Analysis
The tool performs statistical analysis to compare poly(A) site distributions between WT and MUT groups. The analysis includes the following steps:

1) Extract Poly(A) Sites:
 - Identify poly(A) sites from BAM files.
 - Filter reads to only consider poly(A) sites located after the stop codon.
2) Per-Transcript Statistical Comparison:
 - Use Mann-Whitney U test to compare the distances from poly(A) sites to the stop codon between WT and MUT groups.
 - Apply FDR correction to account for multiple testing.
3) Global Statistical Analysis:
 - Calculate the Wasserstein distance between WT and MUT poly(A) sites.
 - Perform additional global statistical tests, including Mann-Whitney U test.
 - Generate summary statistics, including count, mean, median, standard deviation, minimum, maximum, mode, and interquartile range (IQR).


Logging
The script logs its progress and important information to a log file specified by the --log argument (default: script.log). The log file contains details about the processing of BAM files, the number of poly(A) sites extracted, and the results of statistical tests.

## Example

```bash
python extract_polyA_sites_mapped_to_transcriptome_with_UTR.py --bam transcriptome_tests/test.bam --fasta transcriptome_tests/transcript_with_UTR.fa --polyA_length 7 --fdr 0.05 --log analysis.log --log-level INFO --output results.tsv --reference_transcript transcriptome_tests/transcript.fa
```


# Generating the test files and example bam file. 

You can use the gffread tool from the gffread suite (part of the Cufflinks/gffcompare tools) to extract both the CDS (coding sequence) and the transcripts with UTRs from a GFF/GTF file. 
```bash
gffread input.gff -g genome.fa -x cds_output.fa -w transcripts_with_UTR_output.fa

```

## Map the FASTA sequences to the reference and convert to BAM directly
```bash
minimap2 -a -x sr -k 7 -w 1 -A 1 -B 2 -O 2,16 -E 4,1 -s 40 --secondary=no transcript_with_UTR.fa reads.fa | samtools view -b -o tmp.bam
```

## Sort the BAM file
```bash
samtools sort -o test.bam tmp.bam

```

## Index the sorted BAM file
```bash
samtools index test.bam
```

testing command

```bash
python extract_polyA_sites_mapped_to_transcriptome_with_UTR.py --bam transcriptome_tests/test.bam --reference_transcript transcriptome_tests/transcript.fa --fasta transcriptome_tests/transcript_with_UTR.fa --group WT --output output.test

```


# extract_UTR_lengths_stop_to_polyA.py

This is a script that perfomrs a similar task to the main one but JUST reports the data as 

 transcript_id  read_name   distance_to_stop

The the user can decide what to do with the data. 


##################################


# Genome version (not yet working) -  unlikely to take this further .. 

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

### Example Usage

```bash
python scripts/extract_polyA_sites.py --bam file1.bam file2.bam file3.bam --output polyA_sites.tsv --fasta reference.fasta --gtf annotations.gtf --groups WT WT WT MUT MUT MUT --fdr 0.05 --log script.log
```

