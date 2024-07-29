import unittest
import pysam
import os
from identify_polyA.extract_polyA_sites import extract_polyA_sites

class TestExtractPolyASites(unittest.TestCase):

    def setUp(self):
        # Create a dummy BAM file with test data
        self.bam_filename = 'tests/data/test.bam'
        self.fasta_filename = 'tests/data/test.fasta'
        self.create_test_bam()
        self.create_test_fasta()

        # Create a dummy dictionary for stop codons
        self.stop_codons = {
            'chr1': (1000, '+'),
            'chr2': (1000, '-')
        }

    def tearDown(self):
        # Remove the test files after tests
        if os.path.exists(self.bam_filename):
            os.remove(self.bam_filename)
        if os.path.exists(self.bam_filename + '.bai'):
            os.remove(self.bam_filename + '.bai')
        if os.path.exists(self.fasta_filename):
            os.remove(self.fasta_filename)
        if os.path.exists(self.fasta_filename + '.fai'):
            os.remove(self.fasta_filename + '.fai')

    def create_test_bam(self):
        # Create a BAM file with reads containing poly(A) tails
        header = {'HD': {'VN': '1.0'}, 'SQ': [{'LN': 2000, 'SN': 'chr1'}, {'LN': 2000, 'SN': 'chr2'}]}
        with pysam.AlignmentFile(self.bam_filename, "wb", header=header) as bamfile:
            a = pysam.AlignedSegment()
            a.query_name = "read_28833_29006_6945"
            a.query_sequence = "AGCTAGCTAGCTAGCTAAAAAAAAAAA"  # 27 bases
            a.flag = 0
            a.reference_id = 0
            a.reference_start = 950
            a.mapping_quality = 20
            a.cigar = ((0, 27),)  # Match the sequence length
            a.next_reference_id = -1
            a.next_reference_start = -1
            a.template_length = 0
            a.query_qualities = pysam.qualitystring_to_array("<<<<<<<<<<<<<<<<<<<<<<<<<<<")
            bamfile.write(a)

            b = pysam.AlignedSegment()
            b.query_name = "read_28833_29006_6946"
            b.query_sequence = "TTTTTTTTTTTTTTAGCTAGCTAGCTA"  # 27 bases
            b.flag = 0
            b.reference_id = 1
            b.reference_start = 950
            b.mapping_quality = 20
            b.cigar = ((0, 27),)  # Match the sequence length
            b.next_reference_id = -1
            b.next_reference_start = -1
            b.template_length = 0
            b.query_qualities = pysam.qualitystring_to_array("<<<<<<<<<<<<<<<<<<<<<<<<<<<")
            bamfile.write(b)

        # Index the BAM file
        pysam.index(self.bam_filename)

    def create_test_fasta(self):
        # Create a FASTA file with test data, including specific parts with poly(A) sequences
        with open(self.fasta_filename, 'w') as f:
            f.write(">chr1\n")
            f.write("A" * 950 + "TAGCTAGCTAGCTAGCTAAAAAAAAAAA" + "A" * 1050 + "\n")
            f.write(">chr2\n")
            f.write("A" * 950 + "TTTTTTTTTTTTTTAGCTAGCTAGCTA" + "A" * 1050 + "\n")

        # Index the FASTA file
        pysam.faidx(self.fasta_filename)

    def test_extract_polyA_sites(self):
        # Run the extract_polyA_sites function with the test BAM and FASTA files
        polyA_sites = extract_polyA_sites(self.bam_filename, self.fasta_filename, self.stop_codons, 'WT', min_polya_length=10)

        # Check that the output is as expected
        expected_polyA_sites = [
            ['read_28833_29006_6945', 'chr1', 966, 16, 11, 'AGCTAGCTAGCTAGCTAAAAAAAAAAA', -34, False],
            ['read_28833_29006_6946', 'chr2', 966, 13, 13, 'AGCTAGCTAGCTAAAAAAAAAAAAAAA', -34, True]
        ]
        self.assertEqual(polyA_sites, expected_polyA_sites)

if __name__ == '__main__':
    unittest.main()
