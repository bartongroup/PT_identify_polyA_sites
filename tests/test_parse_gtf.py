#!/usr/bin/env python

"""Tests of GTF parsing functionality"""

import unittest
import os
from identify_polyA.parse_gtf import parse_gff_gft

class TestParseGFFGFT(unittest.TestCase):

    def setUp(self):
        self.file_path = 'tests/input/test.gtf'
        self.sample_data = """\
chr1	Araport11	gene	3631	5899	.	+	.	transcript_id "AT1G01010"; gene_id "AT1G01010";
chr1	Araport11	mRNA	3631	5899	.	+	.	transcript_id "AT1G01010.1"; gene_id "AT1G01010";
chr1	Araport11	CDS	3760	3913	.	+	0	transcript_id "AT1G01010.1"; gene_id "AT1G01010";
chr1	Araport11	CDS	3996	4276	.	+	2	transcript_id "AT1G01010.1"; gene_id "AT1G01010";
chr1	Araport11	CDS	4486	4605	.	+	0	transcript_id "AT1G01010.1"; gene_id "AT1G01010";
chr1	Araport11	CDS	4706	5095	.	+	0	transcript_id "AT1G01010.1"; gene_id "AT1G01010";
chr1	Araport11	CDS	5174	5326	.	+	0	transcript_id "AT1G01010.1"; gene_id "AT1G01010";
chr1	Araport11	CDS	5439	5630	.	+	0	transcript_id "AT1G01010.1"; gene_id "AT1G01010";
chr1	Araport11	exon	3631	3913	.	+	.	transcript_id "AT1G01010.1"; gene_id "AT1G01010";
chr1	Araport11	exon	3996	4276	.	+	.	transcript_id "AT1G01010.1"; gene_id "AT1G01010";
"""
        os.makedirs(os.path.dirname(self.file_path), exist_ok=True)
        with open(self.file_path, 'w') as f:
            f.write(self.sample_data)

    def tearDown(self):
        if os.path.exists(self.file_path):
            os.remove(self.file_path)

    def test_parse_gff_gft_output_format(self):
        features = parse_gff_gft(self.file_path)
        expected_output = [
            ('chr1', 'Araport11', 'gene', 3631, 5899, '.', '+', '.', 'transcript_id "AT1G01010"; gene_id "AT1G01010"'),
            ('chr1', 'Araport11', 'mRNA', 3631, 5899, '.', '+', '.', 'transcript_id "AT1G01010.1"; gene_id "AT1G01010"'),
            ('chr1', 'Araport11', 'CDS', 3760, 3913, '.', '+', '0', 'transcript_id "AT1G01010.1"; gene_id "AT1G01010"'),
            ('chr1', 'Araport11', 'CDS', 3996, 4276, '.', '+', '2', 'transcript_id "AT1G01010.1"; gene_id "AT1G01010"'),
            ('chr1', 'Araport11', 'CDS', 4486, 4605, '.', '+', '0', 'transcript_id "AT1G01010.1"; gene_id "AT1G01010"'),
            ('chr1', 'Araport11', 'CDS', 4706, 5095, '.', '+', '0', 'transcript_id "AT1G01010.1"; gene_id "AT1G01010"'),
            ('chr1', 'Araport11', 'CDS', 5174, 5326, '.', '+', '0', 'transcript_id "AT1G01010.1"; gene_id "AT1G01010"'),
            ('chr1', 'Araport11', 'CDS', 5439, 5630, '.', '+', '0', 'transcript_id "AT1G01010.1"; gene_id "AT1G01010"'),
            ('chr1', 'Araport11', 'exon', 3631, 3913, '.', '+', '.', 'transcript_id "AT1G01010.1"; gene_id "AT1G01010"'),
            ('chr1', 'Araport11', 'exon', 3996, 4276, '.', '+', '.', 'transcript_id "AT1G01010.1"; gene_id "AT1G01010"')
        ]
        self.assertEqual(features, expected_output)

if __name__ == '__main__':
    unittest.main()
