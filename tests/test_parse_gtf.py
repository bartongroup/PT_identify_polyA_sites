#!/usr/bin/env python

"""Tests of GTF parsing functionality"""

import os
import unittest
from identify_polyA.parse_gtf import parse_gff_gft

class TestParseGFFGFT(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.file_path = 'data/test.gtf'  # Path to the test GTF file


    def test_parse_gff_gft_file_exists(self):
        """Test if the GTF file exists"""
        self.assertTrue(os.path.isfile(self.file_path), f"File {self.file_path} does not exist")


    def test_parse_gff_gft_output_format(self):
        """Test the format of the output from parse_gff_gft"""
        features = parse_gff_gft(self.file_path)
        self.assertIsInstance(features, list, "Output is not a list")
        for feature in features:
            self.assertIsInstance(feature, tuple, "Feature is not a tuple")
            self.assertEqual(len(feature), 9, "Feature tuple does not have 9 elements")


    def test_parse_gff_gft_content(self):
        """Test the content of the parsed GTF file"""
        features = parse_gff_gft(self.file_path)
        # Test specific feature content (assuming specific content in your test GTF file)
        example_feature = ('1', 'Araport11', 'gene', 3631, 5899, '.', '+', '.', 'ID=AT1G01010;Name=AT1G01010;Note=NAC domain containing protein 1;symbol=NAC001;full_name=NAC domain containing protein 1;computational_description=NAC domain containing protein 1;locus=2200935;locus_type=protein_coding')
        self.assertIn(example_feature, features, "Example feature is not in the parsed features")

if __name__ == '__main__':
    unittest.main()
