import os
import unittest
from identify_polyA.extract_stop_codon_positions import extract_stop_codon_positions

class TestExtractStopCodonPositions(unittest.TestCase):

    def setUp(self):
        self.gtf_content = """\
1\tAraport11\tgene\t3631\t5899\t.\t+\t.\ttranscript_id "AT1G01010"; gene_id "AT1G01010";
1\tAraport11\tmRNA\t3631\t5899\t.\t+\t.\ttranscript_id "AT1G01010.1"; gene_id "AT1G01010";
1\tAraport11\tCDS\t3760\t3913\t.\t+\t0\ttranscript_id "AT1G01010.1"; gene_id "AT1G01010";
1\tAraport11\tCDS\t3996\t4276\t.\t+\t2\ttranscript_id "AT1G01010.1"; gene_id "AT1G01010";
1\tAraport11\tCDS\t4486\t4605\t.\t+\t0\ttranscript_id "AT1G01010.1"; gene_id "AT1G01010";
1\tAraport11\tCDS\t4706\t5095\t.\t+\t0\ttranscript_id "AT1G01010.1"; gene_id "AT1G01010";
1\tAraport11\tCDS\t5174\t5326\t.\t+\t0\ttranscript_id "AT1G01010.1"; gene_id "AT1G01010";
1\tAraport11\tCDS\t5439\t5630\t.\t+\t0\ttranscript_id "AT1G01010.1"; gene_id "AT1G01010";
"""
        self.test_gtf_file = 'test.gtf'
        with open(self.test_gtf_file, 'w') as f:
            f.write(self.gtf_content)

    def tearDown(self):
        if os.path.exists(self.test_gtf_file):
            os.remove(self.test_gtf_file)

    def test_extract_stop_codon_positions(self):
        stop_codons = extract_stop_codon_positions(self.test_gtf_file)
        expected_stop_codons = {
            "AT1G01010.1": (5630, '+')
        }
        self.assertEqual(stop_codons, expected_stop_codons)

if __name__ == '__main__':
    unittest.main()
