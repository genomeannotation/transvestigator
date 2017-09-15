# coding=utf-8
import io
import unittest
from unittest.mock import Mock, patch, PropertyMock
from src.mrna import Mrna

class TestMrna(unittest.TestCase):

    def test_from_gff_feature_success(self):
        gff_mrna = Mock()
        gff_mrna.type = "mRNA"
        
        mrna = Mrna.from_gff_feature(gff_mrna)
        self.assertTrue(mrna)

    def test_from_gff_features_fails(self):
        gff_mrna = Mock()
        gff_mrna.type = "asdf"
        
        mrna = Mrna.from_gff_feature(gff_mrna)
        self.assertFalse(mrna)

    #### MAKE POSITIVE TESTS ####

    def test_make_positive(self):
        seq_len = 8

        mrna = Mrna(start=1, end=7, strand='-')
        cds = Mock()
        cds.start = 1
        cds.end = 7
        cds.strand = '-'
        exon = Mock()
        exon.start = 1
        exon.end = 7
        cds.strand = '-'

        mrna.children = {'cds':[cds], 'exon':[exon]}

        mrna.make_positive(seq_len)

        self.assertEqual(mrna.start, 2)
        self.assertEqual(mrna.end, 8)
        self.assertEqual(mrna.strand, '+')
        self.assertEqual(cds.start, 2)
        self.assertEqual(cds.end, 8)
        self.assertEqual(cds.strand, '+')
        self.assertEqual(exon.start, 2)
        self.assertEqual(exon.end, 8)
        self.assertEqual(exon.strand, '+')

    #### MATCH CDS AND EXON END TESTS ####

    def test_match_cds_and_exon_end(self):
        mrna = Mrna()
        cds = Mock()
        cds.start = 1
        cds.end = 5
        exon = Mock()
        exon.start = 1
        exon.end = 6

        mrna.children = {'cds':[cds], 'exon':[exon]}

        mrna.match_cds_and_exon_end()
        self.assertEquals(cds.end, 6)

    def test_match_cds_and_exon_end_does_nothing_if_stop_codon_present(self):
        mrna = Mrna()
        cds = Mock()
        cds.start = 1
        cds.end = 5
        exon = Mock()
        exon.start = 1
        exon.end = 6

        mrna.children = {'cds':[cds], 'exon':[exon], 'stop_codon':[Mock()]}

        mrna.match_cds_and_exon_end()
        self.assertEquals(cds.end, 5)

    #### STARTS AND STOPS TESTS ####
    
    def test_create_starts_and_stops_creates_a_start(self):
        mrna = Mrna()
        cds = Mock()
        cds.start = 1
        cds.end = 6
        cds.strand = '+'
        mrna.children = {'cds':[cds]}
        mrna.attributes = {'ID':'m.1234'}
        
        mrna.create_starts_and_stops('ATGNNN')

        self.assertTrue('start_codon' in mrna)

    def test_create_starts_and_stops_creates_a_start_reverse_complement(self):
        mrna = Mrna()
        cds = Mock()
        cds.start = 1
        cds.end = 6
        cds.strand = '-'
        mrna.children = {'cds':[cds]}
        mrna.attributes = {'ID':'m.1234'}
        
        mrna.create_starts_and_stops('NNNCAT')

        self.assertTrue('start_codon' in mrna)

    def test_create_starts_and_stops_creates_a_stop(self):
        mrna = Mrna()
        cds = Mock()
        cds.start = 1
        cds.end = 6
        cds.strand = '+'
        mrna.children = {'cds':[cds]}
        mrna.attributes = {'ID':'m.1234'}
        
        mrna.create_starts_and_stops('NNNTAG')

        self.assertTrue('stop_codon' in mrna)

    def test_create_starts_and_stops_creates_a_stop_reverse_complement(self):
        mrna = Mrna()
        cds = Mock()
        cds.start = 1
        cds.end = 6
        cds.strand = '-'
        mrna.children = {'cds':[cds]}
        mrna.attributes = {'ID':'m.1234'}
        
        mrna.create_starts_and_stops('CTANNN')

        self.assertTrue('stop_codon' in mrna)


###################
# noinspection PyShadowingNames
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestMrna))
    return suite

if __name__ == '__main__':
    unittest.main()
