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
        cds.type = 'cds'
        cds.start = 1
        cds.end = 7
        cds.strand = '-'
        exon = Mock()
        exon.type = 'exon'
        exon.start = 1
        exon.end = 7
        cds.strand = '-'

        mrna.add_child(cds)
        mrna.add_child(exon)

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


###################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestMrna))
    return suite

if __name__ == '__main__':
    unittest.main()
