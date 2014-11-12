import io
import unittest
from unittest.mock import Mock, patch, PropertyMock
from src.gene import Gene

class TestGene(unittest.TestCase):

    def test_from_gff_feature_success(self):
        gff_gene = Mock()
        gff_gene.type = "gene"
        
        tran_gene = Gene.from_gff_feature(gff_gene)
        self.assertTrue(tran_gene)

    def test_from_gff_features_fails(self):
        gff_gene = Mock()
        gff_gene.type = "asdf"
        
        tran_gene = Gene.from_gff_feature(gff_gene)
        self.assertFalse(tran_gene)


###################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestGene))
    return suite

if __name__ == '__main__':
    unittest.main()
