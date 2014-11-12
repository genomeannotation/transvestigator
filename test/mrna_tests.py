import io
import unittest
from unittest.mock import Mock, patch, PropertyMock
from src.mrna import Mrna

class TestMrna(unittest.TestCase):

    def test_from_gff_feature_success(self):
        gff_mrna = Mock()
        gff_mrna.type = "mRNA"
        
        tran_mrna = Mrna.from_gff_feature(gff_mrna)
        self.assertTrue(tran_mrna)

    def test_from_gff_features_fails(self):
        gff_mrna = Mock()
        gff_mrna.type = "asdf"
        
        tran_mrna = Mrna.from_gff_feature(gff_mrna)
        self.assertFalse(tran_mrna)


###################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestMrna))
    return suite

if __name__ == '__main__':
    unittest.main()
