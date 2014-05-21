import unittest
import io
import os
from mock import Mock, patch, PropertyMock
from src.gff_reader import *

class TestGFFFeature(unittest.TestCase):

    def test_add_child(self):
        root = GFFFeature()
        gene0 = GFFFeature()
        gene0.type = 'gene'
        gene1 = GFFFeature()
        gene1.type = 'gene'
        mrna = GFFFeature()
        mrna.type = 'mRNA'

        root.add_child(gene0)
        root.add_child(gene1)
        gene0.add_child(mrna)

        self.assertEqual(root.gene, [gene0, gene1])
        self.assertEqual(root.gene[0].mrna, [mrna])

    def test_has_child(self):
        feature = GFFFeature()
        foo = GFFFeature()
        foo.type = 'foo'
        self.assertFalse(feature.has_child('foo'))
        feature.add_child(foo)
        self.assertTrue(feature.has_child('foo'))
        

###################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestGFFFeature))
    return suite

if __name__ == '__main__':
    unittest.main()
