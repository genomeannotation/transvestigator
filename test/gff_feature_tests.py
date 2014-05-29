import unittest
import io
import os
from mock import Mock, patch, PropertyMock
from src.gff_feature import *

class TestGFFFeature(unittest.TestCase):

    def test_add_child(self):
        root = GFFFeature()
        gene0 = GFFFeature()
        gene0.type = "gene"
        gene1 = GFFFeature()
        gene1.type = "gene"
        mrna = GFFFeature()
        mrna.type = "mRNA"

        root.add_child(gene0)
        root.add_child(gene1)
        gene0.add_child(mrna)

        self.assertEqual(root.gene, [gene0, gene1])
        self.assertEqual(root.gene[0].mrna, [mrna])

    def test_add_annotation_with_no_annotations(self):
        mrna = GFFFeature()
        annokey = "Dbxref"
        annovalue = "Pfam:foo"
        mrna.add_annotation(annokey, annovalue)
        self.assertTrue(mrna.attributes["Dbxref"] == "Pfam:foo")

    def test_add_annotation_with_annotations_already_in_place(self):
        mrna = GFFFeature()
        mrna.attributes["Dbxref"] = "Pfam:bar"
        annokey = "Dbxref"
        annovalue = "Pfam:foo"
        mrna.add_annotation(annokey, annovalue)
        self.assertEquals(mrna.attributes["Dbxref"], "Pfam:bar,Pfam:foo")

    def test_length(self):
        feature = GFFFeature(start=5, end=10)
        self.assertEquals(feature.length(), 6)


###################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestGFFFeature))
    return suite

if __name__ == '__main__':
    unittest.main()
