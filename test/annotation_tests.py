#!/usr/bin/env python
# coding=utf-8

import unittest
import io
from src.annotation import read_annotations, annotate_genes
from src.gff import GFFFeature

class TestAnnotation(unittest.TestCase):

    def test_read_annotations(self):
        anno_input = io.StringIO("foo_mrna\tproduct\tfoo-like protein\n"
                                 "foo_gene\tname\tfgnX\n")
        annos = read_annotations(anno_input)
        self.assertTrue("foo_mrna" in annos)
        self.assertTrue("foo_gene" in annos)
        self.assertEqual(annos["foo_mrna"][0], ("product", "foo-like protein"))

    def test_annotate_genes(self):
        gene = GFFFeature()
        gene['mrna'] = [GFFFeature()]
        gene['mrna'][0].attributes["ID"] = ":)"
        annotations = {":(" : [["DBXREF", "123"]], ":)" : [["DBXREF", "321"]]}
        self.assertTrue("DBXREF" not in gene['mrna'][0].attributes)
        annotate_genes([gene], annotations)
        self.assertTrue("DBXREF" in gene['mrna'][0].attributes)
        self.assertEquals("321", gene['mrna'][0].attributes["DBXREF"])

    def test_annotate_genes_adds_products(self):
        gene = GFFFeature()
        gene['mrna'] = [GFFFeature()]
        gene['mrna'][0].attributes["ID"] = ":)"
        annotations = {":(" : [["DBXREF", "123"]], ":)" : [["product", "foo product"]]}
        self.assertTrue("product" not in gene['mrna'][0].attributes)
        annotate_genes([gene], annotations)
        self.assertTrue("product" in gene['mrna'][0].attributes)
        self.assertEquals("foo product", gene['mrna'][0].attributes["product"])

    def test_annotate_genes_names_genes(self):
        gene = GFFFeature()
        gene.attributes = {"ID": "foo_gene"}
        gene['mrna'] = [GFFFeature()]
        gene['mrna'][0].attributes["ID"] = ":)"
        annotations = {"foo_gene": [("name", "fGnX")]}
        self.assertTrue("DBXREF" not in gene['mrna'][0].attributes)
        annotate_genes([gene], annotations)
        self.assertEqual("fGnX", gene.attributes["Name"])


##########################

# noinspection PyShadowingNames
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestAnnotation))
    return suite

if __name__ == '__main__':
    unittest.main()
