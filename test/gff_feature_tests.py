# coding=utf-8
import unittest
from unittest.mock import Mock

from src.gff_feature import *


class TestGFFFeature(unittest.TestCase):
    def test_add_child(self):
        root = GFFFeature()
        gene0 = GFFFeature()
        gene0.feature_type = "gene"
        gene1 = GFFFeature()
        gene1.feature_type = "gene"
        mrna = GFFFeature()
        mrna.feature_type = "mRNA"

        root.add_child(gene0)
        root.add_child(gene1)
        gene0.add_child(mrna)

        self.assertEqual(root["gene"], [gene0, gene1])
        self.assertEqual(root["gene"][0]["mrna"], [mrna])

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

    def test_fix_feature_lengths_length_ok(self):
        feature = GFFFeature(start=1, end=7)

        seq_len = 8
        feature.fix_feature_lengths(seq_len)

        self.assertEquals(feature.start, 1)
        self.assertEquals(feature.end, 7)

    def test_fix_feature_lengths_calls_child_features(self):
        feature = GFFFeature(start=1, end=7)
        child = Mock()
        feature.children = {'child': [child]}

        seq_len = 8
        feature.fix_feature_lengths(seq_len)

        child.fix_feature_lengths.assertCalledWith(seq_len)

    def test_fix_feature_lengths_1(self):
        feature = GFFFeature(start=1, end=9)

        seq_len = 8
        feature.fix_feature_lengths(seq_len)

        self.assertEquals(feature.start, 1)
        self.assertEquals(feature.end, 6)

    def test_fix_feature_lengths_2(self):
        feature = GFFFeature(start=2, end=10)

        seq_len = 8
        feature.fix_feature_lengths(seq_len)

        self.assertEquals(feature.start, 2)
        self.assertEquals(feature.end, 7)

    def test_fix_feature_lengths_3(self):
        feature = GFFFeature(start=3, end=11)

        seq_len = 8
        feature.fix_feature_lengths(seq_len)

        self.assertEquals(feature.start, 3)
        self.assertEquals(feature.end, 8)

    def test_write(self):
        test = GFFFeature("foo_seq", "geiblabs", "foo", 1, 42, None, '-', 1, {"ID": "1234", "dog": "bazz"})
        expected = "foo_seq\tgeiblabs\tfoo\t1\t42\t.\t-\t1\tID=1234;dog=bazz"
        self.assertEqual(test.write(), expected)


###################
# noinspection PyShadowingNames
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestGFFFeature))
    return suite


if __name__ == '__main__':
    unittest.main()
