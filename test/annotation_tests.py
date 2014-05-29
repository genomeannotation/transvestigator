#!/usr/bin/env python

import unittest
import io
from src.annotation import Annotation, validate_annotation, read_annotations

class TestAnnotation(unittest.TestCase):

    def test_validate_annotation(self):
        bad_anno = Annotation("foo_id", "invalid key!", "foo_value")
        self.assertFalse(validate_annotation(bad_anno))

    def test_read_annotations(self):
        anno_input = io.StringIO("foo_mrna\tproduct\tfoo-like protein\n"
                                 "foo_gene\tname\tfgnX\n")
        anno1 = Annotation("foo_mrna", "product", "foo-like protein")
        anno2 = Annotation("foo_gene", "name", "fgnX")
        expected = [anno1, anno2]
        actual = read_annotations(anno_input)
        self.assertEqual(actual, expected)




##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestAnnotation))
    return suite

if __name__ == '__main__':
    unittest.main()
