#!/usr/bin/env python

import unittest
from src.annotation import Annotation, validate_annotation

class TestAnnotation(unittest.TestCase):

    def test_validate_annotation(self):
        bad_anno = Annotation("foo_id", "invalid key!", "foo_value")
        self.assertFalse(validate_annotation(bad_anno))



##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestAnnotation))
    return suite

if __name__ == '__main__':
    unittest.main()
