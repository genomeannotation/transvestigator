#!/usr/bin/env python

import unittest
from src.annotation import Annotation

class TestAnnotation(unittest.TestCase):

    def test_validate_annotation(self):
        pass



##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestAnnotation))
    return suite

if __name__ == '__main__':
    unittest.main()
