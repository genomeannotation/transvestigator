#!/usr/bin/env python

import unittest
import io
from src.gff_reader import *

class TestGFFReader(unittest.TestCase):

    def setUp(self):
        self.reader = GFFReader()
    
    def test_parse_gff_attributes(self):
        attr = "\t; foo=dog;baz=bub;  \t\n"
        self.assertEquals(parse_gff_attributes(attr), {"foo":"dog", "baz":"bub"})

##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestGFFReader))
    return suite

if __name__ == '__main__':
    unittest.main()
