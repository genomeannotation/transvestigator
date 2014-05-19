#!/usr/bin/env python

import unittest
from src.translator import reverse_complement

class TestTranslate(unittest.TestCase):

    def test_reverse_complement(self):
        self.assertEquals('C', reverse_complement('G'))
        self.assertEquals('CAT', reverse_complement('ATG'))

    def test_reverse_complement_with_bogus_base(self):
        self.assertEquals('CATN', reverse_complement('MATG'))

    def test_reverse_complement_longer_seq(self):
        self.assertEquals('TGTAATCTGTAATCTGTAATCTGTAATCTGTAATC', reverse_complement('GATTACAGATTACAGATTACAGATTACAGATTACA'))

        
##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestTranslate))
    return suite

if __name__ == '__main__':
    unittest.main()
