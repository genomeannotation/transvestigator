#!/usr/bin/env python

import unittest
from src.sequtil import reverse_complement, has_start_codon, has_stop_codon, get_subsequence

class TestTranslate(unittest.TestCase):

    def test_reverse_complement(self):
        self.assertEquals('C', reverse_complement('G'))
        self.assertEquals('CAT', reverse_complement('ATG'))

    def test_reverse_complement_with_bogus_base(self):
        self.assertEquals('CATN', reverse_complement('MATG'))

    def test_reverse_complement_longer_seq(self):
        self.assertEquals('TGTAATCTGTAATCTGTAATCTGTAATCTGTAATC', reverse_complement('GATTACAGATTACAGATTACAGATTACAGATTACA'))

    def test_get_subsequence(self):
        bases = "GATTACA"
        subseq = get_subsequence(bases, 2, 4)
        self.assertEquals("ATT", subseq)

    def test_has_start_codon(self):
        self.assertTrue(has_start_codon('auggattaca'))
        self.assertFalse(has_start_codon('guggattaca')) # currently no support for alternate start codons

    def test_has_stop_codon(self):
        self.assertTrue(has_stop_codon('gattacatag'))
        self.assertTrue(has_stop_codon('gattacataa'))
        self.assertTrue(has_stop_codon('gattacatga'))
        self.assertFalse(has_stop_codon('gattacaact'))

        
##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestTranslate))
    return suite

if __name__ == '__main__':
    unittest.main()
