#!/usr/bin/env python

import unittest
from src.sequtil import reverse_complement, has_start_codon, has_stop_codon, get_subsequence, overlap, nested

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

    #### OVERLAP TESTS ####

    def test_overlap(self):
        indices1 = [1, 10]
        indices2 = [9, 15]
        self.assertTrue(overlap(indices1, indices2))

    def test_overlap_returns_false_when_no_overlap(self):
        indices1 = [1, 10]
        indices2 = [11, 20]
        self.assertFalse(overlap(indices1, indices2))

    def test_overlap_returns_false_when_nested(self):
        indices1 = [1, 10]
        indices2 = [2, 9]
        self.assertFalse(overlap(indices1, indices2))

    def test_overlap_returns_false_when_identical(self): # not sure if this really happens
        indices1 = [1, 10]
        indices2 = [1, 10]
        self.assertFalse(overlap(indices1, indices2))

    def test_overlap_returns_true_when_sharing_one_base(self):
        indices1 = [1, 5]
        indices2 = [5, 10]
        self.assertTrue(overlap(indices1,indices2))
        
    def test_overlap_returns_true_when_sharing_one_base_case2(self):
        indices1 = [5, 10]
        indices2 = [1, 5]
        self.assertTrue(overlap(indices1,indices2))

    def test_list_of_index_pairs_contains_overlap(self):
        pass


    #### NESTED TESTS ####

    def test_nested(self):
        indices1 = [1, 10]
        indices2 = [2, 9]
        self.assertTrue(nested(indices1, indices2))

    def test_nested_returns_false_when_overlapping(self):
        indices1 = [1, 10]
        indices2 = [9, 15]
        self.assertFalse(nested(indices1, indices2))

    def test_nested_returns_true_when_identical(self):
        indices1 = [1, 10]
        indices2 = [1, 10]
        self.assertTrue(nested(indices1, indices2))

    def test_nested_returns_false_when_separate(self):
        indices1 = [1, 10]
        indices2 = [11, 20]
        self.assertFalse(nested(indices1, indices2))

    def test_nested_returns_true_when_flushed_on_one_side(self):
        indices1 = [1, 10]
        indices2 = [1, 5]
        self.assertTrue(nested(indices1, indices2))

        
##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestTranslate))
    return suite

if __name__ == '__main__':
    unittest.main()
