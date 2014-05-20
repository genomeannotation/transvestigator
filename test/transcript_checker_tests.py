#!/usr/bin/env python

import unittest
from mock import Mock
from src.transcript_checker import TranscriptChecker

class TestTranscriptChecker(unittest.TestCase):

    def setUp(self):
        self.checker = TranscriptChecker()

    #### OVERLAP TESTS ####

    def test_overlap(self):
        indices1 = [1, 10]
        indices2 = [9, 15]
        self.assertTrue(self.checker.overlap(indices1, indices2))

    def test_overlap_returns_false_when_no_overlap(self):
        indices1 = [1, 10]
        indices2 = [11, 20]
        self.assertFalse(self.checker.overlap(indices1, indices2))

    def test_overlap_returns_false_when_nested(self):
        indices1 = [1, 10]
        indices2 = [2, 9]
        self.assertFalse(self.checker.overlap(indices1, indices2))

    def test_overlap_returns_false_when_identical(self): # not sure if this really happens
        indices1 = [1, 10]
        indices2 = [1, 10]
        self.assertFalse(self.checker.overlap(indices1, indices2))

    def test_overlap_returns_true_when_sharing_one_base(self):
        indices1 = [1, 5]
        indices2 = [5, 10]
        self.assertTrue(self.checker.overlap(indices1,indices2))
        
    def test_overlap_returns_true_when_sharing_one_base_case2(self):
        indices1 = [5, 10]
        indices2 = [1, 5]
        self.assertTrue(self.checker.overlap(indices1,indices2))

    def test_list_of_index_pairs_contains_overlap(self):
        pass


    #### NESTED TESTS ####

    def test_nested(self):
        indices1 = [1, 10]
        indices2 = [2, 9]
        self.assertTrue(self.checker.nested(indices1, indices2))

    def test_nested_returns_false_when_overlapping(self):
        indices1 = [1, 10]
        indices2 = [9, 15]
        self.assertFalse(self.checker.nested(indices1, indices2))

    def test_nested_returns_true_when_identical(self):
        indices1 = [1, 10]
        indices2 = [1, 10]
        self.assertTrue(self.checker.nested(indices1, indices2))

    def test_nested_returns_false_when_separate(self):
        indices1 = [1, 10]
        indices2 = [11, 20]
        self.assertFalse(self.checker.nested(indices1, indices2))

    def test_nested_returns_true_when_flushed_on_one_side(self):
        indices1 = [1, 10]
        indices2 = [1, 5]
        self.assertTrue(self.checker.nested(indices1, indices2))


##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestTranscriptChecker))
    return suite

if __name__ == '__main__':
    unittest.main()
