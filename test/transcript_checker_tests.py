#!/usr/bin/env python

import unittest
from mock import Mock
from src.transcript_checker import TranscriptChecker

class TestTranscriptChecker(unittest.TestCase):

    def setUp(self):
        self.checker = TranscriptChecker()

    def test_overlap(self):
        pass


##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestTranscriptChecker))
    return suite

if __name__ == '__main__':
    unittest.main()
