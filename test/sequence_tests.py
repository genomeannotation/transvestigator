#!/usr/bin/env python

import unittest
from mock import Mock
from src.sequence import Sequence

class TestSequence(unittest.TestCase):

    def setUp(self):
        self.seq1 = Sequence("seq1", "GATTACA")

    def test_get_subseq(self):
        self.assertEquals("ATTA", self.seq1.get_subseq(2, 5))

    def test_to_fasta(self):
        expected = ">seq1\nGATTACA\n"
        self.assertEquals(expected, self.seq1.to_fasta())

    def test_to_boulder_io(self):
        expected = "SEQUENCE_ID=seq1\nSEQUENCE_TEMPLATE=GATTACA\n=\n"
        self.assertEquals(expected, self.seq1.to_boulder_io())


##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestSequence))
    return suite

if __name__ == '__main__':
    unittest.main()
