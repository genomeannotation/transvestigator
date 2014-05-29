#!/usr/bin/env python

import unittest
import io
from src.sequence import Sequence

class TestSequence(unittest.TestCase):

    def test_to_fasta(self):
        sequence = Sequence("foo_seq", "GATTACA")
        fasta = sequence.to_fasta()
        expected = ">foo_seq\nGATTACA\n"
        self.assertEquals(expected, fasta)

