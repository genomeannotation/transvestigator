#!/usr/bin/env python

import unittest
import io
from src.sequence import Sequence, read_fasta

class TestSequence(unittest.TestCase):

    def test_to_fasta(self):
        sequence = Sequence("foo_seq", "GATTACA")
        fasta = sequence.to_fasta()
        expected = ">foo_seq\nGATTACA\n"
        self.assertEquals(expected, fasta)

    def test_read_fasta_with_line_breaks(self):
        line_breaks = io.StringIO('>seq_1\nGATTACAGATTACAGATTACAGATTACA\nGATTACAGATTACAGATTACAGATTACA\n' +
                                 '>seq_2\nNNNNNNNNGATTACAGATTACAGATTAC\nANNNNNNNNNNN')

        seqs = read_fasta(line_breaks)
        self.assertEquals(2, len(seqs))
        self.assertEquals('NNNNNNNNGATTACAGATTACAGATTACANNNNNNNNNNN', seqs["seq_2"].bases)

    def test_read_fasta_without_line_breaks(self):
        no_line_breaks = io.StringIO('>seq_1\nGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACA\n' +
                                    '>seq_2\nNNNNNNNNGATTACAGATTACAGATTACANNNNNNNNNNN')
        seqs = read_fasta(no_line_breaks)
        self.assertEquals(2, len(seqs))
        self.assertEquals('seq_1', seqs["seq_1"].header)
        self.assertEquals('GATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACA', seqs["seq_1"].bases)
        self.assertEquals('seq_2', seqs["seq_2"].header)
        self.assertEquals('NNNNNNNNGATTACAGATTACAGATTACANNNNNNNNNNN', seqs["seq_2"].bases)
        
##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestSequence))
    return suite

if __name__ == '__main__':
    unittest.main()
