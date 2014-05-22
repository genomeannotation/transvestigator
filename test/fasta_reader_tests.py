#!/usr/bin/env python

import unittest
import io
from src.fasta_reader import read_fasta

class TestFastaReader(unittest.TestCase):

    def test_read_with_line_breaks(self):
        line_breaks = io.StringIO('>seq_1\nGATTACAGATTACAGATTACAGATTACA\nGATTACAGATTACAGATTACAGATTACA\n' +
                                 '>seq_2\nNNNNNNNNGATTACAGATTACAGATTAC\nANNNNNNNNNNN')

        seqs = read_fasta(line_breaks)
        self.assertEquals(2, len(seqs))
        self.assertEquals('NNNNNNNNGATTACAGATTACAGATTACANNNNNNNNNNN', seqs[1].bases)

    def test_read_without_line_breaks(self):
        no_line_breaks = io.StringIO('>seq_1\nGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACA\n' +
                                    '>seq_2\nNNNNNNNNGATTACAGATTACAGATTACANNNNNNNNNNN')
        seqs = read_fasta(no_line_breaks)
        self.assertEquals(2, len(seqs))
        self.assertEquals('seq_1', seqs[0].header)
        self.assertEquals('GATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACA', seqs[0].bases)
        self.assertEquals('seq_2', seqs[1].header)
        self.assertEquals('NNNNNNNNGATTACAGATTACAGATTACANNNNNNNNNNN', seqs[1].bases)
        


##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestFastaReader))
    return suite

if __name__ == '__main__':
    unittest.main()
