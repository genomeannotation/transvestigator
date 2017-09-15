#!/usr/bin/env python

import unittest
from unittest.mock import Mock
from src.sequence import Sequence
from src.transcript import Transcript 

class TestTranscript(unittest.TestCase):

    def setUp(self):
        self.transcript = Transcript()
        self.transcript.sequence = Mock()
        self.transcript.sequence.bases = "GATACA"
        self.gene0 = Mock()
        self.gene1 = Mock()
        self.transcript.genes = [self.gene0, self.gene1]

    def test_fix_multiple_genes(self):
        self.gene0.get_cds_length.return_value = 10
        self.gene1.get_cds_length.return_value = 11
        self.transcript.fix_multiple_genes()
        self.assertEquals([self.gene1], self.transcript.genes)

    def test_to_tbl(self):
        self.gene0.to_tbl.return_value = \
        "<1\t>100\tgene\n"\
        "\t\t\tlocus_tag\tfoo_gene\n"\
        "<1\t>100\tCDS\n"\
        "\t\t\tprotein_id\tm.foo\n"\
        "\t\t\tproduct\thypothetical protein\n"
        self.transcript.sequence.header = "foo_seq"
        expected = \
        ">Feature foo_seq\n"\
        "1\t6\tREFERENCE\n"\
        "<1\t>100\tgene\n"\
        "\t\t\tlocus_tag\tfoo_gene\n"\
        "<1\t>100\tCDS\n"\
        "\t\t\tprotein_id\tm.foo\n"\
        "\t\t\tproduct\thypothetical protein\n"
        actual = self.transcript.to_tbl()
        self.assertEquals(expected, actual)


##########################
# noinspection PyShadowingNames
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestTranscript))
    return suite

if __name__ == '__main__':
    unittest.main()
