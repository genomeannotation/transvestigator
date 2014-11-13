#!/usr/bin/env python

import unittest
from unittest.mock import Mock
from src.sequence import Sequence
from src.transcript import gene_to_tbl, Transcript 

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
        self.assertEquals(self.transcript.genes, [self.gene1])


    ### TO TBL TESTS ###

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

        pass

    def test_remove_contig_from_gene_id(self):
        pass

    #### MAKE POSITIVE TESTS ####

    def test_make_positive(self):
        pass

    #### FIX LENGTHS TESTS ####

    def test_fix_feature_lengths(self):
        pass

    #### MATCH CDS AND EXON END TESTS ####

    def test_match_cds_and_exon_end(self):
        pass

    def test_match_cds_and_exon_end_does_nothing_if_stop_codon_present(self):
        pass


    #### STARTS AND STOPS TESTS ####
    
    def test_create_starts_and_stops_creates_a_start(self):
        pass

    def test_create_starts_and_stops_creates_a_start_reverse_complement(self):
        pass

    def test_create_starts_and_stops_creates_a_stop(self):
        pass

    def test_create_starts_and_stops_creates_a_stop_reverse_complement(self):
        pass

##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestTranscript))
    return suite

if __name__ == '__main__':
    unittest.main()
