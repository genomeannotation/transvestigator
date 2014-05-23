#!/usr/bin/env python

import unittest
from mock import Mock
from src.tbl_writer import gff_gene_to_tbl, transcript_to_tbl

class TestTblWriter(unittest.TestCase):

    def setUp(self):
        self.gff_gene0 = Mock()
        self.gff_gene0.start = 1
        self.gff_gene0.end = 100
        self.gff_gene0.attributes = {"ID":"foo_gene"}
        
        self.gff_mrna0 = Mock()
        self.gff_mrna0.start = 1
        self.gff_mrna0.end = 100
        self.gff_mrna0.attributes = {"ID":"m.foo"}

        self.gff_cds0 = Mock()
        self.gff_cds0.start = 1 
        self.gff_cds0.end = 100
        
        self.gff_gene0.mrna = [self.gff_mrna0]
        self.gff_mrna0.cds = [self.gff_cds0]

        del self.gff_mrna0.start_codon
        del self.gff_mrna0.stop_codon

    def test_gff_gene_to_tbl_nostart_nostop(self):
        expected = \
        "<1\t>100\tgene\n"\
        "\t\t\tlocus_tag\tfoo_gene\n"\
        "<1\t>100\tCDS\n"\
        "\t\t\tprotein_id\tm.foo\n"

        tbl = gff_gene_to_tbl(self.gff_gene0)
        self.assertEquals(tbl, expected)

    def test_gff_gene_to_tbl_start_nostop(self):
        expected = \
        "1\t>100\tgene\n"\
        "\t\t\tlocus_tag\tfoo_gene\n"\
        "1\t>100\tCDS\n"\
        "\t\t\tprotein_id\tm.foo\n"

        self.gff_mrna0.start_codon = [Mock()]

        tbl = gff_gene_to_tbl(self.gff_gene0)
        self.assertEquals(tbl, expected)

    def test_gff_gene_to_tbl_nostart_stop(self):
        expected = \
        "<1\t100\tgene\n"\
        "\t\t\tlocus_tag\tfoo_gene\n"\
        "<1\t100\tCDS\n"\
        "\t\t\tprotein_id\tm.foo\n"

        self.gff_mrna0.stop_codon = [Mock()]

        tbl = gff_gene_to_tbl(self.gff_gene0)
        self.assertEquals(tbl, expected)

    def test_gff_gene_to_tbl_start_stop(self):
        expected = \
        "1\t100\tgene\n"\
        "\t\t\tlocus_tag\tfoo_gene\n"\
        "1\t100\tCDS\n"\
        "\t\t\tprotein_id\tm.foo\n"

        self.gff_mrna0.start_codon = [Mock()]
        self.gff_mrna0.stop_codon = [Mock()]

        tbl = gff_gene_to_tbl(self.gff_gene0)
        self.assertEquals(tbl, expected)

    def test_gff_gene_to_tbl_genename(self):
        expected = \
        "<1\t>100\tgene\n"\
        "\t\t\tgene\tf00x4\n"\
        "\t\t\tlocus_tag\tfoo_gene\n"\
        "<1\t>100\tCDS\n"\
        "\t\t\tprotein_id\tm.foo\n"

        self.gff_gene0.attributes["Name"] = "f00x4"

        tbl = gff_gene_to_tbl(self.gff_gene0)
        self.assertEquals(tbl, expected)

    def test_transcript_to_tbl(self):
        expected = \
        ">Feature foo_seq\n"\
        "1\t4\tREFERENCE\n"\
        "\t\t\tNCBI\t12345\n"\
        "<1\t>100\tgene\n"\
        "\t\t\tlocus_tag\tfoo_gene\n"\
        "<1\t>100\tCDS\n"\
        "\t\t\tprotein_id\tm.foo\n"

        transcript = Mock()
        transcript.genes = [self.gff_gene0] 
        transcript.sequence = Mock()
        transcript.sequence.header = "foo_seq"
        transcript.sequence.bases = "ATGC"

        tbl = transcript_to_tbl(transcript)
        self.assertEquals(tbl, expected)


##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestTblWriter))
    return suite

if __name__ == '__main__':
    unittest.main()
