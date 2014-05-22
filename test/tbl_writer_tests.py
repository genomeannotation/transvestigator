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

        self.gff_gene1 = Mock()
        self.gff_gene1.start = 42
        self.gff_gene1.end = 80
        self.gff_gene1.attributes = {"ID":"foo_gene"}
        
        self.gff_mrna1 = Mock()
        self.gff_mrna1.start = 42 
        self.gff_mrna1.end = 80
        self.gff_mrna1.attributes = {"ID":"m.foo"}

        self.gff_cds1 = Mock()
        self.gff_cds1.start = 42 
        self.gff_cds1.end = 80
        
        self.gff_gene1.mrna = [self.gff_mrna1]
        self.gff_mrna1.cds = [self.gff_cds1]

    def test_gff_gene_to_tbl_throws_on_no_mrnas(self):
        gff_gene = Mock()
        gff_gene.attributes = {"ID":"foo_gene"}
        del gff_gene.mrna

        thrown = False
        try:
            gff_gene_to_tbl(gff_gene)
        except Exception as error:
            self.assertEquals(str(error), "can't write tbl entry for foo_gene because it has no mRNA")
            thrown = True
        self.assertTrue(thrown)

    def test_gff_gene_to_tbl_throws_on_multiple_mrnas(self):
        gff_gene = Mock()
        gff_gene.attributes = {"ID":"foo_gene"}
        gff_gene.mrna = [Mock(), Mock()]

        thrown = False
        try:
            gff_gene_to_tbl(gff_gene)
        except Exception as error:
            self.assertEquals(str(error), "can't write tbl entry for foo_gene because it has multiple mRNAs")
            thrown = True
        self.assertTrue(thrown)

    def test_gff_gene_to_tbl_throws_on_no_cds(self):
        gff_gene = Mock()
        gff_gene.attributes = {"ID":"foo_gene"}
        gff_gene.mrna = [Mock()]
        del gff_gene.mrna[0].cds

        thrown = False
        try:
            gff_gene_to_tbl(gff_gene)
        except Exception as error:
            self.assertEquals(str(error), "can't write tbl entry for foo_gene because its mRNA has no CDS")
            thrown = True
        self.assertTrue(thrown)

    def test_gff_gene_to_tbl_nostart_nostop_nogenename(self):
        expected = \
        "<1\t>100\tgene\n"\
        "\t\t\tlocus_tag\tfoo_gene\n"\
        "<1\t>100\tCDS\n"\
        "\t\t\tprotein_id\tm.foo\n"

        del self.gff_mrna0.start_codon
        del self.gff_mrna0.stop_codon

        tbl = gff_gene_to_tbl(self.gff_gene0)
        self.assertEquals(tbl, expected)

    def test_gff_gene_to_tbl_start_nostop_nogenename(self):
        expected = \
        "1\t>100\tgene\n"\
        "\t\t\tlocus_tag\tfoo_gene\n"\
        "1\t>100\tCDS\n"\
        "\t\t\tprotein_id\tm.foo\n"

        self.gff_mrna0.start_codon = [Mock()]
        del self.gff_mrna0.stop_codon

        tbl = gff_gene_to_tbl(self.gff_gene0)
        self.assertEquals(tbl, expected)

    def test_gff_gene_to_tbl_nostart_stop_nogenename(self):
        expected = \
        "<1\t100\tgene\n"\
        "\t\t\tlocus_tag\tfoo_gene\n"\
        "<1\t100\tCDS\n"\
        "\t\t\tprotein_id\tm.foo\n"

        del self.gff_mrna0.start_codon
        self.gff_mrna0.stop_codon = [Mock()]

        tbl = gff_gene_to_tbl(self.gff_gene0)
        self.assertEquals(tbl, expected)

    def test_gff_gene_to_tbl_start_stop_nogenename(self):
        expected = \
        "1\t100\tgene\n"\
        "\t\t\tlocus_tag\tfoo_gene\n"\
        "1\t100\tCDS\n"\
        "\t\t\tprotein_id\tm.foo\n"

        self.gff_mrna0.start_codon = [Mock()]
        self.gff_mrna0.stop_codon = [Mock()]

        tbl = gff_gene_to_tbl(self.gff_gene0)
        self.assertEquals(tbl, expected)

    def test_transcript_to_tbl_throws_on_no_genes(self):
        transcript = Mock()
        transcript.genes = [] 
        transcript.sequence.header = "foo_seq"

        thrown = False
        try:
            tbl = transcript_to_tbl(transcript)
        except Exception as error:
            self.assertEquals(str(error), "can't write tbl entry for foo_seq because it has no genes")
            thrown = True
        self.assertTrue(thrown)

    def test_transcript_to_tbl_throws_on_multiple_genes(self):
        transcript = Mock()
        transcript.genes = [self.gff_gene0, self.gff_gene1] 
        transcript.sequence.header = "foo_seq"

        thrown = False
        try:
            tbl = transcript_to_tbl(transcript)
        except Exception as error:
            self.assertEquals(str(error), "can't write tbl entry for foo_seq because it has multiple genes")
            thrown = True
        self.assertTrue(thrown)

    def test_transcript_to_tbl(self):
        expected = \
        ">Feature foo_seq\n"\
        "1\t4\tREFERENCE\n"\
        "\t\t\tNCBI\t12345\n"\
        "1\t100\tgene\n"\
        "\t\t\tlocus_tag\tfoo_gene\n"\
        "1\t100\tCDS\n"\
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
