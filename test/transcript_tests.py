#!/usr/bin/env python

import unittest
from mock import Mock
from src.transcript import gene_to_tbl, Transcript 

class TestTranscript(unittest.TestCase):

    def create_fake_gene(self):
        gff_gene0 = Mock()
        gff_gene0.start = 1
        gff_gene0.end = 100
        gff_gene0.attributes = {"ID":"foo_gene"}
        
        gff_mrna0 = Mock()
        gff_mrna0.start = 1
        gff_mrna0.end = 100
        gff_mrna0.attributes = {"ID":"m.foo"}

        gff_cds0 = Mock()
        gff_cds0.start = 1 
        gff_cds0.end = 100
        
        gff_gene0.mrna = [gff_mrna0]
        gff_mrna0.cds = [gff_cds0]

        del gff_mrna0.start_codon
        del gff_mrna0.stop_codon

        return gff_gene0

    def test_gene_to_tbl_nostart_nostop(self):
        expected = \
        "<1\t>100\tgene\n"\
        "\t\t\tlocus_tag\tfoo_gene\n"\
        "<1\t>100\tCDS\n"\
        "\t\t\tprotein_id\tm.foo\n"

        gff_gene0 = self.create_fake_gene()

        tbl = gene_to_tbl(gff_gene0)
        self.assertEquals(tbl, expected)

    def test_gene_to_tbl_start_nostop(self):
        expected = \
        "1\t>100\tgene\n"\
        "\t\t\tlocus_tag\tfoo_gene\n"\
        "1\t>100\tCDS\n"\
        "\t\t\tprotein_id\tm.foo\n"

        gff_gene0 = self.create_fake_gene()

        gff_gene0.mrna[0].start_codon = [Mock()]

        tbl = gene_to_tbl(gff_gene0)
        self.assertEquals(tbl, expected)

    def test_gene_to_tbl_nostart_stop(self):
        expected = \
        "<1\t100\tgene\n"\
        "\t\t\tlocus_tag\tfoo_gene\n"\
        "<1\t100\tCDS\n"\
        "\t\t\tprotein_id\tm.foo\n"

        gff_gene0 = self.create_fake_gene()

        gff_gene0.mrna[0].stop_codon = [Mock()]

        tbl = gene_to_tbl(gff_gene0)
        self.assertEquals(tbl, expected)

    def test_gene_to_tbl_start_stop(self):
        expected = \
        "1\t100\tgene\n"\
        "\t\t\tlocus_tag\tfoo_gene\n"\
        "1\t100\tCDS\n"\
        "\t\t\tprotein_id\tm.foo\n"

        gff_gene0 = self.create_fake_gene()

        gff_gene0.mrna[0].start_codon = [Mock()]
        gff_gene0.mrna[0].stop_codon = [Mock()]

        tbl = gene_to_tbl(gff_gene0)
        self.assertEquals(tbl, expected)

    def test_gene_to_tbl_genename(self):
        expected = \
        "<1\t>100\tgene\n"\
        "\t\t\tgene\tf00x4\n"\
        "\t\t\tlocus_tag\tfoo_gene\n"\
        "<1\t>100\tCDS\n"\
        "\t\t\tprotein_id\tm.foo\n"

        gff_gene0 = self.create_fake_gene()

        gff_gene0.attributes["Name"] = "f00x4"

        tbl = gene_to_tbl(gff_gene0)
        self.assertEquals(tbl, expected)

    def test_gene_to_tbl_genename(self):
        expected = \
        "<1\t>100\tgene\n"\
        "\t\t\tlocus_tag\tfoo_gene\n"\
        "<1\t>100\tCDS\n"\
        "\t\t\tprotein_id\tm.foo\n"\
        "\t\t\tdb_xref\tPfam:foo\n"\
        "\t\t\tdb_xref\tPfam:dog\n"\
        "\t\t\tdb_xref\tPfam:baz\n"

        gff_gene0 = self.create_fake_gene()

        gff_gene0.mrna[0].attributes["Dbxref"] = "Pfam:foo,Pfam:dog,Pfam:baz"

        tbl = gene_to_tbl(gff_gene0)
        self.assertEquals(tbl, expected)

    def test_to_tbl(self):
        expected = \
        ">Feature foo_seq\n"\
        "1\t4\tREFERENCE\n"\
        "\t\t\tNCBI\t12345\n"\
        "<1\t>100\tgene\n"\
        "\t\t\tlocus_tag\tfoo_gene\n"\
        "<1\t>100\tCDS\n"\
        "\t\t\tprotein_id\tm.foo\n"

        gff_gene0 = self.create_fake_gene()

        sequence = Mock()
        sequence.header = "foo_seq"
        sequence.bases = "ATGC"

        transcript = Transcript([gff_gene0], sequence) 

        tbl = transcript.to_tbl()
        self.assertEquals(tbl, expected)


##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestTranscript))
    return suite

if __name__ == '__main__':
    unittest.main()
