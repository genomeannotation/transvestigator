#!/usr/bin/env python

import unittest
from mock import Mock
from src.tbl_writer import gff_gene_to_tbl

class TestTblWriter(unittest.TestCase):

    def setUp(self):
        pass

    def test_gff_gene_to_tbl_throws_on_no_mrnas(self):
        gff_gene = Mock()
        gff_gene.attributes = {"ID":"foo_gene"}
        del gff_gene.mrna

        thrown = False
        try:
            gff_gene_to_tbl(gff_gene)
        except Exception as error:
            self.assertEquals(str(error), "can't write tbl entry for foo_gene because it has no mRNAs")
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

        gff_gene = Mock()
        gff_gene.start = 1
        gff_gene.end = 100
        gff_gene.attributes = {"ID":"foo_gene"}
        
        gff_mrna = Mock()
        gff_mrna.start = 1
        gff_mrna.end = 100
        gff_mrna.attributes = {"ID":"m.foo"}

        gff_cds = Mock()
        gff_cds.start = 1 
        gff_cds.end = 100
        
        gff_gene.mrna = [gff_mrna]
        gff_mrna.cds = [gff_cds]
        del gff_mrna.start_codon

        tbl = gff_gene_to_tbl(gff_gene)
        self.assertEquals(tbl, expected)

    def test_gff_gene_to_tbl_start_nostop_nogenename(self):
        expected = \
        "1\t>100\tgene\n"\
        "\t\t\tlocus_tag\tfoo_gene\n"\
        "1\t>100\tCDS\n"\
        "\t\t\tprotein_id\tm.foo\n"

        gff_gene = Mock()
        gff_gene.start = 1
        gff_gene.end = 100
        gff_gene.attributes = {"ID":"foo_gene"}
        
        gff_mrna = Mock()
        gff_mrna.start = 1
        gff_mrna.end = 100
        gff_mrna.attributes = {"ID":"m.foo"}

        gff_cds = Mock()
        gff_cds.start = 1 
        gff_cds.end = 100
        
        gff_gene.mrna = [gff_mrna]
        gff_mrna.cds = [gff_cds]
        gff_mrna.start_codon = [Mock()]

        tbl = gff_gene_to_tbl(gff_gene)
        self.assertEquals(tbl, expected)
        

##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestTblWriter))
    return suite

if __name__ == '__main__':
    unittest.main()
