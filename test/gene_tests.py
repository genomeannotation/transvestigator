import io
import unittest
from unittest.mock import Mock, patch, PropertyMock
from src.gene import Gene

class TestGene(unittest.TestCase):

    def setUp(self):
        self.gene1 = Gene()

    def test_from_gff_feature_success(self):
        self.gene1 = Mock()
        self.gene1.type = "gene"
        
        tran_gene = Gene.from_gff_feature(self.gene1)
        self.assertTrue(tran_gene)

    def test_from_gff_features_fails(self):
        gff_gene = Mock()
        gff_gene.type = "asdf"
        
        tran_gene = Gene.from_gff_feature(gff_gene)
        self.assertFalse(tran_gene)

    def test_remove_contig_from_gene_id(self):
        expected = 'g.123'
        self.gene1.attributes['ID'] = 'contig123|g.123'
        self.gene1.remove_contig_from_gene_id()
        self.assertEquals(expected, self.gene1.attributes['ID'])

"""
    def test_gene_to_tbl_nostart_nostop(self):
        expected = \
        "<1\t>100\tgene\n"\
        "\t\t\tlocus_tag\tfoo_gene\n"\
        "<1\t>100\tCDS\n"\
        "\t\t\tprotein_id\tm.foo\n"\
        "\t\t\tproduct\thypothetical protein\n"

        #gff_gene0 = self.create_fake_gene()

        tbl = gene_to_tbl(gff_gene0)
    #    self.assertEquals(tbl, expected)

    def test_gene_to_tbl_start_nostop(self):
        expected = \
        "1\t>100\tgene\n"\
        "\t\t\tlocus_tag\tfoo_gene\n"\
        "1\t>100\tCDS\n"\
        "\t\t\tprotein_id\tm.foo\n"\
        "\t\t\tproduct\thypothetical protein\n"

        #gff_gene0 = self.create_fake_gene()

        gff_gene0["mrna"][0].start_codon = [Mock()]

        tbl = gene_to_tbl(gff_gene0)
        #self.assertEquals(tbl, expected)

    def test_gene_to_tbl_nostart_stop(self):
        expected = \
        "<1\t100\tgene\n"\
        "\t\t\tlocus_tag\tfoo_gene\n"\
        "<1\t100\tCDS\n"\
        "\t\t\tprotein_id\tm.foo\n"\
        "\t\t\tproduct\thypothetical protein\n"

        #gff_gene0 = self.create_fake_gene()

        gff_gene0["mrna"][0].stop_codon = [Mock()]

        tbl = gene_to_tbl(gff_gene0)
        #self.assertEquals(tbl, expected)

    def test_gene_to_tbl_start_stop(self):
        expected = \
        "1\t100\tgene\n"\
        "\t\t\tlocus_tag\tfoo_gene\n"\
        "1\t100\tCDS\n"\
        "\t\t\tprotein_id\tm.foo\n"\
        "\t\t\tproduct\thypothetical protein\n"

        #gff_gene0 = self.create_fake_gene()

        gff_gene0["mrna"][0].start_codon = [Mock()]
        gff_gene0["mrna"][0].stop_codon = [Mock()]

        tbl = gene_to_tbl(gff_gene0)
        #self.assertEquals(tbl, expected)

    def test_gene_to_tbl_genename(self):
        expected = \
        "<1\t>100\tgene\n"\
        "\t\t\tgene\tf00x4\n"\
        "\t\t\tlocus_tag\tfoo_gene\n"\
        "<1\t>100\tCDS\n"\
        "\t\t\tprotein_id\tm.foo\n"\
        "\t\t\tproduct\thypothetical protein\n"

        #gff_gene0 = self.create_fake_gene()

        gff_gene0.attributes["Name"] = "f00x4"

        tbl = gene_to_tbl(gff_gene0)
        #self.assertEquals(tbl, expected)

    def test_gene_to_tbl_dbxref(self):
        expected = \
        "<1\t>100\tgene\n"\
        "\t\t\tlocus_tag\tfoo_gene\n"\
        "<1\t>100\tCDS\n"\
        "\t\t\tprotein_id\tm.foo\n"\
        "\t\t\tdb_xref\tPfam:foo\n"\
        "\t\t\tdb_xref\tPfam:dog\n"\
        "\t\t\tdb_xref\tPfam:baz\n"\
        "\t\t\tproduct\thypothetical protein\n"\

        #gff_gene0 = self.create_fake_gene()

        gff_gene0["mrna"][0].attributes["Dbxref"] = "Pfam:foo,Pfam:dog,Pfam:baz"

        tbl = gene_to_tbl(gff_gene0)
        #self.assertEquals(tbl, expected)

    def test_fix_phase(self):
        self.gene0.start = 2
        self.mrna0.start = 2
        del self.mrna0.start_codon
        self.cds0.start = 2
        self.cds0.phase = 0
        #self.assertEqual(self.cds0.phase, 0)
        self.transcript.fix_phase()
        #self.assertEqual(self.cds0.phase, 1)
        
    def test_fix_phase_to_two(self):
        self.gene0.start = 3
        self.mrna0.start = 3
        del self.mrna0.start_codon
        self.cds0.start = 3
        self.cds0.phase = 0
        #self.assertEqual(self.cds0.phase, 0)
        self.transcript.fix_phase()
        #self.assertEqual(self.cds0.phase, 2)
        
    def test_fix_phase_does_nothing_when_indices_too_large(self):
        self.gene0.start = 4
        self.mrna0.start = 4
        self.cds0.start = 4
        self.cds0.phase = 0
        #self.assertEqual(self.cds0.phase, 0)
        self.transcript.fix_phase()
        #self.assertEqual(self.cds0.phase, 0)

    def test_fix_phase_works_on_cds_only(self):
        self.gene0.start = 1
        self.mrna0.start = 1
        self.cds0.start = 3
        self.cds0.phase = 0
        del self.mrna0.start_codon
        #self.assertEqual(self.cds0.phase, 0)
        #self.assertEqual(self.cds0.start, 3)
        self.transcript.fix_phase()
        #self.assertEqual(self.cds0.phase, 2)
        #self.assertEqual(self.cds0.start, 1)

    def test_fix_phase_does_nothing_when_not_partial(self):
        self.gene0.start = 2
        self.mrna0.start = 2
        self.cds0.start = 2
        self.cds0.phase = 0
        #self.assertEqual(self.cds0.phase, 0)
        self.transcript.fix_phase()
        #self.assertEqual(self.cds0.phase, 0)

    def test_fix_phase_adjusts_end_on_3prime_partial(self):
        self.gene0.start = 2
        self.mrna0.start = 2
        self.cds0.start = 2
        self.gene0.end = 2
        self.mrna0.end = 2
        self.cds0.end = 2
        del self.mrna0.stop_codon
        #self.assertEqual(self.cds0.end, 2)
        self.transcript.fix_phase()
        #self.assertEqual(self.cds0.end, 6)

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

"""

###################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestGene))
    return suite

if __name__ == '__main__':
    unittest.main()
