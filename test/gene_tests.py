import io
import unittest
from unittest.mock import Mock, patch, PropertyMock
from src.gene import Gene

# Some mock functions
def no_start_no_stop(self, item):
    return False

def start_stop(self, item):
    return True

def start_no_stop(self, item):
    if item == 'start_codon':
        return True
    return False

def no_start_stop(self, item):
    if item == 'stop_codon':
        return True
    return False

###################

class TestGene(unittest.TestCase):

    def setUp(self):
        self.gene1 = Gene()
        self.mrna1 = Mock()
        self.cds1 = Mock()
        self.exon1 = Mock()
        self.gene1.get_mrna = Mock(return_value=self.mrna1)
        self.mrna1.get_cds = Mock(return_value=self.cds1)
        self.mrna1.get_exon = Mock(return_value=self.exon1)

        self.gene1.start = 1
        self.gene1.end = 100
        self.gene1.attributes = {'ID':'foo_gene'}

        self.mrna1.start = 1
        self.mrna1.end = 100
        self.mrna1.attributes = {'ID':'m.foo'}

        self.cds1.start = 1
        self.cds1.end = 100
        self.cds1.phase = 0

        self.exon1.start = 1
        self.exon1.end = 100

    def test_from_gff_feature_success(self):
        gff_gene = Mock()
        gff_gene.type = "gene"
        
        tran_gene = Gene.from_gff_feature(gff_gene)
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

    def test_gene_to_tbl_nostart_nostop(self):
        expected = \
        "<1\t>100\tgene\n"\
        "\t\t\tlocus_tag\tfoo_gene\n"\
        "<1\t>100\tCDS\n"\
        "\t\t\tprotein_id\tm.foo\n"\
        "\t\t\tproduct\thypothetical protein\n"

        self.mrna1.__contains__ = no_start_no_stop

        tbl = self.gene1.to_tbl()
        self.assertEquals(tbl, expected)

    def test_gene_to_tbl_start_nostop(self):
        expected = \
        "1\t>100\tgene\n"\
        "\t\t\tlocus_tag\tfoo_gene\n"\
        "1\t>100\tCDS\n"\
        "\t\t\tprotein_id\tm.foo\n"\
        "\t\t\tproduct\thypothetical protein\n"

        self.mrna1.__contains__ = start_no_stop

        tbl = self.gene1.to_tbl()
        self.assertEquals(tbl, expected)

    def test_gene_to_tbl_nostart_stop(self):
        expected = \
        "<1\t100\tgene\n"\
        "\t\t\tlocus_tag\tfoo_gene\n"\
        "<1\t100\tCDS\n"\
        "\t\t\tprotein_id\tm.foo\n"\
        "\t\t\tproduct\thypothetical protein\n"

        self.mrna1.__contains__ = no_start_stop

        tbl = self.gene1.to_tbl()
        self.assertEquals(tbl, expected)

    def test_gene_to_tbl_start_stop(self):
        expected = \
        "1\t100\tgene\n"\
        "\t\t\tlocus_tag\tfoo_gene\n"\
        "1\t100\tCDS\n"\
        "\t\t\tprotein_id\tm.foo\n"\
        "\t\t\tproduct\thypothetical protein\n"

        self.mrna1.__contains__ = start_stop

        tbl = self.gene1.to_tbl()
        self.assertEquals(tbl, expected)

    def test_gene_to_tbl_genename(self):
        expected = \
        "<1\t>100\tgene\n"\
        "\t\t\tgene\tf00x4\n"\
        "\t\t\tlocus_tag\tfoo_gene\n"\
        "<1\t>100\tCDS\n"\
        "\t\t\tprotein_id\tm.foo\n"\
        "\t\t\tproduct\thypothetical protein\n"

        self.gene1.attributes["Name"] = "f00x4"

        self.mrna1.__contains__ = no_start_no_stop

        tbl = self.gene1.to_tbl()
        self.assertEquals(tbl, expected)

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

        self.gene1.get_mrna().attributes["Dbxref"] = "Pfam:foo,Pfam:dog,Pfam:baz"

        self.mrna1.__contains__ = no_start_no_stop

        tbl = self.gene1.to_tbl()
        self.assertEquals(tbl, expected)

    ### FIX PHASE TESTS ###

    def test_fix_phase(self):
        self.gene1.start = 2
        self.mrna1.start = 2
        self.mrna1.__contains__ = no_start_stop
        self.cds1.start = 2
        self.cds1.phase = 0
        self.assertEqual(self.cds1.phase, 0)
        self.gene1.fix_phase("ATGC")
        self.assertEqual(self.cds1.phase, 1)
        
    def test_fix_phase_to_two(self):
        self.gene1.start = 3
        self.mrna1.start = 3
        self.mrna1.__contains__ = no_start_stop
        self.cds1.start = 3
        self.cds1.phase = 0
        self.assertEqual(self.cds1.phase, 0)
        self.gene1.fix_phase("ATGC")
        self.assertEqual(self.cds1.phase, 2)
        
    def test_fix_phase_does_nothing_when_indices_too_large(self):
        self.gene1.start = 4
        self.mrna1.start = 4
        self.mrna1.__contains__ = start_stop
        self.cds1.start = 4
        self.cds1.phase = 0
        self.assertEqual(self.cds1.phase, 0)
        self.gene1.fix_phase("ATGC")
        self.assertEqual(self.cds1.phase, 0)

    def test_fix_phase_works_on_cds_only(self):
        self.gene1.start = 1
        self.mrna1.start = 1
        self.mrna1.__contains__ = no_start_stop
        self.cds1.start = 3
        self.cds1.phase = 0
        self.assertEqual(self.cds1.phase, 0)
        self.assertEqual(self.cds1.start, 3)
        self.gene1.fix_phase("ATGC")
        self.assertEqual(self.cds1.phase, 2)
        self.assertEqual(self.cds1.start, 1)

    def test_fix_phase_does_nothing_when_not_partial(self):
        self.gene1.start = 2
        self.mrna1.start = 2
        self.mrna1.__contains__ = start_stop
        self.cds1.start = 2
        self.cds1.phase = 0
        self.assertEqual(self.cds1.phase, 0)
        self.gene1.fix_phase("ATGC")
        self.assertEqual(self.cds1.phase, 0)

    def test_fix_phase_adjusts_end_on_3prime_partial(self):
        self.gene1.start = 2
        self.mrna1.start = 2
        self.mrna1.__contains__ = start_no_stop
        self.cds1.start = 2
        self.gene1.end = 2
        self.mrna1.end = 2
        self.cds1.end = 2
        self.assertEqual(self.cds1.end, 2)
        self.gene1.fix_phase("ATGC")
        self.assertEqual(self.cds1.end, 4)

    #### MAKE POSITIVE TESTS ####

    def test_make_positive(self):
        seq_len = 8

        gene = Gene(start=1, end=7, strand='-')

        mrna = Mock()
        mrna.type = 'mrna'
        mrna.make_positive = Mock()
        gene.add_child(mrna)

        gene.make_positive(seq_len)

        self.assertEqual(gene.start, 2)
        self.assertEqual(gene.end, 8)
        self.assertEqual(gene.strand, '+')
        mrna.make_positive.assertCalledWith(seq_len)

"""
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
