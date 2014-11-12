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

        self.mrna0 = Mock()
        self.mrna1 = Mock()

        self.cds0 = Mock()
        self.cds0.length = Mock(return_value=5)
        self.cds1 = Mock()
        self.cds1.length = Mock(return_value=6)

        self.transcript.genes = [self.gene0, self.gene1]
        self.gene0.__getitem__ = Mock(return_value=[self.mrna0])
        self.gene0.__contains__ = Mock()
        self.gene0.__contains__.assertCalledWith("mrna")
        self.gene1.__getitem__ = Mock(return_value=[self.mrna1])
        self.gene1.__iter__ = Mock(return_value=[self.mrna1].__iter__())
        #self.gene1["mrna"] = [self.mrna1]
        self.mrna0.__getitem__ = Mock(return_value=[self.cds0])
        self.mrna1.__getitem__ = Mock(return_value=[self.cds1])
        #self.mrna0["cds"] = [self.cds0]
        #self.mrna1["cds"] = [self.cds1]

    def test_fix_multiple_genes(self):
        self.transcript.fix_multiple_genes()
        self.assertEquals(self.transcript.genes, [self.gene1])

    def test_fix_phase(self):
        self.gene0.start = 2
        self.mrna0.start = 2
        del self.mrna0.start_codon
        self.cds0.start = 2
        self.cds0.phase = 0
        self.assertEqual(self.cds0.phase, 0)
        self.transcript.fix_phase()
        self.assertEqual(self.cds0.phase, 1)
        
    def test_fix_phase_to_two(self):
        self.gene0.start = 3
        self.mrna0.start = 3
        del self.mrna0.start_codon
        self.cds0.start = 3
        self.cds0.phase = 0
        self.assertEqual(self.cds0.phase, 0)
        self.transcript.fix_phase()
        self.assertEqual(self.cds0.phase, 2)
        
    def test_fix_phase_does_nothing_when_indices_too_large(self):
        self.gene0.start = 4
        self.mrna0.start = 4
        self.cds0.start = 4
        self.cds0.phase = 0
        self.assertEqual(self.cds0.phase, 0)
        self.transcript.fix_phase()
        self.assertEqual(self.cds0.phase, 0)

    def test_fix_phase_works_on_cds_only(self):
        self.gene0.start = 1
        self.mrna0.start = 1
        self.cds0.start = 3
        self.cds0.phase = 0
        del self.mrna0.start_codon
        self.assertEqual(self.cds0.phase, 0)
        self.assertEqual(self.cds0.start, 3)
        self.transcript.fix_phase()
        self.assertEqual(self.cds0.phase, 2)
        self.assertEqual(self.cds0.start, 1)

    def test_fix_phase_does_nothing_when_not_partial(self):
        self.gene0.start = 2
        self.mrna0.start = 2
        self.cds0.start = 2
        self.cds0.phase = 0
        self.assertEqual(self.cds0.phase, 0)
        self.transcript.fix_phase()
        self.assertEqual(self.cds0.phase, 0)

    def test_fix_phase_adjusts_end_on_3prime_partial(self):
        self.gene0.start = 2
        self.mrna0.start = 2
        self.cds0.start = 2
        self.gene0.end = 2
        self.mrna0.end = 2
        self.cds0.end = 2
        del self.mrna0.stop_codon
        self.assertEqual(self.cds0.end, 2)
        self.transcript.fix_phase()
        self.assertEqual(self.cds0.end, 6)


    ### TO TBL TESTS ###

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
        gff_cds0.phase = 0
        
        gff_gene0["mrna"] = [gff_mrna0]
        gff_mrna0["cds"] = [gff_cds0]

        del gff_mrna0.start_codon
        del gff_mrna0.stop_codon

        return gff_gene0

    def test_gene_to_tbl_nostart_nostop(self):
        expected = \
        "<1\t>100\tgene\n"\
        "\t\t\tlocus_tag\tfoo_gene\n"\
        "<1\t>100\tCDS\n"\
        "\t\t\tprotein_id\tm.foo\n"\
        "\t\t\tproduct\thypothetical protein\n"

        gff_gene0 = self.create_fake_gene()

        tbl = gene_to_tbl(gff_gene0)
        self.assertEquals(tbl, expected)

    def test_gene_to_tbl_start_nostop(self):
        expected = \
        "1\t>100\tgene\n"\
        "\t\t\tlocus_tag\tfoo_gene\n"\
        "1\t>100\tCDS\n"\
        "\t\t\tprotein_id\tm.foo\n"\
        "\t\t\tproduct\thypothetical protein\n"

        gff_gene0 = self.create_fake_gene()

        gff_gene0["mrna"][0].start_codon = [Mock()]

        tbl = gene_to_tbl(gff_gene0)
        self.assertEquals(tbl, expected)

    def test_gene_to_tbl_nostart_stop(self):
        expected = \
        "<1\t100\tgene\n"\
        "\t\t\tlocus_tag\tfoo_gene\n"\
        "<1\t100\tCDS\n"\
        "\t\t\tprotein_id\tm.foo\n"\
        "\t\t\tproduct\thypothetical protein\n"

        gff_gene0 = self.create_fake_gene()

        gff_gene0["mrna"][0].stop_codon = [Mock()]

        tbl = gene_to_tbl(gff_gene0)
        self.assertEquals(tbl, expected)

    def test_gene_to_tbl_start_stop(self):
        expected = \
        "1\t100\tgene\n"\
        "\t\t\tlocus_tag\tfoo_gene\n"\
        "1\t100\tCDS\n"\
        "\t\t\tprotein_id\tm.foo\n"\
        "\t\t\tproduct\thypothetical protein\n"

        gff_gene0 = self.create_fake_gene()

        gff_gene0["mrna"][0].start_codon = [Mock()]
        gff_gene0["mrna"][0].stop_codon = [Mock()]

        tbl = gene_to_tbl(gff_gene0)
        self.assertEquals(tbl, expected)

    def test_gene_to_tbl_genename(self):
        expected = \
        "<1\t>100\tgene\n"\
        "\t\t\tgene\tf00x4\n"\
        "\t\t\tlocus_tag\tfoo_gene\n"\
        "<1\t>100\tCDS\n"\
        "\t\t\tprotein_id\tm.foo\n"\
        "\t\t\tproduct\thypothetical protein\n"

        gff_gene0 = self.create_fake_gene()

        gff_gene0.attributes["Name"] = "f00x4"

        tbl = gene_to_tbl(gff_gene0)
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

        gff_gene0 = self.create_fake_gene()

        gff_gene0["mrna"][0].attributes["Dbxref"] = "Pfam:foo,Pfam:dog,Pfam:baz"

        tbl = gene_to_tbl(gff_gene0)
        self.assertEquals(tbl, expected)

    def test_gene_to_tbl_codon_start(self):
        expected = \
        "<1\t>100\tgene\n"\
        "\t\t\tlocus_tag\tfoo_gene\n"\
        "<1\t>100\tCDS\n"\
        "\t\t\tcodon_start\t2\n"\
        "\t\t\tprotein_id\tm.foo\n"\
        "\t\t\tproduct\thypothetical protein\n"\

        gff_gene0 = self.create_fake_gene()

        gff_gene0["mrna"][0]["cds"][0].phase = 1

        tbl = gene_to_tbl(gff_gene0)
        self.assertEquals(tbl, expected)

    def test_to_tbl(self):
        expected = \
        ">Feature foo_seq\n"\
        "1\t4\tREFERENCE\n"\
        "<1\t>100\tgene\n"\
        "\t\t\tlocus_tag\tfoo_gene\n"\
        "<1\t>100\tCDS\n"\
        "\t\t\tprotein_id\tm.foo\n"\
        "\t\t\tproduct\thypothetical protein\n"

        gff_gene0 = self.create_fake_gene()

        sequence = Mock()
        sequence.header = "foo_seq"
        sequence.bases = "ATGC"

        transcript = Transcript([gff_gene0], sequence) 

        tbl = transcript.to_tbl()
        self.assertEquals(tbl, expected)

    def test_remove_contig_from_gene_id(self):
        expected = 'g.123'
        gff_gene0 = self.create_fake_gene()
        gff_gene0.attributes['ID'] = 'contig123|g.123'

        sequence = Mock()
        sequence.header = "foo_seq"
        sequence.bases = "ATGC"

        transcript = Transcript([gff_gene0], sequence) 
        transcript.remove_contig_from_gene_id()

        self.assertEquals(transcript.genes[0].attributes['ID'], expected)

    #### MAKE POSITIVE TESTS ####

    def test_make_positive(self):
        seq = Mock()
        seq.header = "foo_seq"
        seq.bases = "ACACACTT"

        gene = Mock()
        gene.start = 1
        gene.end = 7
        gene.strand = '-'
        gene["mrna"] = [Mock()]
        gene["mrna"][0].start = 2
        gene["mrna"][0].end = 6
        gene["mrna"][0]["cds"] = [Mock()]
        gene["mrna"][0]["cds"][0].start = 2
        gene["mrna"][0]["cds"][0].end = 6
        gene["mrna"][0]["exon"] = [Mock()]
        gene["mrna"][0]["exon"][0].start = 2
        gene["mrna"][0]["exon"][0].end = 6
        
        tran = Transcript([gene], seq)
        tran.make_positive()

        self.assertEqual(tran.sequence.bases, "AAGTGTGT")
        self.assertEqual(tran.genes[0].start, 2)
        self.assertEqual(tran.genes[0].end, 8)
        self.assertEqual(tran.genes[0].strand, '+')
        self.assertEqual(tran.genes[0]["mrna"][0].start, 3)
        self.assertEqual(tran.genes[0]["mrna"][0].end, 7)
        self.assertEqual(tran.genes[0]["mrna"][0].strand, '+')
        self.assertEqual(tran.genes[0]["mrna"][0]["cds"][0].start, 3)
        self.assertEqual(tran.genes[0]["mrna"][0]["cds"][0].end, 7)
        self.assertEqual(tran.genes[0]["mrna"][0]["cds"][0].strand, '+')
        self.assertEqual(tran.genes[0]["mrna"][0]["exon"][0].start, 3)
        self.assertEqual(tran.genes[0]["mrna"][0]["exon"][0].end, 7)
        self.assertEqual(tran.genes[0]["mrna"][0]["exon"][0].strand, '+')

    #### FIX LENGTHS TESTS ####

    def test_fix_feature_lengths(self):
        seq = Mock()
        seq.header = "foo_seq"
        seq.bases = "ACACACTT"

        gene = Mock()
        gene.start = 1
        gene.end = 7
        gene["mrna"] = [Mock()]
        gene["mrna"][0].start = 1
        gene["mrna"][0].end = 9
        gene["mrna"][0]["cds"] = [Mock()]
        gene["mrna"][0]["cds"][0].start = 2
        gene["mrna"][0]["cds"][0].end = 10
        gene["mrna"][0]["exon"] = [Mock()]
        gene["mrna"][0]["exon"][0].start = 3
        gene["mrna"][0]["exon"][0].end = 11
        
        tran = Transcript([gene], seq)
        tran.fix_feature_lengths()

        self.assertEqual(tran.sequence.bases, "ACACACTT")
        self.assertEqual(tran.genes[0].start, 1)
        self.assertEqual(tran.genes[0].end, 7)
        self.assertEqual(tran.genes[0]["mrna"][0].start, 1)
        self.assertEqual(tran.genes[0]["mrna"][0].end, 6)
        self.assertEqual(tran.genes[0]["mrna"][0]["cds"][0].start, 2)
        self.assertEqual(tran.genes[0]["mrna"][0]["cds"][0].end, 7)
        self.assertEqual(tran.genes[0]["mrna"][0]["exon"][0].start, 3)
        self.assertEqual(tran.genes[0]["mrna"][0]["exon"][0].end, 8)

    #### MATCH CDS AND EXON END TESTS ####

    def test_match_cds_and_exon_end(self):
        seq = Sequence("foo_seq", "ATGNNN") 
        gene = Mock()
        gene["mrna"] = [Mock()]
        gene["mrna"][0]["cds"] = [Mock()]
        gene["mrna"][0]["cds"][0].start = 1
        gene["mrna"][0]["cds"][0].end = 5
        gene["mrna"][0]["exon"] = [Mock()]
        gene["mrna"][0]["exon"][0].start = 1
        gene["mrna"][0]["exon"][0].end = 6
        gene["mrna"][0].attributes = {"ID": "foo_mrna", "Parent": "foo_gene"}

        # method doesn't touch indices if a stop codon is present;
        # mock will return True to the 'hasattr' question unless I do this
        delattr(gene["mrna"][0], "stop_codon")
        
        tran = Transcript([gene], seq)
        self.assertEquals(5, gene["mrna"][0]["cds"][0].end)
        tran.match_cds_and_exon_end()
        self.assertEquals(6, gene["mrna"][0]["cds"][0].end)

    def test_match_cds_and_exon_end_does_nothing_if_stop_codon_present(self):
        seq = Sequence("foo_seq", "ATGNNN") 
        gene = Mock()
        gene["mrna"] = [Mock()]
        gene["mrna"][0]["cds"] = [Mock()]
        gene["mrna"][0]["cds"][0].start = 1
        gene["mrna"][0]["cds"][0].end = 4
        gene["mrna"][0]["exon"] = [Mock()]
        gene["mrna"][0]["exon"][0].start = 1
        gene["mrna"][0]["exon"][0].end = 6
        gene["mrna"][0].attributes = {"ID": "foo_mrna", "Parent": "foo_gene"}
        # Note: not necessary to mock stop_codon feature, 
        # since hasattr(Mock, "wateva") returns True
        tran = Transcript([gene], seq)
        self.assertEquals(4, gene["mrna"][0]["cds"][0].end)
        tran.match_cds_and_exon_end()
        self.assertEquals(4, gene["mrna"][0]["cds"][0].end) # cds end is unchanged


    #### STARTS AND STOPS TESTS ####
    
    def test_create_starts_and_stops_creates_a_start(self):
        seq = Sequence("foo_seq", "ATGNNN") 
        gene = Mock()
        gene["mrna"] = [Mock()]
        gene["mrna"][0]["cds"] = [Mock()]
        gene["mrna"][0]["cds"][0].start = 1
        gene["mrna"][0]["cds"][0].end = 6
        gene["mrna"][0].attributes = {"ID": "foo_mrna", "Parent": "foo_gene"}
        tran = Transcript([gene], seq)
        self.assertEquals(0, len(gene["mrna"][0].mock_calls))
        tran.create_starts_and_stops()
        self.assertEquals(1, len(gene["mrna"][0].mock_calls))
        self.assertEquals("add_child", gene["mrna"][0].mock_calls[0][0])

    def test_create_starts_and_stops_creates_a_start_reverse_complement(self):
        seq = Sequence("foo_seq", "NNNCAT") 
        gene = Mock()
        gene["mrna"] = [Mock()]
        gene["mrna"][0]["cds"] = [Mock()]
        gene["mrna"][0]["cds"][0].start = 1
        gene["mrna"][0]["cds"][0].end = 6
        gene["mrna"][0]["cds"][0].strand = "-"
        gene["mrna"][0].attributes = {"ID": "foo_mrna", "Parent": "foo_gene"}
        tran = Transcript([gene], seq)
        self.assertEquals(0, len(gene["mrna"][0].mock_calls))
        tran.create_starts_and_stops()
        self.assertEquals(1, len(gene["mrna"][0].mock_calls))
        self.assertEquals("add_child", gene["mrna"][0].mock_calls[0][0])

    def test_create_starts_and_stops_creates_a_stop(self):
        seq = Sequence("foo_seq", "NNNTAG") 
        gene = Mock()
        gene["mrna"] = [Mock()]
        gene["mrna"][0]["cds"] = [Mock()]
        gene["mrna"][0]["cds"][0].start = 1
        gene["mrna"][0]["cds"][0].end = 6
        gene["mrna"][0].attributes = {"ID": "foo_mrna", "Parent": "foo_gene"}
        tran = Transcript([gene], seq)
        self.assertEquals(0, len(gene["mrna"][0].mock_calls))
        tran.create_starts_and_stops()
        self.assertEquals(1, len(gene["mrna"][0].mock_calls))
        self.assertEquals("add_child", gene["mrna"][0].mock_calls[0][0])

    def test_create_starts_and_stops_creates_a_stop_reverse_complement(self):
        seq = Sequence("foo_seq", "CTANNN") 
        gene = Mock()
        gene["mrna"] = [Mock()]
        gene["mrna"][0]["cds"] = [Mock()]
        gene["mrna"][0]["cds"][0].start = 1
        gene["mrna"][0]["cds"][0].end = 6
        gene["mrna"][0]["cds"][0].strand = "-"
        gene["mrna"][0].attributes = {"ID": "foo_mrna", "Parent": "foo_gene"}
        tran = Transcript([gene], seq)
        self.assertEquals(0, len(gene["mrna"][0].mock_calls))
        tran.create_starts_and_stops()
        self.assertEquals(1, len(gene["mrna"][0].mock_calls))
        self.assertEquals("add_child", gene["mrna"][0].mock_calls[0][0])


##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestTranscript))
    return suite

if __name__ == '__main__':
    unittest.main()
