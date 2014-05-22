import unittest
from mock import Mock, patch, PropertyMock
from src.transcript_fixer import fix_transcript 

class TestTranscriptFixer(unittest.TestCase):

    def test_fix_transcript_throws_on_no_mrna(self):
        transcript = Mock()
        transcript.sequence.header = "foo_seq"
        gene = Mock()
        gene.attributes = {"ID":"foo_gene"}
        del gene.mrna
        transcript.genes = [gene]

        thrown = False
        try:
            fix_transcript(transcript)
        except Exception as error:
            self.assertEqual(str(error), "can't fix transcript foo_seq because gene foo_gene has no mRNA")
            thrown = True
        self.assertTrue(thrown)
        
    def test_fix_transcript_throws_on_multiple_mrnas(self):
        transcript = Mock()
        transcript.sequence.header = "foo_seq"
        gene = Mock()
        gene.attributes = {"ID":"foo_gene"}
        gene.mrna = [Mock(), Mock()]
        transcript.genes = [gene]

        thrown = False
        try:
            fix_transcript(transcript)
        except Exception as error:
            self.assertEqual(str(error), "can't fix transcript foo_seq because gene foo_gene has multiple mRNAs")
            thrown = True
        self.assertTrue(thrown)

    def test_fix_transcript_throws_on_no_cds(self):
        transcript = Mock()
        transcript.sequence.header = "foo_seq"
        gene = Mock()
        mrna = Mock()
        mrna.attributes = {"ID":"foo_mrna"}
        del mrna.cds
        gene.mrna = [mrna]
        transcript.genes = [gene]

        thrown = False
        try:
            fix_transcript(transcript)
        except Exception as error:
            self.assertEqual(str(error), "can't fix transcript foo_seq because mRNA foo_mrna has no CDS")
            thrown = True
        self.assertTrue(thrown)
        
    def test_fix_transcript_throws_on_multiple_cds(self):
        transcript = Mock()
        transcript.sequence.header = "foo_seq"
        gene = Mock()
        gene.attributes = {"ID":"foo_gene"}
        mrna = Mock()
        mrna.attributes = {"ID":"foo_mrna"}
        mrna.cds = [Mock(), Mock()]
        gene.mrna = [mrna]
        transcript.genes = [gene]

        thrown = False
        try:
            fix_transcript(transcript)
        except Exception as error:
            self.assertEqual(str(error), "can't fix transcript foo_seq because mRNA foo_mrna has multiple CDSs")
            thrown = True
        self.assertTrue(thrown)

    def test_fix_transcript(self):
        transcript = Mock()
        gene0 = Mock()
        gene1 = Mock()

        mrna0 = Mock()
        mrna1 = Mock()

        cds0 = Mock()
        cds0.start = 1
        cds0.end = 5
        cds1 = Mock()
        cds1.start = 1
        cds1.end = 6

        transcript.genes = [gene0, gene1]
        gene0.mrna = [mrna0]
        gene1.mrna = [mrna1]
        mrna0.cds = [cds0]
        mrna1.cds = [cds1]

        fix_transcript(transcript)

        self.assertEquals(transcript.genes, [gene0])


###################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestTranscriptFixer))
    return suite

if __name__ == '__main__':
    unittest.main()
