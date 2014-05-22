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


###################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestTranscriptFixer))
    return suite

if __name__ == '__main__':
    unittest.main()
