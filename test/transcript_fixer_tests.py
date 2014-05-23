import unittest
from mock import Mock, patch, PropertyMock
from src.transcript_fixer import fix_transcript 

class TestTranscriptFixer(unittest.TestCase):

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

        self.assertEquals(transcript.genes, [gene1])


###################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestTranscriptFixer))
    return suite

if __name__ == '__main__':
    unittest.main()
