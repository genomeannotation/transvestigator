import unittest
from unittest.mock import Mock, patch, PropertyMock
from src.transcript import Transcript
from src.transcript_builder import build_transcript_dictionary 

class TestTranscriptBuilder(unittest.TestCase):

    def test_build_transcript_dictionary(self):
        seq0 = Mock()
        seq0.header = "seq0" 
        seq1 = Mock()
        seq1.header = "seq1"
        
        gene0 = Mock()
        gene0.type = "gene"
        gene0.seqid = "seq0"
        gene0.children = {'mrna':[]}
        gene1 = Mock()
        gene1.type = "gene"
        gene1.seqid = "seq0"
        gene1.children = {'mrna':[]}
        gene2 = Mock()
        gene2.type = "gene"
        gene2.seqid = "seq1"
        gene2.children = {'mrna':[]}

        transcripts = build_transcript_dictionary({"seq0":seq0, "seq1":seq1}, [gene0, gene1, gene2])

        self.assertEquals(2, len(transcripts["seq0"].genes))
        self.assertEquals(1, len(transcripts["seq1"].genes))


###################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestTranscriptBuilder))
    return suite

if __name__ == '__main__':
    unittest.main()
