import unittest
from mock import Mock, patch, PropertyMock
from src.types import Transcript
from src.transcript_builder import build_transcript_dictionary 

class TestTranscriptBuilder(unittest.TestCase):

    def test_build_transcript_dictionary(self):
        seq0 = Mock()
        seq0.header = "seq0" 
        seq1 = Mock()
        seq1.header = "seq1"
        
        gene0 = Mock()
        gene0.seqid = "seq0"
        gene1 = Mock()
        gene1.seqid = "seq0"
        gene2 = Mock()
        gene2.seqid = "seq1"

        expected = {"seq0":Transcript([gene0, gene1], seq0), "seq1":Transcript([gene2], seq1)}

        transcripts = build_transcript_dictionary({"seq0":seq0, "seq1":seq1}, [gene0, gene1, gene2])

        self.assertEquals(expected["seq0"].genes, transcripts["seq0"].genes)
        self.assertEquals(expected["seq1"].sequence, transcripts["seq1"].sequence)


###################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestTranscriptBuilder))
    return suite

if __name__ == '__main__':
    unittest.main()
