import unittest
from mock import Mock, patch, PropertyMock
from src.transcript_builder import build_transcript_dictionary 

class TestTranscriptBuilder(unittest.TestCase):

    def test_add_child(self):
        self.assertTrue(True)


###################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestTranscriptBuilder))
    return suite

if __name__ == '__main__':
    unittest.main()
