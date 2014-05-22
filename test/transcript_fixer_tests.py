import unittest
from mock import Mock, patch, PropertyMock
from src.transcript_fixer import fix_transcript 

class TestTranscriptFixer(unittest.TestCase):

    def test_fix_transcript(self):
        self.assertTrue(True)


###################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestTranscriptFixer))
    return suite

if __name__ == '__main__':
    unittest.main()
