import unittest
from mock import Mock, patch, PropertyMock
from src.transcript_fixer import fix_transcript, fix_phase

class TestTranscriptFixer(unittest.TestCase):

    def setUp(self):
        self.transcript = Mock()
        self.gene0 = Mock()
        self.gene1 = Mock()

        self.mrna0 = Mock()
        self.mrna1 = Mock()

        self.cds0 = Mock()
        self.cds0.length = Mock(return_value=5)
        self.cds1 = Mock()
        self.cds1.length = Mock(return_value=6)

        self.transcript.genes = [self.gene0, self.gene1]
        self.gene0.mrna = [self.mrna0]
        self.gene1.mrna = [self.mrna1]
        self.mrna0.cds = [self.cds0]
        self.mrna1.cds = [self.cds1]

    def test_fix_transcript(self):
        fix_transcript(self.transcript)
        self.assertEquals(self.transcript.genes, [self.gene1])

    def test_fix_phase(self):
        self.gene0.start = 2
        self.mrna0.start = 2
        self.cds0.start = 2
        self.cds0.phase = 0
        self.assertEqual(self.cds0.phase, 0)
        fix_phase(self.transcript)
        self.assertEqual(self.cds0.phase, 1)
        
    def test_fix_phase_to_two(self):
        self.gene0.start = 3
        self.mrna0.start = 3
        self.cds0.start = 3
        self.cds0.phase = 0
        self.assertEqual(self.cds0.phase, 0)
        fix_phase(self.transcript)
        self.assertEqual(self.cds0.phase, 2)
        
    def test_fix_phase_does_nothing_when_indices_too_large(self):
        self.gene0.start = 4
        self.mrna0.start = 4
        self.cds0.start = 4
        self.cds0.phase = 0
        self.assertEqual(self.cds0.phase, 0)
        fix_phase(self.transcript)
        self.assertEqual(self.cds0.phase, 0)
        

###################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestTranscriptFixer))
    return suite

if __name__ == '__main__':
    unittest.main()
