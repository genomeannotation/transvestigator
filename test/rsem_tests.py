import io
import unittest
from mock import Mock, patch, PropertyMock
from src.rsem import read_rsem

class TestRsem(unittest.TestCase):

    def test_read_rsem(self):
        rsem = io.StringIO(\
        "transcript_id\tgene_id\tlength\teffective_length\texpected_count\tTPM\tFPKM\tIsoPct\n"\
        "c10000_g1_i1\tc10000_g1\t3938\t3717.75\t8342.00\t19.55\t14.87\t100.00\n"\
        "c10001_g1_i1\tc10001_g1\t774\t553.75\t20.24\t0.32\t0.24\t100.00\n"\
        )
        expected = [\
        ("c10000_g1_i1", "c10000_g1", 3938, 3717.75, 8342.00, 19.55, 14.87, 100.00),\
        ("c10001_g1_i1", "c10001_g1", 774, 553.75, 20.24, 0.32, 0.24, 100.00),\
        ]
        self.assertEquals(read_rsem(rsem), expected)


###################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestRsem))
    return suite

if __name__ == '__main__':
    unittest.main()
