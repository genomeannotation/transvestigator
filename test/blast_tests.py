import io
import unittest
from unittest.mock import Mock, patch, PropertyMock
from src.blast import read_blast

class TestBlast(unittest.TestCase):

    def test_read_blast(self):
        blast_input = io.StringIO(
                "c25_g1_i1\tsp|P28755|SODC_CERCA\t81.25\t96\t18\t0\t760\t473\t1\t96\t6e-68\t167\n"
                "c25_g1_i2\tsp|P28755|SODC_CERCA\t84.97\t153\t23\t0\t699\t241\t1\t153\t6e-91\t273\n"
                "c42_g1_i1\tsp|Q9QXZ7|NR2E3_MOUSE\t60.27\t73\t29\t0\t2\t220\t49\t121\t2e-22\t95.9\n"
        )
        expected = [
            ("c25_g1_i1", "sp|P28755|SODC_CERCA", 81.25, 96, 18, 0, 760, 473, 1, 96, 6e-68, 167),
            ("c25_g1_i2", "sp|P28755|SODC_CERCA", 84.97, 153, 23, 0, 699, 241, 1, 153, 6e-91, 273),
            ("c42_g1_i1", "sp|Q9QXZ7|NR2E3_MOUSE", 60.27, 73, 29, 0, 2, 220, 49, 121, 2e-22, 95.9),
        ]
        actual = read_blast(blast_input)
        self.assertEquals(actual, expected)


###################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestBlast))
    return suite

if __name__ == '__main__':
    unittest.main()
