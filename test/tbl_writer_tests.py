#!/usr/bin/env python

import unittest

class TestTblWriter(unittest.TestCase):

    def setUp(self):
        pass
        

##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestTblWriter))
    return suite

if __name__ == '__main__':
    unittest.main()
