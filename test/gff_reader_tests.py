#!/usr/bin/env python

import unittest
import io
from src.gff_reader import *

class TestGFFReader(unittest.TestCase):

    def setUp(self):
        pass

    def test_parse_gff_attributes(self):
        attr = '\t; foo=dog;baz=bub;  \t\n'
        self.assertEquals(parse_gff_attributes(attr), {'foo':'dog', 'baz':'bub'})

    def test_read_gff(self):
        return
        
        test_gff = io.StringIO(\
        'seq0\tfoo_feature\t.\t.\t.\t.\t.\tID=foo_feat;Parent=foo_cds\n'+\
        'seq0\tGeibLabs\tgene\t1\t42\t.\t+\t.\tID=foo_gene\n'+\
        'seq0\tGeibLabs\tmRNA\t1\t42\t.\t+\t.\tID=foo_mrna;Parent=foo_gene\n'+\
        'seq0\tGeibLabs\texon\t1\t42\t0.9\t+\t.\tID=foo_exon0;Parent=foo_mrna\n'+\
        'seq0\tGeibLabs\texon\t1\t42\t0.9\t+\t.\tID=foo_exon1;Parent=foo_mrna\n'+\
        'seq0\tGeibLabs\tCDS\t1\t42\t.\t+\t0\tID=foo_cds;Parent=foo_mrna\n'\
        )

        gff = read_gff(test_gff)

        self.assertEquals(len(gff.gene), 1)

##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestGFFReader))
    return suite

if __name__ == '__main__':
    unittest.main()
