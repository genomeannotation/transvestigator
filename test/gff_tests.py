#!/usr/bin/env python
# coding=utf-8

import unittest
import io
from src.gff import *


class TestGFF(unittest.TestCase):
    def setUp(self):
        pass

    def test_parse_gff_attributes(self):
        attr = "\t; foo=dog;baz=bub;  \t\n"
        self.assertEquals(parse_gff_attributes(attr), {"foo": "dog", "baz": "bub"})

    def test_parse_gff_attributes_multiple_identical_attributes(self):
        attr = "ID=foo_feature;foo=dog;foo=baz;foo=buzz"
        self.assertEquals(parse_gff_attributes(attr), {"ID": "foo_feature", "foo": "dog,baz,buzz"})

    def test_read_gff(self):
        test_gff = io.StringIO(
            "seq0\tGeibLabs\tfoo_feature\t.\t.\t.\t.\t.\tID=foo_feat;Parent=foo_cds\n"
            "seq0\tGeibLabs\tgene\t1\t42\t.\t+\t.\tID=foo_gene\n"
            "seq0\tGeibLabs\tmRNA\t1\t42\t.\t+\t.\tID=foo_mrna;Parent=foo_gene\n"
            "seq0\tGeibLabs\texon\t1\t42\t0.9\t+\t.\tID=foo_exon0;Parent=foo_mrna\n"
            "seq0\tGeibLabs\texon\t1\t24\t0.9\t+\t.\tID=foo_exon1;Parent=foo_mrna\n"
            "seq0\tGeibLabs\tCDS\t1\t42\t.\t+\t2\tID=foo_cds;Parent=foo_mrna\n")

        gff = read_gff(test_gff)

        self.assertEquals(len(gff["gene"]), 1)
        self.assertEquals(gff["gene"][0].source, "GeibLabs")
        self.assertEquals(gff["gene"][0].start, 1)
        self.assertEquals(gff["gene"][0].end, 42)
        self.assertEquals(gff["gene"][0].score, None)
        self.assertEquals(gff["gene"][0].strand, "+")
        self.assertEquals(gff["gene"][0].phase, 0)
        self.assertEquals(gff["gene"][0].attributes["ID"], "foo_gene")

        self.assertEquals(len(gff["gene"][0]["mrna"]), 1)
        self.assertEquals(len(gff["gene"][0]["mrna"][0]["exon"]), 2)
        self.assertEquals(len(gff["gene"][0]["mrna"][0]["cds"]), 1)

        self.assertEquals(gff["gene"][0]["mrna"][0]["exon"][0].start, 1)
        self.assertEquals(gff["gene"][0]["mrna"][0]["exon"][0].end, 42)
        self.assertEquals(gff["gene"][0]["mrna"][0]["exon"][0].score, 0.9)
        self.assertEquals(gff["gene"][0]["mrna"][0]["exon"][0].strand, "+")
        self.assertEquals(gff["gene"][0]["mrna"][0]["exon"][0].phase, 0)

        self.assertEquals(gff["gene"][0]["mrna"][0]["exon"][1].start, 1)
        self.assertEquals(gff["gene"][0]["mrna"][0]["exon"][1].end, 24)

        self.assertEquals(gff["gene"][0]["mrna"][0]["cds"][0].start, 1)
        self.assertEquals(gff["gene"][0]["mrna"][0]["cds"][0].end, 42)
        self.assertEquals(gff["gene"][0]["mrna"][0]["cds"][0].phase, 2)

        self.assertEquals(len(gff["gene"][0]["mrna"][0]["cds"][0]["foo_feature"]), 1)

    def test_read_gff_throws_on_missing_id(self):
        test_gff = io.StringIO("seq0\tGeibLabs\tfoo_feature\t.\t.\t.\t.\t.\t.\n")

        thrown = False
        try:
            read_gff(test_gff)
        except GFFError as error:
            self.assertEquals(str(error), "at line 0: feature has no ID attribute")
            thrown = True
        self.assertTrue(thrown)


##########################
# noinspection PyShadowingNames
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestGFF))
    return suite


if __name__ == "__main__":
    unittest.main()
