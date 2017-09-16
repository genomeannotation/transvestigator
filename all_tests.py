#!/usr/bin/env python
# coding=utf-8

# import all the lovely files
import unittest
from test import sequtil_tests, sequence_tests, gff_tests, gff_feature_tests, transcript_tests, \
    transcript_builder_tests, annotation_tests, rsem_tests, blast_tests, config_tests, gene_tests, mrna_tests

# get suites from test modules

suite = unittest.TestSuite()

suite.addTest(sequtil_tests.suite())
suite.addTest(sequence_tests.suite())
suite.addTest(gff_tests.suite())
suite.addTest(gff_feature_tests.suite())
suite.addTest(transcript_tests.suite())
suite.addTest(transcript_builder_tests.suite())
suite.addTest(annotation_tests.suite())
suite.addTest(rsem_tests.suite())
suite.addTest(blast_tests.suite())
suite.addTest(config_tests.suite())
suite.addTest(gene_tests.suite())
suite.addTest(mrna_tests.suite())

# run suite
unittest.TextTestRunner(verbosity=2).run(suite)
