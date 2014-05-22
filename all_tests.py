#!/usr/bin/env python

# import all the lovely files
import unittest
import test.translator_tests
import test.fasta_tests
import test.transcript_checker_tests
import test.gff_reader_tests
import test.gff_feature_tests
import test.annotation_extractor_tests
import test.tbl_writer_tests
import test.transcript_builder_tests

# get suites from test modules
suites = [
test.translator_tests.suite(),\
test.fasta_tests.suite(),\
test.transcript_checker_tests.suite(),\
test.gff_reader_tests.suite(),\
test.gff_feature_tests.suite(),\
test.annotation_extractor_tests.suite(),\
test.tbl_writer_tests.suite(),\
test.transcript_builder_tests.suite()\
]

# collect suites in a TestSuite object
suite = unittest.TestSuite()
for s in suites:
    suite.addTest(s)

# run suite
unittest.TextTestRunner(verbosity=2).run(suite)
