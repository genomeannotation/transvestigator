#!/usr/bin/env python

# import all the lovely files
import unittest
import test.sequtil_tests
import test.transcript_checker_tests
import test.gff_tests
import test.gff_feature_tests
import test.annotation_extractor_tests
import test.tbl_writer_tests
import test.transcript_builder_tests
import test.transcript_fixer_tests

# get suites from test modules
suites = [
test.sequtil_tests.suite(),\
test.transcript_checker_tests.suite(),\
test.gff_tests.suite(),\
test.gff_feature_tests.suite(),\
test.annotation_extractor_tests.suite(),\
test.tbl_writer_tests.suite(),\
test.transcript_builder_tests.suite(),\
test.transcript_fixer_tests.suite(),\
]

# collect suites in a TestSuite object
suite = unittest.TestSuite()
for s in suites:
    suite.addTest(s)

# run suite
unittest.TextTestRunner(verbosity=2).run(suite)
