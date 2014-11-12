#!/usr/bin/env python

# import all the lovely files
import unittest
import test.sequtil_tests
import test.sequence_tests
import test.gff_tests
import test.gff_feature_tests
import test.transcript_tests
import test.transcript_builder_tests
import test.annotation_tests
import test.rsem_tests
import test.blast_tests
import test.config_tests
import test.gene_tests
import test.mrna_tests

# get suites from test modules
suites = [
test.sequtil_tests.suite(),\
test.sequence_tests.suite(),\
test.gff_tests.suite(),\
test.gff_feature_tests.suite(),\
test.transcript_tests.suite(),\
test.transcript_builder_tests.suite(),\
test.annotation_tests.suite(),\
test.rsem_tests.suite(),\
test.blast_tests.suite(),\
test.config_tests.suite(),\
test.gene_tests.suite(),\
test.mrna_tests.suite(),\
]

# collect suites in a TestSuite object
suite = unittest.TestSuite()
for s in suites:
    suite.addTest(s)

# run suite
unittest.TextTestRunner(verbosity=2).run(suite)
