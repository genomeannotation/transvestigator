#!/usr/bin/env python

# import all the lovely files
import unittest
import test.sequence_tests
import test.translator_tests
import test.fasta_reader_tests
import test.transcript_checker_tests

# get suites from test modules
suites = [
test.sequence_tests.suite(),\
test.translator_tests.suite(),\
test.fasta_reader_tests.suite(),\
test.transcript_checker_tests.suite()\
]

# collect suites in a TestSuite object
suite = unittest.TestSuite()
for s in suites:
    suite.addTest(s)

# run suite
unittest.TextTestRunner(verbosity=2).run(suite)
