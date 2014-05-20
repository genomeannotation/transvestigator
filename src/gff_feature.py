#!/usr/bin/env python

class GFFFeature:

    def __init__(self):
        self.seqid = None
        self.source = None
        self.type = None
        self.start = None
        self.end = None
        self.score = None
        self.strand = None
        self.phase = None
        self.attributes = {} # Empty dictionary
