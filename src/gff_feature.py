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

    def add_child(self, child):
        type = child.type.lower() # For API, use lowercase names
        if hasattr(self, type): # Already have children of this type, append it to the list
            getattr(self, type).append(child)
        else: # No children with this type yet, make new list
            setattr(self, type, [child])
