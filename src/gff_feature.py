#!/usr/bin/env python

class GFFFeature:

    def __init__(self, seqid=None, source=None, type=None, start=None, end=None, score=None, strand=None, phase=0, attributes=None):
        self.seqid = seqid
        self.source = source
        self.type = type
        self.start = start
        self.end = end
        self.score = score
        self.strand = strand
        self.phase = phase
        if not attributes:
            self.attributes = {} # Empty dictionary
        else:
            self.attributes = attributes

    def add_child(self, child):
        child_type = child.type.lower() # For API, use lowercase names
        if hasattr(self, child_type): # Already have children of this type, append it to the list
            getattr(self, child_type).append(child)
        else: # No children with this type yet, make new list
            setattr(self, child_type, [child])

    def add_annotation(self, key, value):
        if key in self.attributes:
            self.attributes[key] += ","+value
        else:
            self.attributes[key] = value

    def length(self):
        return self.end-self.start+1
