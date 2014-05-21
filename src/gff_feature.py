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
        child_type = child.type.lower() # For API, use lowercase names
        if hasattr(self, child_type): # Already have children of this type, append it to the list
            getattr(self, child_type).append(child)
        else: # No children with this type yet, make new list
            setattr(self, child_type, [child])

    def has_child(self, child_type):
        if hasattr(self, child_type):
            return True
        return False
