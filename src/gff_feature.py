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
        self.children = []

    def add_child(self, child):
        child_type = child.type.lower() # For API, use lowercase names
        if hasattr(self, child_type): # Already have children of this type, append it to the list
            getattr(self, child_type).append(child)
        else: # No children with this type yet, make new list
            setattr(self, child_type, [child])
        self.children.append(child)

    def add_annotation(self, key, value):
        if key in self.attributes:
            self.attributes[key] += ","+value
        else:
            self.attributes[key] = value

    def length(self):
        return self.end-self.start+1

    def write(self):
        out = ""
        if self.seqid != None:
            out += self.seqid+"\t"
        else:
            out += ".\t"
        if self.source != None:
            out += self.source+"\t"
        else:
            out += ".\t"
        if self.type != None:
            out += self.type+"\t"
        else:
            out += ".\t"
        if self.start != None:
            out += str(self.start)+"\t"
        else:
            out += ".\t"
        if self.end != None:
            out += str(self.end)+"\t"
        else:
            out += ".\t"
        if self.score != None:
            out += str(self.score)+"\t"
        else:
            out += ".\t"
        if self.strand != None:
            out += self.strand+"\t"
        else:
            out += ".\t"
        if self.phase != None:
            out += str(self.phase)+"\t"
        else:
            out += ".\t"
        out += ";".join([key+"="+val for key, val in sorted(self.attributes.items())])
        return out
