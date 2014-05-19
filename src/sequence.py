#!/usr/bin/env python

class Sequence:

    def __init__(self, header="", bases=""):
        self.header = header
        self.bases = bases

    def to_fasta(self):
        result = '>' + self.header + '\n'
        result += self.bases + '\n'
        return result

    def get_subseq(self, start=1, stop=None):
        """Returns subseq from 1-based index "start" to "stop", inclusive."""
        if not stop:
            stop = len(self.bases)
        if stop > len(self.bases):
            return ""
        return self.bases[start-1:stop]

    def to_boulder_io(self):
        result = "SEQUENCE_ID="
        result += self.header + "\n"
        result += "SEQUENCE_TEMPLATE="
        result += self.bases + "\n=\n"
        return result
    
