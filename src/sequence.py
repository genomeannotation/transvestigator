#!/usr/bin/env python


class Sequence:

    def __init__(self, header="", bases=""):
        self.header = header
        self.bases = bases

    def to_fasta(self):
        result = '>' + self.header + '\n'
        result += self.bases + '\n'
        return result

