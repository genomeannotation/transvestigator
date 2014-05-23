#!/usr/bin/env python

class Sequence:

    def __init__(self, header="", bases=""):
        self.header = header
        self.bases = bases


class Transcript:

    def __init__(self, genes=None, sequence=None):
        if not genes:
            self.genes = []
        else:
            self.genes = genes
        if not sequence:
            self.sequence = Sequence()
        else:
            self.sequence = sequence


