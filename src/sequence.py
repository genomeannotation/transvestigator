#!/usr/bin/env python


class Sequence:

    def __init__(self, header="", bases=""):
        self.header = header
        self.bases = bases

    def to_fasta(self):
        result = '>' + self.header + '\n'
        result += self.bases + '\n'
        return result

####################################################

def read_fasta(io_buffer):
    header = ''
    bases = ''
    seqs = {}
    for line in io_buffer:
        if line[0] == '>':
            if len(header) > 0:
                # Save the data
                seqs[header] = Sequence(header, bases)
            header = line[1:].strip().split()[0] # Get the next header
            bases = ''
        else:
            bases += line.strip()
    # Add the last sequence
    seqs[header] = Sequence(header, bases)
    return seqs

