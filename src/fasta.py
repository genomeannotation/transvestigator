#!/usr/bin/env python

from src.types import Sequence

def read_fasta(io_buffer):
    header = ''
    bases = ''
    seqs = []
    for line in io_buffer:
        if line[0] == '>':
            if len(header) > 0:
                # Save the data
                seqs.append(Sequence(header, bases))
            header = line[1:].strip().split()[0] # Get the next header
            bases = ''
        else:
            bases += line.strip()
    # Add the last sequence
    seqs.append(Sequence(header, bases))
    return seqs

def sequence_to_fasta(seq):
    result = '>' + seq.header + '\n'
    result += seq.bases + '\n'
    return result

def get_subsequence(seq, begin, end):
    return seq.bases[begin-1:end]