#!/usr/bin/env python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4

def reverse_complement(seq):
    bases = ['a', 'c', 'g', 't', 'n', 'A', 'C', 'G', 'T', 'N']
    complements = ['t', 'g', 'c', 'a', 'n', 'T', 'G', 'C', 'A', 'N']
    rev_comp_dict = dict(zip(bases, complements))
    # Convert mixed or illegal bases to 'N'
    for i, base in enumerate(seq):
        if base not in 'actgnACTGN':
            seq = seq[0:i] + 'N' + seq[i+1:]
    return ''.join([rev_comp_dict.get(base) for base in reversed(seq)])

def has_start_codon(seq):
    return seq[:3].lower() == 'aug' or seq[:3].lower() == 'atg'

def has_stop_codon(seq):
    last3 = seq[-3:].lower()
    if last3 == 'tag':
        return True
    elif last3 == 'taa':
        return True
    elif last3 == 'tga':
        return True
    else:
        return False

def get_subsequence(bases, begin, end):
    return bases[begin-1:end]

def overlap(indices1, indices2):
    # Case 1:
    #   indices1    --------
    #   indices2  ------
    if indices1[0] > indices2[0] and indices1[0] <= indices2[1]:
        return True
    # Case 2:
    #   indices1  --------
    #   indices2      ------
    elif indices1[1] >= indices2[0] and indices1[1] < indices2[1]:
        return True
    else:
        return False

def nested(indices1, indices2):
    if indices1[0] >= indices2[0] and indices1[1] <= indices2[1]:
        return True
    elif indices2[0] >= indices1[0] and indices2[1] <= indices1[1]:
        return True
    else:
        return False
