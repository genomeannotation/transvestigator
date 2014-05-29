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

