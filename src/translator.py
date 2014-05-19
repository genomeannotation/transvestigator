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

