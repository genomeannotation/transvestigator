#!/usr/bin/env python

def fix_transcript(transcript):
    longest = None
    length = 0
    for gene in transcript.genes:
        this_length = gene.mrna[0].cds[0].length()
        if this_length > length:
            length = this_length
            longest = gene
    if longest:
        transcript.genes = [longest]
