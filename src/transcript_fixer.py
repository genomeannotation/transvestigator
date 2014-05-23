#!/usr/bin/env python

from src.gff_feature import gff_feature_length

def fix_transcript(transcript):
    longest = None
    length = 0
    for gene in transcript.genes:
        this_length = gff_feature_length(gene.mrna[0].cds[0])
        if this_length > length:
            length = this_length
            longest = gene
    if longest:
        transcript.genes = [longest]
