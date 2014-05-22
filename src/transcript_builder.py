#!/usr/bin/env python

from src.types import Transcript 

def build_transcript_dictionary(seqs, genes):
    transcripts = {}

    for seq in seqs:
        for gene in genes:
            if seq.header == gene.seqid:
                if seq.header in transcripts:
                    transcripts[seq.header].genes.append(gene)
                else:
                    transcripts[seq.header] = Transcript([gene], seq)

    return transcripts
