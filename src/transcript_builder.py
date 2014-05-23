#!/usr/bin/env python

from src.types import Transcript 

def build_transcript_dictionary(seqs, genes):
    transcripts = {}

    for gene in genes:
        if gene.seqid in transcripts:
            transcripts[gene.seqid].genes.append(gene)
        else:
            transcripts[gene.seqid] = Transcript([gene], seqs[gene.seqid])

    return transcripts
