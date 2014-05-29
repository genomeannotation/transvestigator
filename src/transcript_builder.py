#!/usr/bin/env python

from src.transcript import Transcript 

def build_transcript_dictionary(seqs, genes):
    transcripts = {}

    for gene in genes:
        if gene.seqid in transcripts:
            transcripts[gene.seqid].genes.append(gene)
        else:
            if gene.seqid in seqs:
                transcripts[gene.seqid] = Transcript([gene], seqs[gene.seqid])
            else:
                print("Gene "+gene.attributes["ID"]+" is on sequence "+gene.seqid+" which does not exist. Skipping...")

    return transcripts
