#!/usr/bin/env python

def fix_transcript(transcript):
    for gene in transcript.genes:
        if not hasattr(gene, "mrna"):
            raise Exception("can't fix transcript "+transcript.sequence.header+" because gene "+gene.attributes["ID"]+" has no mRNA")
        if len(gene.mrna) > 1:
            raise Exception("can't fix transcript "+transcript.sequence.header+" because gene "+gene.attributes["ID"]+" has multiple mRNAs")
