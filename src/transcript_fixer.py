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

def fix_phase(transcript):
    """Changes start indices and phase values for CDSs starting at 2 or 3.

    Adjusts start index for gene, mRNA and CDS to 1 and adds
    appropriate phase value for the CDS; we theorize that this
    is necessary to eliminate errors from the NCBI TSA submission.
    """
    for gene in transcript.genes:
        # Verify we have a valid gene here
        if not gene.mrna:
            return
        if not gene.mrna[0].cds:
            return
        gene_start = gene.start
        mrna_start = gene.mrna[0].start
        cds_start = gene.mrna[0].cds[0].start
        # All three indices should be the same, or this task gets complicated
        if not (gene_start == mrna_start and mrna_start == cds_start):
            return

        # Adjust phase if our feature start on base 2 or 3
        if gene_start == 2:
            gene.start = 1
            gene.mrna[0].start = 1
            gene.mrna[0].cds[0].start = 1
            gene.mrna[0].cds[0].phase = 1
        elif gene_start == 3:
            gene.start = 1
            gene.mrna[0].start = 1
            gene.mrna[0].cds[0].start = 1
            gene.mrna[0].cds[0].phase = 2



