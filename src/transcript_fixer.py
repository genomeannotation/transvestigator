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

    Adjusts start index for partial gene, mRNA and CDS to 1 and adds
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

        # Adjust phase if our feature start on base 2 or 3
        if not hasattr(gene.mrna[0], "start_codon"):
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
            if gene.mrna[0].cds[0].start == 2:
                gene.mrna[0].cds[0].start = 1
                gene.mrna[0].cds[0].phase = 1
            elif gene.mrna[0].cds[0].start == 3:
                gene.mrna[0].cds[0].start = 1
                gene.mrna[0].cds[0].phase = 2
        # Adjust end if partial
        if not hasattr(gene.mrna[0], "stop_codon"):
            gene.end = len(transcript.sequence.bases)
            gene.mrna[0].end = len(transcript.sequence.bases)
            gene.mrna[0].cds[0].end = len(transcript.sequence.bases)
