#!/usr/bin/env python

from src.gff_feature import GFFFeature
from src.sequtil import has_start_codon, has_stop_codon, reverse_complement
from src.sequtil import get_subsequence

def create_starts_and_stops(transcript):
    for gene in transcript.genes:
        cds = gene.mrna[0].cds[0]
        subseq = get_subsequence(transcript.sequence.bases, cds.start, cds.end)
        if cds.strand == '-':
            subseq = reverse_complement(subseq)
        if has_start_codon(subseq):
            seqid = cds.seqid
            source = cds.source
            type = "start_codon"
            codon_start = cds.start
            codon_end = cds.start + 2
            score = None
            strand = cds.strand
            phase = cds.phase
            mrna_id = gene.mrna[0].attributes["ID"]
            attributes = {"ID": mrna_id+":start", "Parent": mrna_id}
            start_codon = GFFFeature(seqid, source, type, codon_start, codon_end, score,
                                    strand, phase, attributes)
            gene.mrna[0].add_child(start_codon)
        if has_stop_codon(subseq):
            seqid = cds.seqid
            source = cds.source
            type = "stop_codon"
            codon_start = cds.end - 2
            codon_end = cds.start
            score = None
            strand = cds.strand
            phase = cds.phase
            mrna_id = gene.mrna[0].attributes["ID"]
            attributes = {"ID": mrna_id+":stop", "Parent": mrna_id}
            stop_codon = GFFFeature(seqid, source, type, codon_start, codon_end, score,
                                    strand, phase, attributes)
            gene.mrna[0].add_child(stop_codon)


class TranscriptChecker:

    def __init__(self):
        self.transcripts = {}

    def overlap(self, indices1, indices2):
        # Case 1:
        #   indices1    --------
        #   indices2  ------
        if indices1[0] > indices2[0] and indices1[0] <= indices2[1]:
            return True
        # Case 2:
        #   indices1  --------
        #   indices2      ------
        elif indices1[1] >= indices2[0] and indices1[1] < indices2[1]:
            return True
        else:
            return False

    def nested(self, indices1, indices2):
        if indices1[0] >= indices2[0] and indices1[1] <= indices2[1]:
            return True
        elif indices2[0] >= indices1[0] and indices2[1] <= indices1[1]:
            return True
        else:
            return False
        
    
    def sort_genes(self, gff):
        for gene in gff.gene:
            if gene.seqid in self.transcripts:
                self.transcripts[gene.seqid].append(gene)
            else:
                self.transcripts[gene.seqid] = [gene]
