#!/usr/bin/env python

from src.gff_feature import GFFFeature
from src.translator import has_start_codon, has_stop_codon
from src.fasta import get_subsequence

def create_starts_and_stops(transcript):
    for gene in transcript.genes:
        if not gene.mrna:
            return
        for mrna in gene.mrna:
            if not mrna.cds:
                return
            begin, end = mrna.cds.start, mrna.cds.end
            cds_seq = get_subsequence(transcript.sequence, begin, end)
            if has_start_codon(cds_seq):
                seqid = mrna.seqid
                source = mrna.source
                type = "start_codon"
                codon_start = begin
                codon_end = begin + 2
                score = None
                strand = mrna.strand
                phase = mrna.phase
                attributes = None # is this okay?
                start_codon = GFFFeature(seqid, source, type, codon_start, codon_end, score,
                                         strand, phase, attributes)
                mrna.add_child(start_codon) 
            if has_stop_codon(cds_seq):
                seqid = mrna.seqid
                source = mrna.source
                type = "stop_codon"
                codon_start = end - 2
                codon_end = end
                score = None
                strand = mrna.strand
                phase = mrna.phase
                attributes = None # is this okay?
                stop_codon = GFFFeature(seqid, source, type, codon_start, codon_end, score,
                                         strand, phase, attributes)
                mrna.add_child(stop_codon) 


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
