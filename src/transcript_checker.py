#!/usr/bin/env python

from src.gff_feature import GFFFeature
from src.sequtil import has_start_codon, has_stop_codon, reverse_complement
from src.sequtil import get_subsequence



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
