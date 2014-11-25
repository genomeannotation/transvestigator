#!/usr/bin/env python

from collections import namedtuple
from src.sequence import Sequence
from src.sequtil import reverse_complement
from src.gff_feature import GFFFeature

Rsem = namedtuple('Rsem', ['tpm', 'fpkm', 'isopct'])

###################

class Transcript:

    def __init__(self, genes=None, sequence=None):
        if not genes:
            self.genes = []
        else:
            self.genes = genes
        if not sequence:
            self.sequence = Sequence()
        else:
            self.sequence = sequence
        self.rsem = None

    def get_gene(self):
        return self.genes[0]

    def fix_phase(self):
        for gene in self.genes:
            gene.fix_phase(self.sequence.bases)

    def fix_multiple_genes(self):
        longest = None
        length = 0
        for gene in self.genes:
            this_length = gene.get_cds_length()
            if this_length > length:
                length = this_length
                longest = gene
        if longest:
            self.genes = [longest]

    def remove_contig_from_gene_id(self):
        for gene in self.genes:
            gene.remove_contig_from_gene_id()

    def make_positive(self):
        """Only really works if transcript contains a single gene;
        otherwise results not guaranteed."""
        if not self.genes or self.genes[0].strand == "+":
            return
        seq_len = len(self.sequence.bases)
        self.sequence.bases = reverse_complement(self.sequence.bases)
        for gene in self.genes:
            gene.make_positive(seq_len)

    def fix_feature_lengths(self):
        seq_len = len(self.sequence.bases)
        for gene in self.genes:
            gene.fix_feature_lengths(seq_len)

    def match_cds_and_exon_end(self):
        """Check each mRNA's exon/CDS. If no stop codon, make their ends equal.
        This is a blind attempt to avoid PartialProblem errors from the NCBI."""
        for gene in self.genes:
            gene.match_cds_and_exon_end()

    def create_starts_and_stops(self):
        for gene in self.genes:
            gene.create_starts_and_stops(self.sequence.bases)

    def passes_filtering(self):
        if self.rsem and (self.rsem.tpm < 0.5 or self.rsem.isopct < 5.0):
            return False
        return True

    def to_tbl(self):
        tbl = ""
        tbl += ">Feature "+self.sequence.header+"\n"
        tbl += "1\t"+str(len(self.sequence.bases))+"\tREFERENCE\n"
        tbl += self.genes[0].to_tbl()
        return tbl
