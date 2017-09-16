#!/usr/bin/env python
# coding=utf-8

from collections import namedtuple

from src.sequence import Sequence
from src.sequtil import reverse_complement

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
        self.MAXIMUM_GENE_LENGTH = 50  # I made this up :/

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

    def name_too_long(self, gene_name):
        return len(gene_name) >= self.MAXIMUM_GENE_LENGTH

    def fix_long_gene_names(self):
        for gene in self.genes:
            if "Name" in gene.attributes:
                gene_name = gene.attributes["Name"]
                if self.name_too_long(gene_name):
                    del gene.attributes["Name"]

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

    def stats(self):
        """Returns a string of tab-separated boolean values:
        transcript_id\tcomplete\tpfam_domain\tgo_annotation\tgene_name
        """
        # add transcript id
        result = self.sequence.header + "\t"
        # find other values
        complete = True
        pfam = False
        go = False
        gene_name = False
        for gene in self.genes:
            if not gene.is_complete():
                complete = False
            if "Name" in gene.attributes:
                gene_name = True
            if "Dbxref" in gene.get_mrna().attributes:
                for dbxref in gene.get_mrna().attributes["Dbxref"].split(","):
                    if "PFAM" in dbxref:
                        pfam = True
                    if "GO" in dbxref:
                        go = True
        # append other values to result string
        result += str(complete) + "\t"
        result += str(pfam) + "\t"
        result += str(go) + "\t"
        result += str(gene_name)
        result += "\n"
        return result

    def to_tbl(self):
        tbl = ""
        tbl += ">Feature " + self.sequence.header + "\n"
        tbl += "1\t" + str(len(self.sequence.bases)) + "\tREFERENCE\n"
        tbl += self.genes[0].to_tbl()
        return tbl
