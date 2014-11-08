#!/usr/bin/env python

from collections import namedtuple
from src.sequence import Sequence
from src.sequtil import get_subsequence, reverse_complement, has_start_codon, has_stop_codon
from src.gff_feature import GFFFeature

def gene_to_tbl(gene):
    # Check for starts and stops
    has_start = False
    has_stop = False
    if hasattr(gene.mrna[0], "start_codon"):
        has_start = True
    if hasattr(gene.mrna[0], "stop_codon"):
        has_stop = True
    # Create tbl entry
    tbl = ""
    if not has_start:
        tbl += "<"
    tbl += str(gene.start)+"\t"
    if not has_stop:
        tbl += ">"
    tbl += str(gene.end)+"\tgene\n"
    # Gene name if it has one
    if "Name" in gene.attributes:
        tbl += "\t\t\tgene\t"+gene.attributes["Name"]+"\n"
    # Locus tag
    tbl += "\t\t\tlocus_tag\t"+gene.attributes["ID"]+"\n"
    if not has_start:
        tbl += "<"
    tbl += str(gene.mrna[0].cds[0].start)+"\t"
    if not has_stop:
        tbl += ">"
    tbl += str(gene.mrna[0].cds[0].end)+"\tCDS\n"
    # Codon start
    if gene.mrna[0].cds[0].phase != 0:
        tbl += "\t\t\tcodon_start\t"+str(gene.mrna[0].cds[0].phase+1)+"\n"
    # Protein id
    tbl += "\t\t\tprotein_id\t"+gene.mrna[0].attributes["ID"]+"\n"
    # Dbxref if it has any
    if "Dbxref" in gene.mrna[0].attributes:
        for dbxref in gene.mrna[0].attributes["Dbxref"].split(","):
            tbl += "\t\t\tdb_xref\t"+dbxref+"\n"
    # product if it has any
    if "product" in gene.mrna[0].attributes:
        tbl += "\t\t\tproduct\t"+gene.mrna[0].attributes["product"]+"\n"
    else: # no product, write hypothetical protein
        tbl += "\t\t\tproduct\thypothetical protein\n"
    # Ontology_term if it has any
    if "Ontology_term" in gene.mrna[0].attributes:
        for term in gene.mrna[0].attributes["Ontology_term"].split(","):
            tbl += "\t\t\tOntology_term\t"+term+"\n"
    return tbl

###################

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

    def fix_phase(self):
        """Changes start indices and phase values for CDSs starting at 2 or 3.

        Adjusts start index for partial gene, mRNA and CDS to 1 and adds
        appropriate phase value for the CDS; we theorize that this
        is necessary to eliminate errors from the NCBI TSA submission.
        """
        for gene in self.genes:
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
                gene.end = len(self.sequence.bases)
                gene.mrna[0].end = len(self.sequence.bases)
                gene.mrna[0].cds[0].end = len(self.sequence.bases)

    def fix_multiple_genes(self):
        print("\n\n\nfoo")
        longest = None
        length = 0
        for gene in self.genes:
            this_length = gene.mrna[0].cds[0].length()
            if this_length > length:
                length = this_length
                longest = gene
        if longest:
            self.genes = [longest]

    def remove_contig_from_gene_id(self):
        for gene in self.genes:
            id_split = gene.attributes['ID'].split('|')
            if len(id_split) == 2:
                gene.attributes['ID'] = id_split[1]

    def make_positive(self):
        if not self.genes or self.genes[0].strand == "+":
            return
        seq_len = len(self.sequence.bases)
        self.sequence.bases = reverse_complement(self.sequence.bases)
        for gene in self.genes:
            gene.start, gene.end = seq_len-gene.end+1, seq_len-gene.start+1
            gene.strand = "+"
            for mrna in gene.mrna:
                mrna.start, mrna.end = seq_len-mrna.end+1, seq_len-mrna.start+1
                mrna.strand = "+"
                for cds in mrna.cds:
                    cds.start, cds.end = seq_len-cds.end+1, seq_len-cds.start+1
                    cds.strand = "+"
                for exon in mrna.exon:
                    exon.start, exon.end = seq_len-exon.end+1, seq_len-exon.start+1
                    exon.strand = "+"

    def fix_feature_lengths(self):
        seq_len = len(self.sequence.bases)
        for gene in self.genes:
            if gene.end > seq_len:
                over = gene.end-seq_len
                gene.end = seq_len-((abs(3-over))%3)
            for mrna in gene.mrna:
                if mrna.end > seq_len:
                    over = mrna.end-seq_len
                    mrna.end = seq_len-((abs(3-over))%3)
                for cds in mrna.cds:
                    if cds.end > seq_len:
                        over = cds.end-seq_len
                        cds.end = seq_len-((abs(3-over))%3)
                for exon in mrna.exon: 
                    if exon.end > seq_len:
                        over = exon.end-seq_len
                        exon.end = seq_len-((abs(3-over))%3)

    def match_cds_and_exon_end(self):
        """Check each mRNA's exon/CDS. If no stop codon, make their ends equal.
        This is a blind attempt to avoid PartialProblem errors from the NCBI."""
        for gene in self.genes:
            for mrna in gene.mrna:
                if hasattr(mrna, "stop_codon"):
                    return
                else:
                    if mrna.cds[0].end != mrna.exon[0].end:
                        mrna.cds[0].end = mrna.exon[0].end

    def create_starts_and_stops(self):
        for gene in self.genes:
            cds = gene.mrna[0].cds[0]
            subseq = get_subsequence(self.sequence.bases, cds.start, cds.end)
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

    def passes_filtering(self):
        if self.rsem == None or self.rsem.tpm < 0.5 or self.rsem.isopct < 5.0:
            return False
        return True

    def to_tbl(self):
        tbl = ""
        tbl += ">Feature "+self.sequence.header+"\n"
        tbl += "1\t"+str(len(self.sequence.bases))+"\tREFERENCE\n"
        tbl += gene_to_tbl(self.genes[0])
        return tbl
