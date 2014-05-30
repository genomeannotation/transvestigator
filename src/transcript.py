#!/usr/bin/env python

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
    tbl += "\t\t\tprotein_id\t"+gene.mrna[0].attributes["ID"]+"\n"
    # Dbxref if it has any
    if "Dbxref" in gene.mrna[0].attributes:
        for dbxref in gene.mrna[0].attributes["Dbxref"].split(","):
            tbl += "\t\t\tdb_xref\t"+dbxref+"\n"
    return tbl

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

    def to_tbl(self):
        tbl = ""
        tbl += ">Feature "+self.sequence.header+"\n"
        tbl += "1\t"+str(len(self.sequence.bases))+"\tREFERENCE\n"
        tbl += "\t\t\tNCBI\t12345\n"
        tbl += gene_to_tbl(self.genes[0])
        return tbl