#!/usr/bin/env python

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

    def to_tbl(self):
        tbl = ""
        tbl += ">Feature "+self.sequence.header+"\n"
        tbl += "1\t"+str(len(self.sequence.bases))+"\tREFERENCE\n"
        tbl += "\t\t\tNCBI\t12345\n"
        tbl += gene_to_tbl(self.genes[0])
        return tbl
