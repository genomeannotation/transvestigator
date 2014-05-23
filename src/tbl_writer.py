#!/usr/bin/env python

def gff_gene_to_tbl(gff_gene):
    # Check for starts and stops
    has_start = False
    has_stop = False
    if hasattr(gff_gene.mrna[0], "start_codon"):
        has_start = True
    if hasattr(gff_gene.mrna[0], "stop_codon"):
        has_stop = True
    # Create tbl entry
    tbl = ""
    if not has_start:
        tbl += "<"
    tbl += str(gff_gene.start)+"\t"
    if not has_stop:
        tbl += ">"
    tbl += str(gff_gene.end)+"\tgene\n"
    tbl += "\t\t\tlocus_tag\t"+gff_gene.attributes["ID"]+"\n"
    if not has_start:
        tbl += "<"
    tbl += str(gff_gene.mrna[0].cds[0].start)+"\t"
    if not has_stop:
        tbl += ">"
    tbl += str(gff_gene.mrna[0].cds[0].end)+"\tCDS\n"
    tbl += "\t\t\tprotein_id\t"+gff_gene.mrna[0].attributes["ID"]+"\n"
    return tbl

def transcript_to_tbl(transcript):
    tbl = ""
    tbl += ">Feature "+transcript.sequence.header+"\n"
    tbl += "1\t"+str(len(transcript.sequence.bases))+"\tREFERENCE\n"
    tbl += "\t\t\tNCBI\t12345\n"
    tbl += gff_gene_to_tbl(transcript.genes[0])
    return tbl
