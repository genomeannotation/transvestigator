#!/usr/bin/env python

def gff_gene_to_tbl(gff_gene):
    if not gff_gene.has_child("mRNA"):
        raise Exception("can't write tbl entry for "+gff_gene.attributes["ID"]+" because it has no mRNAs")
    if len(gff_gene.mrna) > 1:
        raise Exception("can't write tbl entry for "+gff_gene.attributes["ID"]+" because it has multiple mRNAs")
    if not gff_gene.mrna[0].has_child("CDS"):
        raise Exception("can't write tbl entry for "+gff_gene.attributes["ID"]+" because its mRNA has no CDS")
    tbl = ""
    tbl += "<"+str(gff_gene.start)+"\t>"+str(gff_gene.end)+"\tgene\n"
    tbl += "<"+str(gff_gene.mrna[0].cds[0].start)+"\t>"+str(gff_gene.mrna[0].cds[0].end)+"\tCDS\n"
    tbl += "\t\t\tprotein_id\t"+gff_gene.mrna[0].attributes["ID"]+"\n"
    return tbl
