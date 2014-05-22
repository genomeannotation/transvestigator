#!/usr/bin/env python

def gff_gene_to_tbl(gff_gene):
    if not hasattr(gff_gene, "mrna"):
        raise Exception("can't write tbl entry for "+gff_gene.attributes["ID"]+" because it has no mRNAs")
    if len(gff_gene.mrna) > 1:
        raise Exception("can't write tbl entry for "+gff_gene.attributes["ID"]+" because it has multiple mRNAs")
    if not hasattr(gff_gene.mrna[0], "cds"):
        raise Exception("can't write tbl entry for "+gff_gene.attributes["ID"]+" because its mRNA has no CDS")
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
