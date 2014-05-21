#!/usr/bin/env python

def gff_gene_to_tbl(gff_gene):
    if not gff_gene.has_child('mRNA'):
        raise Exception("can't write tbl entry for "+gff_gene.attributes["ID"]+" because it has no mRNAs")
    if len(gff_gene.mrna) > 1:
        raise Exception("can't write tbl entry for "+gff_gene.attributes["ID"]+" because it has multiple mRNAs")
    return ""
