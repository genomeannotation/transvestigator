#!/usr/bin/env python

def gff_gene_to_tbl(gff_gene):
    if not gff_gene.has_child('mRNA'):
        raise Exception("can't write tbl entry for "+gff_gene.attributes["ID"]+" because it has no mRNAs")
    return ""
