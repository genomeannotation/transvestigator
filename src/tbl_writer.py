#!/usr/bin/env python

def gff_gene_to_tbl(gff_gene):
    if not gff_gene.mrna:
        raise Exception("can't write tbl for "+gff_gene.attributes["ID"]+" because it has no mRNAs")
    return ""
