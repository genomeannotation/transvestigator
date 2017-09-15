#!/usr/bin/env python

import sys

def read_annotations(io_buffer):
    """Reads input buffer of annotations, returns dictionary of id: [(key, value)].

    Input file should be in the format 'id\tkey\tvalue'.
    """
    error_message = "Error on read_annotations -- skipping invalid line:\n"
    annotations = {}
    for line in io_buffer:
        splitline = line.strip().split("\t")
        if len(splitline) != 3:
            sys.stderr.write(error_message + line + "\n")
            continue
        identifier = splitline[0]
        key = splitline[1]
        value = splitline[2]
        if not validate_key(key):
            sys.stderr.write(error_message + line + "\n")
            continue
        if identifier in annotations:
            annotations[identifier].append((key, value))
        else:
            annotations[identifier] = [(key, value)]
    return annotations

def validate_key(key):
    """Checks annotation key against list of approved keys."""
    permitted_keys = ["name", "Dbxref", "Ontology_term", "product"]
    if key in permitted_keys:
        return True
    return False

def annotate_genes(genes, annotations):
    for gene in genes:
        # TODO do we only annotate mrnas for genes with names?
        if "ID" in gene.attributes:
            gene_id = gene.attributes["ID"]
            if gene_id in annotations:
                gene_anno = annotations[gene_id]
                gene_name = gene_anno[0][1]
                gene.attributes["Name"] = gene_name
        for mrna in gene["mrna"]:
            if "ID" in mrna.attributes:
                mrna_id = mrna.attributes["ID"]
                if mrna_id in annotations:
                    for anno in annotations[mrna_id]:
                        mrna.add_annotation(anno[0], anno[1])
