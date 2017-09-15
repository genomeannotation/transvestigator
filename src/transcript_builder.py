#!/usr/bin/env python
# coding=utf-8

from src.transcript import Transcript 
from src.gene import Gene
from src.mrna import Mrna

def build_transcript_dictionary(seqs, genes):
    transcripts = {}

    for gene in genes:
        gene = Gene.from_gff_feature(gene)
        if gene == None:
            print("Could not convert GFFFeature to Gene, skipping")
            continue

        new_mrnas = []
        for mrna in gene['mrna']:
            mrna = Mrna.from_gff_feature(mrna)
            if mrna == None:
                print("Could not convert GFFFeature to Mrna, skipping")
                continue
            new_mrnas.append(mrna)
        gene.children['mrna'] = new_mrnas

        if gene.seqid in transcripts:
            transcripts[gene.seqid].genes.append(gene)
        else:
            if gene.seqid in seqs:
                transcripts[gene.seqid] = Transcript([gene], seqs[gene.seqid])
            else:
                print("Gene "+gene.attributes["ID"]+" is on sequence "+gene.seqid+" which does not exist. Skipping...")

    return transcripts
