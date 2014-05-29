#!/usr/bin/env python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4

from src.gff import read_gff, annotate_genes
from src.sequence import read_fasta
from src.ipr import read_ipr
from src.transcript_builder import build_transcript_dictionary
from src.transcript_fixer import fix_transcript

fastapath = "foo.fasta"
gffpath = "foo.gff"
annopath = "foo.anno"
tblpath = "foo.tbl"

def main():
    with open(fastapath, "r") as fastafile:
        seqs = read_fasta(fastafile)
    if not seqs:
        exit()
    with open(gffpath, "r") as gfffile:
        gff = read_gff(gfffile)
    if not gff:
        exit()
    with open(annopath, "r") as annofile:
        annos = read_ipr(annofile)
    if not annos:
        exit() # TODO maybe no annos is okay...
    annotate_genes(gff.gene, annos)
    transcript_dict = build_transcript_dictionary(seqs, gff.gene)
    with open(tblpath, "w") as tblfile:
        for transcript in transcript_dict.values():
            transcript.create_starts_and_stops()
            fix_transcript(transcript)
            tblfile.write(transcript.to_tbl())







###########################################

if __name__ == '__main__':
    main()
