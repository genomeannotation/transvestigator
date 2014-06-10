#!/usr/bin/env python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4

from src.gff import read_gff 
from src.annotation import read_annotations, annotate_genes
from src.sequence import read_fasta
from src.ipr import read_ipr
from src.transcript_builder import build_transcript_dictionary
from src.transcript_fixer import fix_transcript

fastapath = "foo.fasta"
gffpath = "foo.gff"
annopath = "foo.anno"
tblpath = "foo.tbl"
outfastapath = "foo.out.fsa"

def read_transcript_blacklist(io_buffer):
    blacklist = []
    for line in io_buffer:
        transcript = line.strip("\n\t ")
        if transcript:
            blacklist.append(transcript)
    return blacklist

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
        annos = read_annotations(annofile)
    if not annos:
        exit() # TODO maybe no annos is okay...
    with open("transcript_blacklist", "r") as transcript_blacklist_file:
        transcript_blacklist = read_transcript_blacklist(transcript_blacklist_file)
    if not transcript_blacklist:
        exit()
    annotate_genes(gff.gene, annos)
    transcript_dict = build_transcript_dictionary(seqs, gff.gene)
    with open(tblpath, "w") as tblfile:
        for transcript in transcript_dict.values():
            if transcript.sequence.header in transcript_blacklist:
                continue
            transcript.fix_feature_lengths()
            transcript.make_positive()
            transcript.create_starts_and_stops()
            fix_transcript(transcript)
            tblfile.write(transcript.to_tbl())
    with open(outfastapath, "w") as outfastafile:
        for transcript in transcript_dict.values():
            outfastafile.write(transcript.sequence.to_fasta())


###########################################

if __name__ == '__main__':
    main()
