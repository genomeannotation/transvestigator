#!/usr/bin/env python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4

import sys
from src.gff import read_gff 
from src.annotation import read_annotations, annotate_genes
from src.sequence import read_fasta
from src.transcript_builder import build_transcript_dictionary
from src.transcript_fixer import fix_transcript, fix_phase

fastapath = "transcriptome.fasta"
gffpath = "transcriptome.gff"
annopath = "transcriptome.anno"
tblpath = "transcriptome.new.tbl"
outfastapath = "transcriptome.new.fsa"

def read_transcript_blacklist(io_buffer):
    blacklist = []
    for line in io_buffer:
        transcript = line.strip("\n\t ")
        if transcript:
            blacklist.append(transcript)
    return blacklist

def main():
    path = ""
    if len(sys.argv) > 1:
        path = sys.argv[1] + "/"
    sys.stderr.write("Reading files ... ")
    with open(path + fastapath, "r") as fastafile:
        seqs = read_fasta(fastafile)
    if not seqs:
        sys.stderr.write("Error reading fasta; exiting now\n")
        exit()
    with open(path + gffpath, "r") as gfffile:
        gff = read_gff(gfffile)
    if not gff:
        sys.stderr.write("Error reading gff; exiting now\n")
        exit()
    with open(path + annopath, "r") as annofile:
        annos = read_annotations(annofile)
    if not annos:
        sys.stderr.write("Error reading annotations; exiting now\n")
        exit() # TODO maybe no annos is okay...
    with open(path + "transcript_blacklist", "r") as transcript_blacklist_file:
        transcript_blacklist = read_transcript_blacklist(transcript_blacklist_file)
    if not transcript_blacklist:
        sys.stderr.write("Error reading blacklist; exiting now\n")
        exit()
    sys.stderr.write("done.\n\n")
    sys.stderr.write("Annotating genes ... ")
    annotate_genes(gff.gene, annos)
    sys.stderr.write("done.\n\n")
    sys.stderr.write("Mapping gff data to transcripts ... ")
    transcript_dict = build_transcript_dictionary(seqs, gff.gene)
    sys.stderr.write("done.\n\n")
    sys.stderr.write("Writing .tbl file ... ")
    with open(tblpath, "w") as tblfile:
        for transcript in transcript_dict.values():
            if transcript.sequence.header in transcript_blacklist:
                continue
            transcript.fix_feature_lengths()
            transcript.make_positive()
            transcript.create_starts_and_stops()
            fix_transcript(transcript)
            fix_phase(transcript)
            tblfile.write(transcript.to_tbl())
    sys.stderr.write("done.\n\n")
    sys.stderr.write("Writing .fsa file ... ")
    with open(outfastapath, "w") as outfastafile:
        for transcript in transcript_dict.values():
            if transcript.sequence.header in transcript_blacklist:
                continue
            outfastafile.write(transcript.sequence.to_fasta())
    sys.stderr.write("done.\n\n")


###########################################

if __name__ == '__main__':
    main()
