#!/usr/bin/env python
import os
import sys
from src.gff import read_gff, write_gff
from src.annotation import read_annotations, annotate_genes
from src.sequence import read_fasta
from src.transcript_builder import build_transcript_dictionary
from src.transcript_fixer import fix_transcript, fix_phase

fastapath = "transcriptome.fasta"  # Required
gffpath = "transcriptome.gff"  # Required
annopath = "transcriptome.anno"  # Optional
blacklistpath = "transcriptome.blacklist"  # Optional
tblpath = "transcriptome.new.tbl"
outgffpath = "transcriptome.new.gff"
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
    annos = None
    transcript_blacklist = None
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
    if verify_path(path + annopath):
        with open(path + annopath, "r") as annofile:
            annos = read_annotations(annofile)
        if not annos:
            sys.stderr.write("Error reading annotations; exiting now\n")
            exit() 
    if verify_path(blacklistpath):
        with open(path + blacklistpath, "r") as transcript_blacklist_file:
            transcript_blacklist = read_transcript_blacklist(transcript_blacklist_file)
        if not transcript_blacklist:
            sys.stderr.write("Error reading blacklist; exiting now\n")
            exit()
    sys.stderr.write("done.\n\n")
    if annos:
        sys.stderr.write("Annotating genes ... ")
        annotate_genes(gff.gene, annos)
        sys.stderr.write("done.\n\n")
    sys.stderr.write("Mapping gff data to transcripts ... ")
    transcript_dict = build_transcript_dictionary(seqs, gff.gene)
    sys.stderr.write("done.\n\n")
    sys.stderr.write("Writing .tbl file ... ")
    with open(tblpath, "w") as tblfile:
        for transcript in transcript_dict.values():
            if transcript_blacklist and\
                    transcript.sequence.header in transcript_blacklist:
                continue
            transcript.fix_feature_lengths()
            transcript.make_positive()
            transcript.create_starts_and_stops()
            fix_transcript(transcript)
            fix_phase(transcript)
            tblfile.write(transcript.to_tbl())
    sys.stderr.write("done.\n\n")
    sys.stderr.write("Writing new .gff file...")
    with open(outgffpath, "w") as outgfffile:
        for transcript in transcript_dict.values():
            if transcript_blacklist and\
                    transcript.sequence.header in transcript_blacklist:
                continue
            write_gff(outgfffile, transcript.genes[0])
    sys.stderr.write("done.\n\n")
    sys.stderr.write("Writing .fsa file ... ")
    with open(outfastapath, "w") as outfastafile:
        for transcript in transcript_dict.values():
            if transcript_blacklist and\
                    transcript.sequence.header in transcript_blacklist:
                continue
            outfastafile.write(transcript.sequence.to_fasta())
    sys.stderr.write("done.\n\n")


def verify_path(path):
    if not os.path.isfile(path):
        return False
    return True


###########################################

if __name__ == '__main__':
    main()
