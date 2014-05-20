#!/usr/bin/env python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4

from src.gff_reader import read_gff
from src.gff_feature import GFFFeature
from src.transcript_checker import TranscriptChecker

def main():
    with open("foo.gff", "r") as gfffile:
        gff = read_gff(gfffile)
    if not gff:
        exit()
    checker = TranscriptChecker()
    checker.sort_genes(gff)
    for transcript, geneslist in checker.transcripts.items():
        print(transcript + ": " + str(len(geneslist)))


###########################################

if __name__ == '__main__':
    main()
