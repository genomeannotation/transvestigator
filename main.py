#!/usr/bin/env python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4

from src.gff_reader import read_gff
from src.gff_feature import GFFFeature
from src.transcript_checker import TranscriptChecker

gffpath = "foo.gff"
annopath = "foo.anno"

def main():
    with open(gffpath, "r") as gfffile:
        gff = read_gff(gfffile)
    if not gff:
        exit()
    with open(annopath, "r") as annofile:
        annos = read_ipr(annofile)



###########################################

if __name__ == '__main__':
    main()
