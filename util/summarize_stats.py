#!/usr/bin/env python
# coding=utf-8
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4

import sys

def main():
    if len(sys.argv) != 2:
        sys.stderr.write("usage: summarize_stats.py <file.stats>\n")
        sys.exit()

    total_transcripts = 0
    total_complete = 0
    total_pfam = 0
    total_go = 0
    total_gene_name = 0

    with open(sys.argv[1], 'r') as statsfile:
        for line in statsfile:
            if line.startswith("transcript_id"):
                # skip header line
                continue
            fields = line.strip().split("\t")
            if len(fields) != 5:
                sys.stderr.write("error splitting this line: " + line)
                sys.stderr.write("skipping ...\n")
                continue
            # get values
            complete = fields[1]
            pfam = fields[2]
            go = fields[3]
            gene_name = fields[4]
            # update counts
            total_transcripts += 1
            if complete == "True":
                total_complete += 1
            if pfam == "True":
                total_pfam += 1
            if go == "True":
                total_go += 1
            if gene_name == "True":
                total_gene_name += 1

    complete_percent = (float(total_complete) / total_transcripts) * 100
    pfam_percent = (float(total_pfam) / total_transcripts) * 100
    go_percent = (float(total_go) / total_transcripts) * 100
    gene_name_percent = (float(total_gene_name) / total_transcripts) * 100
    print("Number of transcripts:\t" + str(total_transcripts) + "\n")
    print("Count(%) complete:\t" + str(total_complete) + "(" + str(complete_percent) + ")")
    print("Count(%) with PFAM:\t" + str(total_pfam) + "(" + str(pfam_percent) + ")")
    print("Count(%) with GO:\t" + str(total_go) + "(" + str(go_percent) + ")")
    print("Count(%) with gene name:\t" + str(total_gene_name) + "(" + str(gene_name_percent) + ")")

#####################################################################################

if __name__ == '__main__':
    main()
