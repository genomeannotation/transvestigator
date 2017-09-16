#!/usr/bin/env python
# coding=utf-8
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4

import sys


def print_usage():
    sys.stderr.write("usage: remove_features.py <remove=file.txt> ")
    sys.stderr.write(" [fasta=file.fasta] ")
    sys.stderr.write("[tbl=file.tbl] [gff=file.gff]")
    sys.stderr.write("[out=prefix] \n")
    sys.stderr.write("file.txt is a list of feature ids to remove\n")


def main():
    if len(sys.argv) < 2:
        print_usage()
        sys.exit()

    # Parse command line arguments
    tbl_filename = ""
    gff_filename = ""
    fasta_filename = ""
    remove_filename = ""
    out_prefix = ""
    for arg in sys.argv[1:]:
        fields = arg.split("=")
        if len(fields) != 2:
            print_usage()
            sys.exit()
        if fields[0] == "remove":
            remove_filename = fields[1]
        if fields[0] == "fasta":
            fasta_filename = fields[1]
        if fields[0] == "tbl":
            tbl_filename = fields[1]
        elif fields[0] == "gff":
            gff_filename = fields[1]
        elif fields[0] == "out":
            out_prefix = fields[1]
    if not remove_filename:
        sys.stderr.write("Error: no 'remove' file specified.\n\n")
        print_usage()
        sys.exit()
    if not out_prefix:
        out_prefix = "features_removed"

    # Read list of features to remove
    to_remove = []
    with open(remove_filename, 'r') as remove_file:
        for line in remove_file:
            to_remove.append(line.strip())

    # Update fasta if provided
    if fasta_filename:
        with open(fasta_filename, 'r') as fasta_file, \
                open(out_prefix + ".fasta", 'w') as fasta_out:
            current_header = ""
            current_seq = ""
            for line in fasta_file:
                if line.startswith(">"):
                    # Header line
                    # Reached the end of a seq
                    if current_seq:
                        fasta_out.write(current_seq + "\n")
                        current_seq = ""
                    header = line.strip()[1:]
                    if header in to_remove:
                        current_header = ""
                    else:
                        current_header = header
                        fasta_out.write(line)
                else:
                    # Sequence line
                    if current_header:
                        current_seq += line.strip()
            # Handle last line
            if current_seq:
                fasta_out.write(current_seq + "\n")

    # Update tbl if provided
    if tbl_filename:
        with open(tbl_filename, 'r') as tbl_file, \
                open(out_prefix + ".tbl", 'w') as tbl_out:
            valid_feature = False
            for line in tbl_file:
                if line.startswith(">"):
                    # Feature line
                    feature_id = line.strip().split()[1]
                    if feature_id in to_remove:
                        valid_feature = False
                    else:
                        valid_feature = True
                if valid_feature:
                    tbl_out.write(line)

    # Update gff if provided
    if gff_filename:
        with open(gff_filename, 'r') as gff_file, \
                open(out_prefix + ".gff", 'w') as gff_out:
            for line in gff_file:
                feature_id = line.strip().split("\t")[0]
                if feature_id not in to_remove:
                    gff_out.write(line)


#####################################################################################

if __name__ == '__main__':
    main()
