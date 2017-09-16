#!/usr/bin/env python
# coding=utf-8
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4

import sys

MAX_GENE_NAME_LENGTH = 30


def print_usage():
    sys.stderr.write("usage: replace_long_gene_names.py ")
    sys.stderr.write("[tbl=file.tbl] [gff=file.gff]")
    sys.stderr.write("[out=prefix] \n")


def main():
    if len(sys.argv) < 2:
        print_usage()
        sys.exit()

    # Parse command line arguments
    tbl_filename = ""
    gff_filename = ""
    out_prefix = ""
    for arg in sys.argv[1:]:
        fields = arg.split("=")
        if len(fields) != 2:
            print_usage()
            sys.exit()
        if fields[0] == "tbl":
            tbl_filename = fields[1]
        elif fields[0] == "gff":
            gff_filename = fields[1]
        elif fields[0] == "out":
            out_prefix = fields[1]
    if not out_prefix:
        out_prefix = "long_gene_names_removed"

    # Update tbl if provided
    if tbl_filename:
        with open(tbl_filename, 'r') as tbl_file, \
                open(out_prefix + ".tbl", 'w') as tbl_out:
            for line in tbl_file:
                fields = line.split("\t")
                if len(fields) >= 5 and fields[3] == "gene" and len(fields[4]) > MAX_GENE_NAME_LENGTH:
                    continue
                else:
                    tbl_out.write(line)

    # Update gff if provided
    if gff_filename:
        # noinspection PyUnusedLocal
        with open(gff_filename, 'r') as gff_file, \
                open(out_prefix + ".gff", 'w') as gff_out:
            # noinspection PyUnusedLocal
            for line in gff_file:
                # ToDo:
                pass


#####################################################################################

if __name__ == '__main__':
    main()
