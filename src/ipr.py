#!/usr/bin/env python

def read_ipr(io_buffer):
    """Returns a dictionary of feature ids to lists, each containing "Dbxref" and annotation."""
    iprs = {}
    for line in io_buffer:
        columns = line.split("\t")
        if len(columns)>1:
            if columns[0] in iprs:
                iprs[columns[0]].append(["Dbxref", columns[3]+":"+columns[4]])
            else:
                iprs[columns[0]] = [["Dbxref", columns[3]+":"+columns[4]]]
    return iprs
