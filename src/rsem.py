# coding=utf-8
from collections import namedtuple

Rsem = namedtuple('Rsem', ['transcript_id', 'gene_id', 'length', 'effective_length', 'expected_count', 'tpm', 'fpkm', 'isopct'])

def read_rsem(io_buffer):
    data = []
    columns = io_buffer.readline()
    for line in io_buffer:
        cols = line.strip().split("\t")
        data.append(Rsem(cols[0], cols[1], int(cols[2]), float(cols[3]), float(cols[4]), float(cols[5]), float(cols[6]), float(cols[7])))
    return data
