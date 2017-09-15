# coding=utf-8
import sys
from collections import namedtuple

BlastHit = namedtuple('BlastHit', ['query_id', 'subject_id', 'percent_identity', 
                        'alignment_length', 'mismatch_count', 'gap_open_count', 
                        'query_start', 'query_end', 'subject_start',
                        'subject_end', 'e_value', 'bit_score'])

def read_blast(io_buffer):
    data = []
    for line in io_buffer:
        cols = line.strip().split()
        if len(cols) < 12:
            sys.stderr.write("Invalid blast line: " + line)
            sys.stderr.write("Skipping... \n")
            continue
        data.append(BlastHit(cols[0], cols[1], float(cols[2]), int(cols[3]), 
                    int(cols[4]), int(cols[5]), int(cols[6]), int(cols[7]),
                    int(cols[8]), int(cols[9]), float(cols[10]), float(cols[11])))
    return data
