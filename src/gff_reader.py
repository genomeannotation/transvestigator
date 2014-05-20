#!/usr/bin/env python

def parse_gff_attributes(attr):
    attr = attr.strip(' \t\n;').split(';') # Sanitize and split
    key_vals = [a.split('=') for a in attr]
    return dict(zip([kv[0] for kv in key_vals], [kv[1] for kv in key_vals]))

class GFFReader:
    def __init__(self):
        self.features = []
        self.orphans = []

    def read(self, reader):
        pass        
