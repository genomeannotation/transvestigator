#!/usr/bin/env python

from src.gff_feature import GFFFeature

class GFFException(Exception):
   
    """GFF custom exception class that prints line numbers with error message.
    """

    def __init__(self, line_num, message):
        Exception.__init__(self, "GFF error:"+str(line_num)+": "+message)

###################

def parse_gff_attributes(attr):
    
    """Takes GFF attribute column and returns attributes as a dictionary.
    """

    attr = attr.strip(' \t\n;').split(';') # Sanitize and split
    key_vals = [a.split('=') for a in attr]
    return dict(zip([kv[0] for kv in key_vals], [kv[1] for kv in key_vals]))

###################

def read_gff(io_buffer):

    """Reads a GFF file and returns the root GFFFeature.
    """

    root = GFFFeature()
    features = {} # Dictionary of ID to feature
    orphans = []  # List of orphans

    for line_number, line in enumerate(io_buffer):
        # Skip comments and empty lines
        if not line or line[0] == '#':
            continue

        columns = line.strip(' \t\n').split('\t')
        
        # Build feature
        feature = GFFFeature()
        if columns[0] != '.':
            feature.seqid = columns[0]
        if columns[1] != '.':
            feature.source = columns[1]
        if columns[2] != '.':
            feature.type = columns[2]
        if columns[3] != '.':
            feature.start = int(columns[3])
        if columns[4] != '.':
            feature.end = int(columns[4])
        if columns[5] != '.':
            feature.score = float(columns[5])
        if columns[6] != '.':
            feature.strand = columns[6]
        if columns[7] != '.':
            feature.phase = int(columns[7])
        if columns[8] != '.':
            feature.attributes = parse_gff_attributes(columns[8])

        # Make sure feature has ID
        if not 'ID' in feature.attributes:
            raise GFFException(line_number, "feature has no ID attribute"+str(feature.attributes)+columns[8])

        # Add feature to GFF tree
        if not 'Parent' in feature.attributes: # No parent, add to root
            root.add_child(feature)
        elif feature.attributes['ID'] in features: # Has parent, parent has been created, add to parent
            features[feature.attributes['ID']].add_child(feature)
        else: # Has parent, but it hasn't been created yet. It's an orphan
            orphans.append(feature)

        # Add the feature to our dictionary of features
        features[feature.attributes['ID']] = feature

    # Deal with orphans until there are no more orphans
    while orphans:
        orphan = orphans[0]
        # Check if orphan's parent is in the feature list yet
        if orphan.attributes['ID'] in features:
            features[orphan.attributes['ID']].add_child(orphan)
            orphans.remove(orphan)

    return root
