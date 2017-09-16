#!/usr/bin/env python
# coding=utf-8


class GFFFeature:
    def __init__(self, seqid=None, source=None, feature_type=None, start=None, end=None,
                 score=None, strand=None, phase=0, attributes=None, children=None):
        self.seqid = seqid
        self.source = source
        self.feature_type = feature_type
        self.start = start
        self.end = end
        self.score = score
        self.strand = strand
        self.phase = phase
        if not attributes:
            self.attributes = {}  # Empty dictionary
        else:
            self.attributes = attributes
        if children is None:
            self.children = {}
        else:
            self.children = children

    def add_child(self, child):
        child_type = child.type.lower()
        if child.type in self.children:
            self.children[child_type].append(child)
        else:
            self.children[child_type] = [child]

    def __getitem__(self, index):
        return self.children[index]

    def __setitem__(self, index, value):
        self.children[index] = value

    def __iter__(self):
        return self.children.__iter__()

    def __contains__(self, item):
        return item in self.children

    def add_annotation(self, key, value):
        if key in self.attributes:
            self.attributes[key] += "," + value
        else:
            self.attributes[key] = value

    def length(self):
        return self.end - self.start + 1

    def fix_feature_lengths(self, seq_len):
        """Trims sequences that extend beyond the end of a sequence, ensuring to only
        remove full codons
        """
        if self.end > seq_len:
            over = self.end - seq_len
            self.end = seq_len - ((abs(3 - over)) % 3)
        for features in self.children.values():
            for feature in features:
                feature.fix_feature_lengths(seq_len)

    def write(self):
        out = ""
        if self.seqid is not None:
            out += self.seqid + "\t"
        else:
            out += ".\t"
        if self.source is not None:
            out += self.source + "\t"
        else:
            out += ".\t"
        if self.feature_type is not None:
            out += self.feature_type + "\t"
        else:
            out += ".\t"
        if self.start is not None:
            out += str(self.start) + "\t"
        else:
            out += ".\t"
        if self.end is not None:
            out += str(self.end) + "\t"
        else:
            out += ".\t"
        if self.score is not None:
            out += str(self.score) + "\t"
        else:
            out += ".\t"
        if self.strand is not None:
            out += self.strand + "\t"
        else:
            out += ".\t"
        if self.phase is not None:
            out += str(self.phase) + "\t"
        else:
            out += ".\t"
        out += ";".join([key + "=" + val for key, val in sorted(self.attributes.items())])
        return out
