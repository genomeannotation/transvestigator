#!/usr/bin/env python

import sys
from collections import namedtuple

Annotation = namedtuple('Annotation', 'id key value')

def validate_annotation(anno):
    """Checks that no field is blank and that the 'key' entry on a list of allowed keys."""

    for field in anno:
        if not field:
            return False

    # The following list can be modified as we support more annotations,
    # but for now let's not add annotations like "Zub=foo:skittles"
    permitted_keys = ["name", "Dbxref", "Ontology_term", "product"]
    if anno.key not in permitted_keys:
        return False
    return True

def read_annotations(io_buffer):
    error_message = ("Error on read_annotations-- unable to turn"
                    " the following line into an annotation (will skip line):\n")
    annotations = []
    for line in io_buffer:
        splitline = line.strip().split("\t")
        if len(splitline) != 3:
            sys.stderr.write(error_message + line + "\n")
            continue
        anno = Annotation(splitline[0], splitline[1], splitline[2])
        if not validate_annotation(anno):
            sys.stderr.write(error_message + line + "\n")
            continue
        annotations.append(anno)
    return annotations
