#!/usr/bin/env python

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
