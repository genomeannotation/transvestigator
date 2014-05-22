#!/usr/bin/env python

from collections import namedtuple

Sequence = namedtuple("Sequence", "header bases")
Transcript = namedtuple("Transcript", "genes sequence")

