#!/usr/bin/env python

import sys
import random
import FASTA

for rec in FASTA.Iterator(sys.stdin):
    chars = list(rec.sequence)
    random.shuffle(chars)
    rec.sequence = "".join(chars)
    print rec
