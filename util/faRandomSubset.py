#!/usr/bin/env python

import sys
import random
import FASTA
from optparse import OptionParser

usage = "usage: %prog [options] < fastaInput"
optparser = OptionParser(usage)
optparser.add_option("-s", "--seed", type="int", dest="seed",
                     help="seed to random number generator",
                     default=None)
optparser.add_option("-n", type="int", dest="number",
                     help="number of records to randomly output",
                     default=10, metavar="NUMBER")
(options, args) = optparser.parse_args()
if len(args) != 0:
    optparser.error("incorrect number of arguments")

recs = [rec for rec in FASTA.Iterator(sys.stdin)]

random.seed(options.seed)
random.shuffle(recs)
for rec in recs[:options.number]:
    print rec

