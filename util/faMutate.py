#!/usr/bin/python

import FASTA
import random
import sys
from optik import OptionParser

usage = "usage: %prog [options] < faInput > faOutput"
parser = OptionParser(usage)
parser.add_option("-r", "--rate", type="float", dest="rate",
                  help="probability that a given base mutates",
                  default=0.1)
parser.add_option("-s", "--seed", type="int", dest="seed",
                  help="seed for initializing RNG",
                  default=None)
(options, args) = parser.parse_args()
if len(args) != 0:
    parser.error("incorrect number of arguments")

mutations = {'A': ('C', 'G', 'T'),
             'C': ('A', 'G', 'T'),
             'G': ('A', 'C', 'T'),
             'T': ('A', 'C', 'G')}

def mutateCharacter(char, rate):
    if random.random() < rate:
        return random.choice(mutations[char])
    else:
        return char

def mutateSequence(seq, rate):
    return ''.join(map(mutateCharacter, seq, [rate] * len(seq)))

random.seed(options.seed)

for rec in FASTA.Iterator(sys.stdin):
    rec.sequence = mutateSequence(rec.sequence, options.rate)
    print rec
