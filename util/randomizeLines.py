#!/usr/bin/env python

import sys
import random

lines = sys.stdin.readlines()
random.shuffle(lines)
for line in lines:
    sys.stdout.write(line)
