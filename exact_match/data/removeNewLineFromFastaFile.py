#!/bin/python

import sys
import random

dbFile = open(sys.argv[1])

line = dbFile.readline()
line = line.strip()
print line

for line in dbFile:
    line = line.strip()
    sys.stdout.write(line)

print ""
