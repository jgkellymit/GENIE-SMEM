#!/bin/python

import sys
import random

dbFile = open(sys.argv[1])

dbFile.readline()
ind = 0

n = int(sys.argv[2])
length = int(sys.argv[3])

dna = ""
for line in dbFile:
    dna = dna + line.strip()
dnalen = len(dna)

random.seed(0)
for i in range(n):
    ind = random.randint(0, dnalen - length)

    print ">" + sys.argv[1] + "_" + str(ind)
    print dna[ind:ind + length]
