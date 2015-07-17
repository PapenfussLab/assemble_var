#! /usr/bin/env python2
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
import khmer
import sys
import screed

K = 32
ht = khmer.new_hashbits(K, 1, 1)

x = [0] * 255
y = [0] * 255

ht.load_stop_tags(sys.argv[1])
for n, record in enumerate(screed.open(sys.argv[2])):
    if n % 10000 == 0:
        sys.stderr.write('... %d\n' % n)

    s, p = ht.trim_on_stoptags(record.sequence)

    if len(s) == len(record.sequence):
        continue

    if p == 0:
        p = 31
    else:
        p += 1

    x[p] += 1
    y[len(record.sequence)] += 1

for i, (n, m) in enumerate(zip(x, y)):
    if m:
        print '%d,%d,%d' % (i, n, m)
