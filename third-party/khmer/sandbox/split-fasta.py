#! /usr/bin/env python2
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
import sys
import screed

filename = sys.argv[1]
prefix = sys.argv[2]
size = int(float(sys.argv[3]))          # e.g. 1e9

division = -1
for n, record in enumerate(screed.open(filename)):
    if n % 100000 == 0:
        print '...', n

    if n % size == 0:
        division += 1
        new_name = '%s.%04d.fa' % (prefix, division)
        print 'opening', new_name
        fp = open(new_name, 'w')

    fp.write('>%s\n%s\n' % (record['name'], record['sequence']))
