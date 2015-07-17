#! /usr/bin/env python2
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
import sys
import screed
import gzip
import os.path

next_partition = 2
filenum = 0
for filename in sys.argv[1:]:
    filenum += 1
    outfp = gzip.open('group%03d.fa.gz' % filenum, 'w')

    old_to_new = {}
    for n, record in enumerate(screed.open(filename)):
        if n > 0 and n % 10000 == 0:
            print '...', os.path.basename(filename), n
        partition = record.name.split()[-1]
        name = record.name.split()[0]

        new_part = old_to_new.get(partition)
        if new_part is None:
            new_part = next_partition
            next_partition += 1
            old_to_new[partition] = new_part

        outfp.write('>%s\t%d\n%s\n' % (name, new_part, record.sequence))
    outfp.close()
    print 'renumbered %d partitions in %s' % (len(old_to_new), filename)
