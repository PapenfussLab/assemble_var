#! /usr/bin/env python2
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
import khmer
import sys
import gc
import glob

K = 32

subset_filenames = sys.argv[1:]

ht = khmer.new_hashbits(K, 1, 1)
for filename in subset_filenames:
    print '--'
    print 'partition map:', filename
    subset = ht.load_subset_partitionmap(filename)
    n_part, n_orphan = ht.subset_count_partitions(subset)
    print 'num partitions:', n_part
    print 'num orphans:', n_orphan

    (dist, n_unassigned) = ht.subset_partition_size_distribution(subset)
    for (size, count) in dist:
        print size, count
    print '%d unassigned tags' % n_unassigned

    print '--'
