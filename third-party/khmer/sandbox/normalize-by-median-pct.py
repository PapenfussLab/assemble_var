#! /usr/bin/env python2
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
"""
Eliminate reads with median k-mer abundance higher than
DESIRED_COVERAGE.  Output sequences will be placed in 'infile.keep'.

% python scripts/normalize-by-median.py [ -C <cutoff> ] <data1> <data2> ...

Use '-h' for parameter help.
"""

import sys
import screed
import os
import khmer
from itertools import izip
from khmer.counting_args import build_construct_args, DEFAULT_MIN_HASHSIZE
import argparse

DEFAULT_DESIRED_COVERAGE = 5

# Iterate a collection in arbitrary batches
# from: http://stackoverflow.com/questions/4628290/pairs-from-single-list


def batchwise(t, size):
    it = iter(t)
    return izip(*[it] * size)

# Returns true if the pair of records are properly pairs


def validpair(r0, r1):
    return r0.name[-1] == "1" and \
        r1.name[-1] == "2" and \
        r0.name[0:-1] == r1.name[0:-1]


def main():
    parser = build_construct_args()
    parser.add_argument('-C', '--cutoff', type=int, dest='cutoff',
                        default=DEFAULT_DESIRED_COVERAGE)
    parser.add_argument('-p', '--paired', action='store_true')
    parser.add_argument('-s', '--savehash', dest='savehash', default='')
    parser.add_argument('-l', '--loadhash', dest='loadhash',
                        default='')
    parser.add_argument('-R', '--report-to-file', dest='report_file',
                        type=argparse.FileType('w'))
    parser.add_argument('input_filenames', nargs='+')

    args = parser.parse_args()

    if not args.quiet:
        if args.min_hashsize == DEFAULT_MIN_HASHSIZE and not args.loadhash:
            print>>sys.stderr, "** WARNING: hashsize is default!  You absodefly want to increase this!\n** Please read the docs!"

        print>>sys.stderr, '\nPARAMETERS:'
        print>>sys.stderr, ' - kmer size =    %d \t\t(-k)' % args.ksize
        print>>sys.stderr, ' - n hashes =     %d \t\t(-N)' % args.n_hashes
        print>>sys.stderr, ' - min hashsize = %-5.2g \t(-x)' % args.min_hashsize
        print>>sys.stderr, ' - paired =	      %s \t\t(-p)' % args.paired
        print>>sys.stderr, ''
        print>>sys.stderr, 'Estimated memory usage is %.2g bytes (n_hashes x min_hashsize)' % (
            args.n_hashes * args.min_hashsize)
        print>>sys.stderr, '-' * 8

    K = args.ksize
    HT_SIZE = args.min_hashsize
    N_HT = args.n_hashes
    DESIRED_COVERAGE = args.cutoff
    report_fp = args.report_file
    filenames = args.input_filenames

    # In paired mode we read two records at a time
    batch_size = 1
    if args.paired:
        batch_size = 2

    if args.loadhash:
        print 'loading hashtable from', args.loadhash
        ht = khmer.load_counting_hash(args.loadhash)
    else:
        print 'making hashtable'
        ht = khmer.new_counting_hash(K, HT_SIZE, N_HT)

    total = 0
    discarded = 0

    for input_filename in filenames:
        output_name = os.path.basename(input_filename) + '.keepmedpct'
        outfp = open(output_name, 'w')

        n = -1
        for n, batch in enumerate(batchwise(screed.open(input_filename), batch_size)):
            if n > 0 and n % 100000 == 0:
                print '... kept', total - discarded, 'of', total, ', or', \
                    int(100. - discarded / float(total) * 100.), '%'
                print '... in file', input_filename

                if report_fp:
                    print>>report_fp, total, total - discarded, \
                        1. - (discarded / float(total))
                    report_fp.flush()

            total += batch_size

            # If in paired mode, check that the reads are properly interleaved
            if args.paired:
                if not validpair(batch[0], batch[1]):
                    print >>sys.stderr, 'Error: Improperly interleaved pairs %s %s' % (
                        batch[0].name, batch[1].name)
                    sys.exit(-1)

            # Emit the batch of reads if any read passes the filter
            # and all reads are longer than K
            passed_filter = False
            passed_length = True
            for record in batch:
                if len(record.sequence) < K:
                    passed_length = False
                    continue

                seq = record.sequence.replace('N', 'A')
                med, avg, dev = ht.get_median_count(seq)

                pct = 0.
                if avg:
                    pct = dev / avg * 100

                if med < DESIRED_COVERAGE and pct < 100:
                    ht.consume(seq)
                    passed_filter = True

            # Emit records if any passed
            if passed_length and passed_filter:
                for record in batch:
                    if hasattr(record, 'accuracy'):
                        outfp.write('@%s\n%s\n+\n%s\n' % (record.name,
                                                          record.sequence,
                                                          record.accuracy))
                    else:
                        outfp.write('>%s\n%s\n' %
                                    (record.name, record.sequence))
            else:
                discarded += batch_size

        if -1 < n:
            print 'DONE with', input_filename, '; kept', total - discarded, 'of',\
                total, 'or', int(100. - discarded / float(total) * 100.), '%'
            print 'output in', output_name
        else:
            print 'SKIPPED empty file', input_filename

    if args.savehash:
        print 'Saving hashfile through', input_filename
        print '...saving to', args.savehash
        ht.save(args.savehash)

    # Change 0.2 only if you really grok it.  HINT: You don't.
    fp_rate = khmer.calc_expected_collisions(ht)
    print 'fp rate estimated to be %1.3f' % fp_rate

    if fp_rate > 0.20:
        print >>sys.stderr, "**"
        print >>sys.stderr, "** ERROR: the counting hash is too small for"
        print >>sys.stderr, "** this data set.  Increase hashsize/num ht."
        print >>sys.stderr, "**"
        print >>sys.stderr, "** Do not use these results!!"
        sys.exit(-1)

if __name__ == '__main__':
    main()
