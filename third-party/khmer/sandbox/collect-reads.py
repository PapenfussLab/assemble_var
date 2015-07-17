#! /usr/bin/env python2
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
# pylint: disable=missing-docstring,invalid-name
"""
Build a counting Bloom filter from the given sequences, save in <htname>.
Stop collecting reads when the average coverage gets above -C (default 50).
Place reads into -o output_file.

% collect-reads.py <htname> <data1> [ <data2> <...> ]

Use '-h' for parameter help.
"""

import sys
import textwrap
import khmer
from khmer.khmer_args import build_counting_args, report_on_config, info
from khmer.file import check_file_status, check_space
from khmer.file import check_space_for_hashtable
import argparse
import screed


def output_single(read):
    if hasattr(read, 'accuracy'):
        return "@%s\n%s\n+\n%s\n" % (read.name, read.sequence, read.accuracy)
    else:
        return ">%s\n%s\n" % (read.name, read.sequence)


def get_parser():
    epilog = """
    The memory usage of this script with :option:`-b` will be about
    1.15x the product of the :option:`-x` and :option:`-N` numbers.

    Example::

        collect-reads.py -k 20 -x 5e7 out.kh data/100k-filtered.fa
    """

    parser = build_counting_args("Collect reads until a given avg coverage.",
                                 epilog=textwrap.dedent(epilog))
    parser.add_argument('output_countingtable_filename', help="The name of the"
                        " file to write the k-mer counting table to.")
    parser.add_argument('input_sequence_filename', nargs='+',
                        help="The names of one or more FAST[AQ] input "
                        "sequence files.")
    parser.add_argument('--report-total-kmers', '-t', action='store_true',
                        help="Prints the total number of k-mers to stderr")
    parser.add_argument('-C', '--coverage', type=int,
                        help='Collect reads until this coverage, then exit.')
    parser.add_argument('-o', '--output', type=argparse.FileType('w'),
                        help='Write collect reads into this file.')
    parser.add_argument('-b', '--no-bigcount', dest='bigcount', default=True,
                        action='store_false',
                        help='Do not count k-mers past 255')
    return parser


def main():

    info('collect-reads.py', ['counting'])
    args = get_parser().parse_args()
    report_on_config(args)

    base = args.output_countingtable_filename
    filenames = args.input_sequence_filename

    for name in args.input_sequence_filename:
        check_file_status(name)

    check_space(args.input_sequence_filename)
    check_space_for_hashtable(args.n_tables * args.min_tablesize)

    print 'Saving k-mer counting table to %s' % base
    print 'Loading sequences from %s' % repr(filenames)
    if args.output:
        print 'Outputting sequences to', args.output

    print 'making k-mer counting table'
    htable = khmer.new_counting_hash(args.ksize, args.min_tablesize,
                                     args.n_tables)
    htable.set_use_bigcount(args.bigcount)

    total_coverage = 0.
    n = 0

    for index, filename in enumerate(filenames):
        for record in screed.open(filename):
            seq = record.sequence.upper()
            if 'N' in seq:
                seq = seq.replace('N', 'G')

            try:
                med, _, _ = htable.get_median_count(seq)
            except ValueError:
                continue

            total_coverage += med
            n += 1

            if total_coverage / float(n) > args.coverage:
                print 'reached target average coverage:', \
                      total_coverage / float(n)
                break

            htable.consume(seq)
            if args.output:
                args.output.write(output_single(record))

            if n % 100000 == 0:
                print '...', index, filename, n, total_coverage / float(n)

        if total_coverage / float(n) > args.coverage:
            break

    print 'Collected %d reads' % (n,)

    if args.report_total_kmers:
        print >> sys.stderr, 'Total number of k-mers: {0}'.format(
            htable.n_occupied())

    print 'saving', base
    htable.save(base)

    info_fp = open(base + '.info', 'w')
    info_fp.write('through end: %s\n' % filenames[-1])

    # Change 0.2 only if you really grok it.  HINT: You don't.
    fp_rate = khmer.calc_expected_collisions(htable)
    print 'fp rate estimated to be %1.3f' % fp_rate
    print >> info_fp, 'fp rate estimated to be %1.3f' % fp_rate

    if fp_rate > 0.20:
        print >> sys.stderr, "**"
        print >> sys.stderr, ("** ERROR: the k-mer counting table is too small"
                              " this data set.  Increase tablesize/# tables.")
        print >> sys.stderr, "**"
        sys.exit(1)

    print 'DONE.'

if __name__ == '__main__':
    main()

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
