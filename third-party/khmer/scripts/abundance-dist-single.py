#! /usr/bin/env python2
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2010-2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
# pylint: disable=invalid-name,missing-docstring
"""
Produce the k-mer abundance distribution for the given file, without
loading a prebuilt k-mer counting table.

% python scripts/abundance-dist-single.py <data> <histout>

Use '-h' for parameter help.
"""
import os
import sys
import khmer
import threading
import textwrap
from khmer.khmer_args import (build_counting_args, add_threading_args,
                              report_on_config, info)
from khmer.file import (check_file_status, check_space,
                        check_space_for_hashtable)


def get_parser():
    epilog = '''
    Note that with :option:`-b` this script is constant memory; in exchange,
    k-mer counts will stop at 255. The memory usage of this script with
    :option:`-b` will be about 1.15x the product of the :option:`-x` and
    :option:`-N` numbers.

    To count k-mers in multiple files use :program:`load_into_counting.py` and
    :program:`abundance_dist.py`.
    '''
    parser = build_counting_args(
        descr="Calculate the abundance distribution of k-mers from a "
        "single sequence file.", epilog=textwrap.dedent(epilog))
    add_threading_args(parser)

    parser.add_argument('input_sequence_filename', help='The name of the input'
                        ' FAST[AQ] sequence file.')
    parser.add_argument('output_histogram_filename', help='The name of the '
                        'output histogram file. The columns are: (1) k-mer '
                        'abundance, (2) k-mer count, (3) cumulative count, '
                        '(4) fraction of total distinct k-mers.')
    parser.add_argument('-z', '--no-zero', dest='output_zero', default=True,
                        action='store_false',
                        help='Do not output 0-count bins')
    parser.add_argument('-b', '--no-bigcount', dest='bigcount', default=True,
                        action='store_false',
                        help='Do not count k-mers past 255')
    parser.add_argument('-s', '--squash', dest='squash_output', default=False,
                        action='store_true',
                        help='Overwrite output file if it exists')
    parser.add_argument('--savetable', default='', metavar="filename",
                        help="Save the k-mer counting table to the specified "
                        "filename.")
    parser.add_argument('--report-total-kmers', '-t', action='store_true',
                        help="Prints the total number of k-mers to stderr")
    return parser


def main():  # pylint: disable=too-many-locals,too-many-branches
    info('abundance-dist-single.py', ['counting'])
    args = get_parser().parse_args()
    report_on_config(args)

    check_file_status(args.input_sequence_filename)
    check_space([args.input_sequence_filename])
    if args.savetable:
        check_space_for_hashtable(args.n_tables * args.min_tablesize)

    if (not args.squash_output and
            os.path.exists(args.output_histogram_filename)):
        print >> sys.stderr, 'ERROR: %s exists; not squashing.' % \
            args.output_histogram_filename
        sys.exit(1)
    else:
        hist_fp = open(args.output_histogram_filename, 'w')

    print 'making k-mer counting table'
    counting_hash = khmer.new_counting_hash(args.ksize, args.min_tablesize,
                                            args.n_tables,
                                            args.threads)
    counting_hash.set_use_bigcount(args.bigcount)

    print 'building k-mer tracking table'
    tracking = khmer.new_hashbits(counting_hash.ksize(), args.min_tablesize,
                                  args.n_tables)

    print 'kmer_size:', counting_hash.ksize()
    print 'k-mer counting table sizes:', counting_hash.hashsizes()
    print 'outputting to', args.output_histogram_filename

    khmer.get_config().set_reads_input_buffer_size(args.threads * 64 * 1024)

    # start loading
    rparser = khmer.ReadParser(args.input_sequence_filename, args.threads)
    threads = []
    print 'consuming input, round 1 --', args.input_sequence_filename
    for _ in xrange(args.threads):
        thread = \
            threading.Thread(
                target=counting_hash.consume_fasta_with_reads_parser,
                args=(rparser, )
            )
        threads.append(thread)
        thread.start()

    for thread in threads:
        thread.join()

    if args.report_total_kmers:
        print >> sys.stderr, 'Total number of k-mers: {0}'.format(
            counting_hash.n_occupied())

    abundance_lists = []

    def __do_abundance_dist__(read_parser):
        abundances = counting_hash.abundance_distribution_with_reads_parser(
            read_parser, tracking)
        abundance_lists.append(abundances)

    print 'preparing hist from %s...' % args.input_sequence_filename
    rparser = khmer.ReadParser(args.input_sequence_filename, args.threads)
    threads = []
    print 'consuming input, round 2 --', args.input_sequence_filename
    for _ in xrange(args.threads):
        thread = \
            threading.Thread(
                target=__do_abundance_dist__,
                args=(rparser, )
            )
        threads.append(thread)
        thread.start()

    for thread in threads:
        thread.join()

    assert len(abundance_lists) == args.threads, len(abundance_lists)
    abundance = {}
    for abundance_list in abundance_lists:
        for i, count in enumerate(abundance_list):
            abundance[i] = abundance.get(i, 0) + count

    total = sum(abundance.values())

    if 0 == total:
        print >> sys.stderr, \
            "ERROR: abundance distribution is uniformly zero; " \
            "nothing to report."
        print >> sys.stderr, "\tPlease verify that the input files are valid."
        sys.exit(1)

    sofar = 0
    for _, i in sorted(abundance.items()):
        if i == 0 and not args.output_zero:
            continue

        sofar += i
        frac = sofar / float(total)

        print >> hist_fp, _, i, sofar, round(frac, 3)

        if sofar == total:
            break

    if args.savetable:
        print 'Saving k-mer counting table ', args.savetable
        print '...saving to', args.savetable
        counting_hash.save(args.savetable)

    print >> sys.stderr, 'wrote to: ' + args.output_histogram_filename

if __name__ == '__main__':
    main()

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
