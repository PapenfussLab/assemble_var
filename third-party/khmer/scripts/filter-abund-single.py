#! /usr/bin/env python2
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
# pylint: disable=missing-docstring,invalid-name
"""
Trim sequences at k-mers of the given abundance for the given file,
without loading a prebuilt counting table.  Output sequences will be
placed in 'infile.abundfilt'.

% python scripts/filter-abund-single.py <data>

Use '-h' for parameter help.
"""
import os
import sys
import khmer
import threading
import textwrap
from khmer.thread_utils import ThreadedSequenceProcessor, verbose_loader
from khmer.khmer_args import (build_counting_args, report_on_config,
                              add_threading_args, info)
from khmer.file import (check_file_status, check_space,
                        check_space_for_hashtable)
#
DEFAULT_CUTOFF = 2


def get_parser():
    epilog = """
    Trimmed sequences will be placed in ${input_sequence_filename}.abundfilt.

    This script is constant memory.

    To trim reads based on k-mer abundance across multiple files, use
    :program:`load-into-counting.py` and :program:`filter-abund.py`.

    Example::

        filter-abund-single.py -k 20 -x 5e7 -C 2 data/100k-filtered.fa
    """
    parser = build_counting_args(
        descr="Trims sequences at a minimum k-mer abundance "
        "(in memory version).", epilog=textwrap.dedent(epilog))
    add_threading_args(parser)

    parser.add_argument('--cutoff', '-C', default=DEFAULT_CUTOFF, type=int,
                        help="Trim at k-mers below this abundance.")
    parser.add_argument('--savetable', metavar="filename", default='',
                        help="If present, the name of the file to save the "
                        "k-mer counting table to")
    parser.add_argument('datafile', metavar='input_sequence_filename',
                        help="FAST[AQ] sequence file to trim")
    parser.add_argument('--report-total-kmers', '-t', action='store_true',
                        help="Prints the total number of k-mers to stderr")
    return parser


def main():
    info('filter-abund-single.py', ['counting'])
    args = get_parser().parse_args()
    check_file_status(args.datafile)
    check_space([args.datafile])
    if args.savetable:
        check_space_for_hashtable(args.n_tables * args.min_tablesize)
    report_on_config(args)

    config = khmer.get_config()
    config.set_reads_input_buffer_size(args.threads * 64 * 1024)

    print 'making k-mer counting table'
    htable = khmer.new_counting_hash(args.ksize, args.min_tablesize,
                                     args.n_tables,
                                     args.threads)

    # first, load reads into hash table
    rparser = khmer.ReadParser(args.datafile, args.threads)
    threads = []
    print 'consuming input, round 1 --', args.datafile
    for _ in xrange(args.threads):
        cur_thread = \
            threading.Thread(
                target=htable.consume_fasta_with_reads_parser,
                args=(rparser, )
            )
        threads.append(cur_thread)
        cur_thread.start()

    for _ in threads:
        _.join()

    if args.report_total_kmers:
        print >> sys.stderr, 'Total number of k-mers: {0}'.format(
            htable.n_occupied())

    fp_rate = khmer.calc_expected_collisions(htable)
    print 'fp rate estimated to be %1.3f' % fp_rate

    # now, trim.

    # the filtering function.
    def process_fn(record):
        name = record['name']
        seq = record['sequence']
        if 'N' in seq:
            return None, None

        trim_seq, trim_at = htable.trim_on_abundance(seq, args.cutoff)

        if trim_at >= args.ksize:
            return name, trim_seq

        return None, None

    # the filtering loop
    print 'filtering', args.datafile
    outfile = os.path.basename(args.datafile) + '.abundfilt'
    outfp = open(outfile, 'w')

    tsp = ThreadedSequenceProcessor(process_fn)
    tsp.start(verbose_loader(args.datafile), outfp)

    print 'output in', outfile

    if args.savetable:
        print 'Saving k-mer counting table filename', args.savetable
        print '...saving to', args.savetable
        htable.save(args.savetable)
    print('wrote to: ' + outfile)

if __name__ == '__main__':
    main()
