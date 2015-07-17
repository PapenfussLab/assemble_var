#! /usr/bin/env python2
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
# pylint: disable=invalid-name, missing-docstring
"""
Partition a graph.

% python scripts/partition-graph.py <base>

This will output many <base>.subset.N.pmap files.

Use '-h' for parameter help.
"""

import threading
import Queue
import gc
import os.path
import argparse
import khmer
from khmer.khmer_args import (add_threading_args, info)
from khmer.file import check_file_status, check_space

# Debugging Support
import re
import platform
if "Linux" == platform.system():
    def __debug_vm_usage(msg):
        print "===> DEBUG: " + msg
        for vmstat in re.findall(r".*Vm.*", file("/proc/self/status").read()):
            print vmstat
else:
    def __debug_vm_usage(msg):  # pylint: disable=unused-argument
        pass

DEFAULT_SUBSET_SIZE = int(1e5)
DEFAULT_N_THREADS = 4


def worker(queue, basename, stop_big_traversals):
    while 1:
        try:
            (htable, index, start, stop) = queue.get(False)
        except Queue.Empty:
            print 'exiting'
            return

        outfile = basename + '.subset.%d.pmap' % (index,)
        if os.path.exists(outfile):
            print 'SKIPPING', outfile, ' -- already exists'
            continue

        print 'starting:', basename, index

        # pay attention to stoptags when partitioning; take command line
        # direction on whether or not to exhaustively traverse.
        subset = htable.do_subset_partition(start, stop, True,
                                            stop_big_traversals)

        print 'saving:', basename, index
        htable.save_subset_partitionmap(subset, outfile)
        del subset
        gc.collect()


def get_parser():
    epilog = """
    The resulting partition maps are saved as '${basename}.subset.#.pmap'
    files.
    """
    parser = argparse.ArgumentParser(
        description="Partition a sequence graph based upon waypoint "
        "connectivity", epilog=epilog,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('basename', help="basename of the input k-mer presence"
                        " table + tagset files")
    parser.add_argument('--stoptags', '-S', metavar='filename', default='',
                        help="Use stoptags in this file during partitioning")
    parser.add_argument('--subset-size', '-s', default=DEFAULT_SUBSET_SIZE,
                        type=float, help='Set subset size (usually 1e5-1e6 is '
                        'good)')
    parser.add_argument('--no-big-traverse', action='store_true',
                        default=False, help='Truncate graph joins at big '
                        'traversals')
    parser.add_argument('--version', action='version', version='%(prog)s '
                        + khmer.__version__)
    add_threading_args(parser)
    return parser


def main():
    info('partition-graph.py', ['graph'])
    args = get_parser().parse_args()
    basename = args.basename

    filenames = [basename + '.pt', basename + '.tagset']
    for _ in filenames:
        check_file_status(_)

    check_space(filenames)

    print '--'
    print 'SUBSET SIZE', args.subset_size
    print 'N THREADS', args.threads
    if args.stoptags:
        print 'stoptag file:', args.stoptags
    print '--'

    print 'loading ht %s.pt' % basename
    htable = khmer.load_hashbits(basename + '.pt')
    htable.load_tagset(basename + '.tagset')

    # do we want to load stop tags, and do they exist?
    if args.stoptags:
        print 'loading stoptags from', args.stoptags
        htable.load_stop_tags(args.stoptags)

    # do we want to exhaustively traverse the graph?
    stop_big_traversals = args.no_big_traverse
    if stop_big_traversals:
        print '** This script brakes for lumps: stop_big_traversals is true.'
    else:
        print '** Traverse all the things: stop_big_traversals is false.'

    #
    # now, partition!
    #

    # divide the tags up into subsets
    divvy = htable.divide_tags_into_subsets(int(args.subset_size))
    n_subsets = len(divvy)
    divvy.append(0)

    # build a queue of tasks:
    worker_q = Queue.Queue()

    # break up the subsets into a list of worker tasks
    for _ in range(0, n_subsets):
        start = divvy[_]
        end = divvy[_ + 1]
        worker_q.put((htable, _, start, end))

    print 'enqueued %d subset tasks' % n_subsets
    open('%s.info' % basename, 'w').write('%d subsets total\n' % (n_subsets))

    n_threads = args.threads
    if n_subsets < n_threads:
        n_threads = n_subsets

    # start threads!
    print 'starting %d threads' % n_threads
    print '---'

    threads = []
    for _ in range(n_threads):
        cur_thrd = threading.Thread(target=worker, args=(worker_q, basename,
                                                         stop_big_traversals))
        threads.append(cur_thrd)
        cur_thrd.start()

    print 'done starting threads'

    # wait for threads
    for _ in threads:
        _.join()

    print '---'
    print 'done making subsets! see %s.subset.*.pmap' % (basename,)

if __name__ == '__main__':
    main()
