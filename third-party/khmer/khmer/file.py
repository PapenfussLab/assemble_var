#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#

'''
File handling/checking utilities for command-line scripts.
'''

import os
import sys


def check_file_status(file_path):
    """
    Check status of file - return if file exists; warn and exit
    if empty, or does not exist
    """
    if not os.path.exists(file_path):
        print >>sys.stderr, "ERROR: Input file %s does not exist; exiting" % \
                            file_path
        sys.exit(1)
    else:
        if os.stat(file_path).st_size == 0:
            print >>sys.stderr, "ERROR: Input file %s is empty; exiting." % \
                                file_path
            sys.exit(1)


def check_space(in_files, _testhook_free_space=None):
    """
    Estimate size of input files passed, then calculate
    disk space available. Exit if insufficient disk space,
    """

    # Get disk free space in Bytes assuming non superuser
    # and assuming all inFiles are in same disk
    in_file = in_files[0]

    dir_path = os.path.dirname(os.path.realpath(in_file))
    target = os.statvfs(dir_path)

    if _testhook_free_space is None:
        free_space = target.f_frsize * target.f_bavail
    else:
        free_space = _testhook_free_space

    # Check input file array, remove corrupt files
    valid_files = [f for f in in_files if os.path.isfile(f)]

    # Get input file size as worst case estimate of
    # output file size
    file_sizes = [os.stat(f).st_size for f in valid_files]
    total_size = sum(file_sizes)

    size_diff = total_size - free_space
    if size_diff > 0:
        print >>sys.stderr, "ERROR: Not enough free space on disk " \
                            "for output files;\n" \
                            "       Need at least %.1f GB more." \
                            % (float(size_diff) / 1e9)
        print >>sys.stderr, "       Estimated output size: %.1f GB" \
                            % (float(total_size) / 1e9,)
        print >>sys.stderr, "       Free space: %.1f GB" \
                            % (float(free_space) / 1e9,)
        sys.exit(1)


def check_space_for_hashtable(hash_size, _testhook_free_space=None):
    """
    Check we have enough size to write a hash table
    """
    cwd = os.getcwd()
    dir_path = os.path.dirname(os.path.realpath(cwd))
    target = os.statvfs(dir_path)
    if _testhook_free_space is None:
        free_space = target.f_frsize * target.f_bavail
    else:
        free_space = _testhook_free_space  # allow us to test this code...

    size_diff = hash_size - free_space
    if size_diff > 0:
        print >>sys.stderr, "ERROR: Not enough free space on disk " \
                            "for saved table files;" \
                            "       Need at least %s GB more." \
                            % (float(size_diff) / 1e9,)
        print >>sys.stderr, "       Table size: %.1f GB" \
                            % (float(hash_size) / 1e9,)
        print >>sys.stderr, "       Free space: %.1f GB" \
                            % (float(free_space) / 1e9,)
        sys.exit(1)


def check_valid_file_exists(in_files):
    """
    In a scenario where we expect multiple input files and
    are OK with some of them being empty or non-existent,
    this check warns to stderr if any input file is empty
    or non-existent
    """
    for in_file in in_files:
        if os.path.exists(in_file):
            if os.stat(in_file).st_size > 0:
                return
            else:
                print >>sys.stderr, 'WARNING: Input file %s is empty' % \
                                    in_file
        else:
            print >>sys.stderr, 'WARNING: Input file %s not found' % \
                                in_file
