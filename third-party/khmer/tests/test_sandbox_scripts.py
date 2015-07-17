#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#

# pylint: disable=C0111,C0103,E1103,W0612

import sys
import os
import shutil
from cStringIO import StringIO
import traceback
import nose

import khmer_tst_utils as utils
import khmer
import khmer.file
import screed


def scriptpath(script):
    return script


def teardown():
    utils.cleanup()


def test_sweep_reads():
    readfile = utils.get_temp_filename('reads.fa')
    contigfile = utils.get_temp_filename('contigs.fp')
    in_dir = os.path.dirname(contigfile)

    shutil.copyfile(utils.get_test_data('test-sweep-reads.fa'), readfile)
    shutil.copyfile(utils.get_test_data('test-sweep-contigs.fp'), contigfile)

    script = scriptpath('sweep-reads.py')
    args = ['-k', '25', '--prefix', 'test', '--label-by-pid',
            contigfile, readfile, 'junkfile.fa']

    status, out, err = utils.runscript(
        script, args, in_dir, fail_ok=True, sandbox=True)

    # check if the bad file was skipped without issue
    assert 'ERROR' in err, err
    assert 'skipping' in err, err

    out1 = os.path.join(in_dir, 'test_0.fa')
    out2 = os.path.join(in_dir, 'test_1.fa')
    mout = os.path.join(in_dir, 'test_multi.fa')
    oout = os.path.join(in_dir, 'test_orphaned.fa')

    print os.listdir(in_dir)

    seqs1 = set([r.name for r in screed.open(out1)])
    seqs2 = set([r.name for r in screed.open(out2)])
    seqsm = set([r.name for r in screed.open(mout)])
    seqso = set([r.name for r in screed.open(oout)])

    print seqs1
    print seqs2
    print seqsm
    print seqso
    assert seqs1 == set(['read1_p0\t0', 'read2_p0\t0'])
    assert seqs2 == set(['read3_p1\t1'])
    assert (seqsm == set(['read4_multi\t0\t1']) or
            seqsm == set(['read4_multi\t1\t0']))
    assert seqso == set(['read5_orphan'])


def test_sweep_reads_fq():
    readfile = utils.get_temp_filename('reads.fa')
    contigfile = utils.get_temp_filename('contigs.fp')
    in_dir = os.path.dirname(contigfile)

    shutil.copyfile(utils.get_test_data('test-sweep-reads.fq'), readfile)
    shutil.copyfile(utils.get_test_data('test-sweep-contigs.fp'), contigfile)

    script = scriptpath('sweep-reads.py')
    args = ['-k', '25', '--prefix', 'test', '--label-by-pid',
            contigfile, readfile, 'junkfile.fa']

    status, out, err = utils.runscript(
        script, args, in_dir, fail_ok=True, sandbox=True)

    # check if the bad file was skipped without issue
    assert 'ERROR' in err, err
    assert 'skipping' in err, err

    out1 = os.path.join(in_dir, 'test_0.fq')
    out2 = os.path.join(in_dir, 'test_1.fq')
    mout = os.path.join(in_dir, 'test_multi.fq')
    oout = os.path.join(in_dir, 'test_orphaned.fq')

    print open(out1).read()

    print os.listdir(in_dir)

    seqs1 = set([r.name for r in screed.open(out1)])
    seqs2 = set([r.name for r in screed.open(out2)])
    seqsm = set([r.name for r in screed.open(mout)])
    seqso = set([r.name for r in screed.open(oout)])

    print seqs1
    print seqs2
    print seqsm
    print seqso
    assert seqs1 == set(['read1_p0\t0', 'read2_p0\t0'])
    assert seqs2 == set(['read3_p1\t1'])
    assert (seqsm == set(['read4_multi\t0\t1']) or
            seqsm == set(['read4_multi\t1\t0']))
    assert seqso == set(['read5_orphan'])

    seqs1 = set([r.accuracy for r in screed.open(out1)])
    seqs2 = set([r.accuracy for r in screed.open(out2)])
    seqsm = set([r.accuracy for r in screed.open(mout)])
    seqso = set([r.accuracy for r in screed.open(oout)])


def test_sweep_reads_2():

    infile = utils.get_temp_filename('seqs.fa')
    inref = utils.get_temp_filename('ref.fa')
    shutil.copyfile(utils.get_test_data('random-20-X2.fa'), infile)
    shutil.copyfile(utils.get_test_data('random-20-a.fa'), inref)
    wdir = os.path.dirname(inref)
    script = scriptpath('sweep-reads.py')
    args = ['-m', '50', '-k', '20', '-l', '9', '-b', '60', '--prefix',
            'test', '--label-by-seq', inref, infile]
    status, out, err = utils.runscript(script, args, wdir, sandbox=True)

    for i in xrange(99):
        p = os.path.join(wdir, 'test_{i}.fa'.format(i=i))
        print p, err, out
        assert os.path.exists(p)
        os.remove(p)
    assert os.path.exists(os.path.join(wdir, 'test.counts.csv'))
    assert os.path.exists(os.path.join(wdir, 'test.dist.txt'))
    assert not os.path.exists(os.path.join(wdir, 'test_multi.fa'))


def test_sweep_reads_3():

    infile = utils.get_temp_filename('seqs.fa')
    shutil.copyfile(utils.get_test_data('random-20-a.fa'), infile)
    wdir = os.path.dirname(infile)
    script = scriptpath('sweep-reads.py')
    args = ['-m', '75', '-k', '20', '-l', '1', '--prefix',
            'test', '--label-by-group', '10', infile, infile]
    status, out, err = utils.runscript(script, args, wdir, sandbox=True)

    for i in xrange(10):
        p = os.path.join(wdir, 'test_{i}.fa'.format(i=i))
        print p, err, out
        assert os.path.exists(p)
        os.remove(p)

    counts_fn = os.path.join(wdir, 'test.counts.csv')
    with open(counts_fn) as cfp:
        for line in cfp:
            _, _, c = line.partition(',')
            assert int(c) in [9, 10]

    assert os.path.exists(counts_fn)
    assert os.path.exists(os.path.join(wdir, 'test.dist.txt'))
    assert not os.path.exists(os.path.join(wdir, 'test_multi.fa'))
