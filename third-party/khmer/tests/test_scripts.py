#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#

# pylint: disable=C0111,C0103,E1103,W0612

import sys
import os
import shutil
from cStringIO import StringIO
import traceback

import khmer_tst_utils as utils
import khmer
import khmer.file
import screed


def scriptpath(script):
    return script


def teardown():
    utils.cleanup()


def test_check_space():
    # @CTB this probably belongs in a new test file, along with other
    # tests of the file.py module.
    khmer.file.check_space(['', utils.get_test_data('test-abund-read-2.fa')])


def test_load_into_counting():
    script = scriptpath('load-into-counting.py')
    args = ['-x', '1e7', '-N', '2', '-k', '20', '-t']

    outfile = utils.get_temp_filename('out.kh')
    infile = utils.get_test_data('test-abund-read-2.fa')

    args.extend([outfile, infile])

    (status, out, err) = utils.runscript(script, args)
    assert 'Total number of k-mers: 95' in err, err
    assert os.path.exists(outfile)


def test_load_into_counting_fail():
    script = scriptpath('load-into-counting.py')
    args = ['-x', '1e2', '-N', '2', '-k', '20']  # use small HT

    outfile = utils.get_temp_filename('out.kh')
    infile = utils.get_test_data('test-abund-read-2.fa')

    args.extend([outfile, infile])

    (status, out, err) = utils.runscript(script, args, fail_ok=True)
    assert status == 1, status
    assert "ERROR:" in err


def _make_counting(infilename, SIZE=1e7, N=2, K=20, BIGCOUNT=True):
    script = scriptpath('load-into-counting.py')
    args = ['-x', str(SIZE), '-N', str(N), '-k', str(K)]

    if not BIGCOUNT:
        args.append('-b')

    outfile = utils.get_temp_filename('out.kh')

    args.extend([outfile, infilename])

    utils.runscript(script, args)
    assert os.path.exists(outfile)

    return outfile


def test_filter_abund_1():
    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)
    counting_ht = _make_counting(infile, K=17)

    script = scriptpath('filter-abund.py')
    args = [counting_ht, infile]
    utils.runscript(script, args, in_dir)

    outfile = infile + '.abundfilt'
    assert os.path.exists(outfile), outfile

    seqs = set([r.sequence for r in screed.open(outfile)])
    assert len(seqs) == 1, seqs
    assert 'GGTTGACGGGGCTCAGGG' in seqs


def test_filter_abund_2():
    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)
    counting_ht = _make_counting(infile, K=17)

    script = scriptpath('filter-abund.py')
    args = ['-C', '1', counting_ht, infile, infile]
    utils.runscript(script, args, in_dir)

    outfile = infile + '.abundfilt'
    assert os.path.exists(outfile), outfile

    seqs = set([r.sequence for r in screed.open(outfile)])
    assert len(seqs) == 2, seqs
    assert 'GGTTGACGGGGCTCAGGG' in seqs

# make sure that FASTQ records are retained.


def test_filter_abund_3_fq_retained():
    infile = utils.get_temp_filename('test.fq')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fq'), infile)
    counting_ht = _make_counting(infile, K=17)

    script = scriptpath('filter-abund.py')
    args = ['-C', '1', counting_ht, infile, infile]
    utils.runscript(script, args, in_dir)

    outfile = infile + '.abundfilt'
    assert os.path.exists(outfile), outfile

    seqs = set([r.sequence for r in screed.open(outfile)])
    assert len(seqs) == 2, seqs
    assert 'GGTTGACGGGGCTCAGGG' in seqs

    # check for 'accuracy' string.
    seqs = set([r.accuracy for r in screed.open(outfile)])
    assert len(seqs) == 2, seqs
    assert '##################' in seqs


def test_filter_abund_1_singlefile():
    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)

    script = scriptpath('filter-abund-single.py')
    args = ['-x', '1e7', '-N', '2', '-k', '17', '-t', infile]
    (status, out, err) = utils.runscript(script, args, in_dir)

    assert 'Total number of k-mers: 98' in err, err

    outfile = infile + '.abundfilt'
    assert os.path.exists(outfile), outfile

    seqs = set([r.sequence for r in screed.open(outfile)])
    assert len(seqs) == 1, seqs
    assert 'GGTTGACGGGGCTCAGGG' in seqs

# test that the -V option does not trim sequences that are low abundance


def test_filter_abund_4_retain_low_abund():
    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)
    counting_ht = _make_counting(infile, K=17)

    script = scriptpath('filter-abund.py')
    args = ['-V', counting_ht, infile]
    utils.runscript(script, args, in_dir)

    outfile = infile + '.abundfilt'
    assert os.path.exists(outfile), outfile

    seqs = set([r.sequence for r in screed.open(outfile)])
    assert len(seqs) == 2, seqs
    assert 'GGTTGACGGGGCTCAGGG' in seqs

# test that the -V option *does* trim sequences that are low abundance


def test_filter_abund_5_trim_high_abund():
    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-3.fa'), infile)
    counting_ht = _make_counting(infile, K=17)

    script = scriptpath('filter-abund.py')
    args = ['-V', counting_ht, infile]
    utils.runscript(script, args, in_dir)

    outfile = infile + '.abundfilt'
    assert os.path.exists(outfile), outfile

    seqs = set([r.sequence for r in screed.open(outfile)])
    assert len(seqs) == 2, seqs

    # trimmed sequence @ error
    assert 'GGTTGACGGGGCTCAGGGGGCGGCTGACTCCGAGAGACAGC' in seqs

# test that -V/-Z setting - should not trip if -Z is set high enough.


def test_filter_abund_6_trim_high_abund_Z():
    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-3.fa'), infile)
    counting_ht = _make_counting(infile, K=17)

    script = scriptpath('filter-abund.py')
    args = ['-V', '-Z', '25', counting_ht, infile]
    utils.runscript(script, args, in_dir)

    outfile = infile + '.abundfilt'
    assert os.path.exists(outfile), outfile

    seqs = set([r.sequence for r in screed.open(outfile)])
    assert len(seqs) == 2, seqs

    # untrimmed seq.
    badseq = 'GGTTGACGGGGCTCAGGGGGCGGCTGACTCCGAGAGACAGCgtgCCGCAGCTGTCGTCAGGG' \
             'GATTTCCGGGCGG'
    assert badseq in seqs       # should be there, untrimmed


def test_filter_stoptags():
    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)
    stopfile = utils.get_temp_filename('stoptags', in_dir)

    # first, copy test-abund-read-2.fa to 'test.fa' in the temp dir.
    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)

    # now, create a file with some stop tags in it --
    K = 18
    kh = khmer.new_hashbits(K, 1, 1)
    kh.add_stop_tag('GTTGACGGGGCTCAGGGG')
    kh.save_stop_tags(stopfile)
    del kh

    # finally, run filter-stoptags.
    script = scriptpath('filter-stoptags.py')
    args = ['-k', str(K), stopfile, infile, infile]
    utils.runscript(script, args, in_dir)

    # verify that the basic output file exists
    outfile = infile + '.stopfilt'
    assert os.path.exists(outfile), outfile

    # it should contain only one unique sequence, because we've trimmed
    # off everything after the beginning of the only long sequence in there.
    seqs = set([r.sequence for r in screed.open(outfile)])
    assert len(seqs) == 1, seqs
    assert 'GGTTGACGGGGCTCAGGG' in seqs, seqs


def test_normalize_by_median_indent():
    infile = utils.get_test_data('paired-mixed.fa.pe')
    hashfile = utils.get_test_data('normC20k20.kh')
    script = scriptpath('normalize-by-median.py')
    args = ['--loadtable', hashfile, infile]
    (status, out, err) = utils.runscript(script, args)
    assert status == 0, (out, err)
    print(out, err)


def test_normalize_by_median():
    CUTOFF = '1'

    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)

    script = scriptpath('normalize-by-median.py')
    args = ['-C', CUTOFF, '-k', '17', '-t', infile]
    (status, out, err) = utils.runscript(script, args, in_dir)

    assert 'Total number of k-mers: 98' in err, err

    outfile = infile + '.keep'
    assert os.path.exists(outfile), outfile

    seqs = [r.sequence for r in screed.open(outfile)]
    assert len(seqs) == 1, seqs
    assert seqs[0].startswith('GGTTGACGGGGCTCAGGGGG'), seqs


def test_normalize_by_median_version():
    script = scriptpath('normalize-by-median.py')
    args = ['--version']
    status, out, err = utils.runscript(script, args)

    errlines = err.splitlines()
    for err in errlines:
        if err.startswith('||') or \
           not err.strip():
            continue
        break

    print errlines
    print err

    assert err.startswith('khmer ')


def test_normalize_by_median_2():
    CUTOFF = '2'

    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)

    script = scriptpath('normalize-by-median.py')
    args = ['-C', CUTOFF, '-k', '17', infile]
    utils.runscript(script, args, in_dir)

    outfile = infile + '.keep'
    assert os.path.exists(outfile), outfile

    seqs = [r.sequence for r in screed.open(outfile)]
    assert len(seqs) == 2, seqs
    assert seqs[0].startswith('GGTTGACGGGGCTCAGGGGG'), seqs
    assert seqs[1] == 'GGTTGACGGGGCTCAGGG', seqs


def test_normalize_by_median_paired():
    CUTOFF = '1'

    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-paired.fa'), infile)

    script = scriptpath('normalize-by-median.py')
    args = ['-C', CUTOFF, '-p', '-k', '17', infile]
    utils.runscript(script, args, in_dir)

    outfile = infile + '.keep'
    assert os.path.exists(outfile), outfile

    seqs = [r.sequence for r in screed.open(outfile)]
    assert len(seqs) == 2, seqs
    assert seqs[0].startswith('GGTTGACGGGGCTCAGGGGG'), seqs
    assert seqs[1].startswith('GGTTGACGGGGCTCAGGG'), seqs


def test_normalize_by_median_impaired():
    CUTOFF = '1'

    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-impaired.fa'), infile)

    script = scriptpath('normalize-by-median.py')
    args = ['-C', CUTOFF, '-p', '-k', '17', infile]
    utils.runscript(script, args, in_dir, fail_ok=True)


def test_normalize_by_median_force():
    CUTOFF = '1'

    corrupt_infile = utils.get_temp_filename('test-corrupt.fq')
    good_infile = utils.get_temp_filename('test-good.fq',
                                          tempdir=os.path.dirname(
                                              corrupt_infile))

    in_dir = os.path.dirname(good_infile)

    shutil.copyfile(utils.get_test_data('test-error-reads.fq'), corrupt_infile)
    shutil.copyfile(utils.get_test_data('test-fastq-reads.fq'), good_infile)

    script = scriptpath('normalize-by-median.py')
    args = ['-f', '-C', CUTOFF, '-k', '17', corrupt_infile, good_infile]

    (status, out, err) = utils.runscript(script, args, in_dir)

    test_ht = khmer.load_counting_hash(corrupt_infile + '.ct.failed')
    test_good_read = 'CAGGCGCCCACCACCGTGCCCTCCAACCTGATGGT'
    test_good_read2 = 'TAGTATCATCAAGGTTCAAGATGTTAATGAATAACAATTGCGCAGCAA'
    assert test_ht.count(test_good_read[:17]) > 0
    assert test_ht.count(test_good_read2[:17]) > 0
    assert os.path.exists(corrupt_infile + '.ct.failed')
    assert '*** Skipping' in err
    assert '** IOErrors' in err


def test_normalize_by_median_no_bigcount():
    infile = utils.get_temp_filename('test.fa')
    hashfile = utils.get_temp_filename('test-out.kh')
    outfile = infile + '.keep'
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)
    counting_ht = _make_counting(infile, K=8)

    script = scriptpath('normalize-by-median.py')
    args = ['-C', '1000', '-k 8', '--savetable', hashfile, infile]

    (status, out, err) = utils.runscript(script, args, in_dir)
    assert status == 0, (out, err)
    print(out, err)

    assert os.path.exists(hashfile), hashfile
    kh = khmer.load_counting_hash(hashfile)

    assert kh.get('GGTTGACG') == 255


def test_normalize_by_median_dumpfrequency():
    CUTOFF = '1'

    infiles = [utils.get_temp_filename('test-0.fq')]
    in_dir = os.path.dirname(infiles[0])
    for x in range(1, 5):
        infiles.append(utils.get_temp_filename('test-{x}.fq'.format(x=x),
                                               tempdir=in_dir))

    for infile in infiles:
        shutil.copyfile(utils.get_test_data('test-fastq-reads.fq'), infile)

    script = scriptpath('normalize-by-median.py')
    args = ['-d', '2', '-C', CUTOFF, '-k', '17']
    args.extend(infiles)

    (status, out, err) = utils.runscript(script, args, in_dir)

    test_ht = khmer.load_counting_hash(os.path.join(in_dir, 'backup.ct'))
    test_good_read = 'CAGGCGCCCACCACCGTGCCCTCCAACCTGATGGT'
    test_good_read2 = 'TAGTATCATCAAGGTTCAAGATGTTAATGAATAACAATTGCGCAGCAA'
    assert test_ht.count(test_good_read[:17]) > 0
    assert test_ht.count(test_good_read2[:17]) > 0

    assert os.path.exists(os.path.join(in_dir, 'backup.ct'))
    assert out.count('Backup: Saving') == 2
    assert 'Nothing' in out


def test_normalize_by_median_empty():
    CUTOFF = '1'

    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-empty.fa'), infile)

    script = scriptpath('normalize-by-median.py')
    args = ['-C', CUTOFF, '-k', '17', infile]
    utils.runscript(script, args, in_dir)

    outfile = infile + '.keep'
    assert os.path.exists(outfile), outfile


def test_normalize_by_median_fpr():
    MIN_TABLESIZE_PARAM = 1

    infile = utils.get_temp_filename('test-fpr.fq')
    in_dir = os.path.dirname(infile)
    shutil.copyfile(utils.get_test_data('test-fastq-reads.fq'), infile)

    script = scriptpath('normalize-by-median.py')
    args = ['-f', '-k 17', '-x ' + str(MIN_TABLESIZE_PARAM), infile]

    (status, out, err) = utils.runscript(script, args, in_dir, fail_ok=True)

    assert os.path.exists(infile + '.keep')
    assert 'fp rate estimated to be' in out, out
    assert '** ERROR: the k-mer counting table is too small' in err, err


def test_count_median():
    infile = utils.get_temp_filename('test.fa')
    outfile = infile + '.counts'

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)
    counting_ht = _make_counting(infile, K=8)

    script = scriptpath('count-median.py')
    args = [counting_ht, infile, outfile]
    utils.runscript(script, args)

    assert os.path.exists(outfile), outfile

    data = [x.strip() for x in open(outfile)]
    data = set(data)
    assert len(data) == 2, data
    assert 'seq 1001 1001.0 0.0 18' in data
    assert '895:1:37:17593:9954/1 1 103.803741455 303.702941895 114' in data

#


def test_load_graph():
    script = scriptpath('load-graph.py')
    args = ['-x', '1e7', '-N', '2', '-k', '20', '-t']

    outfile = utils.get_temp_filename('out')
    infile = utils.get_test_data('random-20-a.fa')

    args.extend([outfile, infile])

    (status, out, err) = utils.runscript(script, args)

    assert 'Total number of k-mers: 3959' in err, err

    ht_file = outfile + '.pt'
    assert os.path.exists(ht_file), ht_file

    tagset_file = outfile + '.tagset'
    assert os.path.exists(tagset_file), tagset_file

    ht = khmer.load_hashbits(ht_file)
    ht.load_tagset(tagset_file)

    # check to make sure we get the expected result for this data set
    # upon partitioning (all in one partition).  This is kind of a
    # roundabout way of checking that load-graph worked :)
    subset = ht.do_subset_partition(0, 0)
    x = ht.subset_count_partitions(subset)
    assert x == (1, 0), x


def test_load_graph_no_tags():
    script = scriptpath('load-graph.py')
    args = ['-x', '1e7', '-N', '2', '-k', '20', '-n']

    outfile = utils.get_temp_filename('out')
    infile = utils.get_test_data('random-20-a.fa')

    args.extend([outfile, infile])

    utils.runscript(script, args)

    ht_file = outfile + '.pt'
    assert os.path.exists(ht_file), ht_file

    tagset_file = outfile + '.tagset'
    assert not os.path.exists(tagset_file), tagset_file

    assert khmer.load_hashbits(ht_file)

    # can't think of a good way to make sure this worked, beyond just
    # loading the ht file...


def test_load_graph_fail():
    script = scriptpath('load-graph.py')
    args = ['-x', '1e3', '-N', '2', '-k', '20']  # use small HT

    outfile = utils.get_temp_filename('out')
    infile = utils.get_test_data('random-20-a.fa')

    args.extend([outfile, infile])

    (status, out, err) = utils.runscript(script, args, fail_ok=True)
    assert status == 1, status
    assert "ERROR:" in err


def test_load_graph_write_fp():
    script = scriptpath('load-graph.py')
    args = ['-x', '1e5', '-N', '2', '-k', '20', '-w']  # use small HT

    outfile = utils.get_temp_filename('out')
    infile = utils.get_test_data('random-20-a.fa')

    args.extend([outfile, infile])

    (status, out, err) = utils.runscript(script, args)

    ht_file = outfile + '.pt'
    assert os.path.exists(ht_file), ht_file

    info_file = outfile + '.info'
    assert os.path.exists(info_file), info_file
    data = [x.strip() for x in open(info_file)]
    data = set(data)
    assert '3959 unique k-mers' in data
    assert 'false positive rate estimated to be 0.002' in data


def test_load_graph_multithread():
    script = scriptpath('load-graph.py')

    outfile = utils.get_temp_filename('test')
    infile = utils.get_test_data('test-reads.fa')

    args = ['-N', '4', '-x', '1e9', '-T', '8', outfile, infile]

    (status, out, err) = utils.runscript(script, args)


def _make_graph(infilename, min_hashsize=1e7, n_hashes=2, ksize=20,
                do_partition=False,
                annotate_partitions=False,
                stop_big_traverse=False):
    script = scriptpath('load-graph.py')
    args = ['-x', str(min_hashsize), '-N', str(n_hashes), '-k', str(ksize)]

    outfile = utils.get_temp_filename('out')
    infile = infilename

    args.extend([outfile, infile])

    utils.runscript(script, args)

    ht_file = outfile + '.pt'
    assert os.path.exists(ht_file), ht_file

    tagset_file = outfile + '.tagset'
    assert os.path.exists(tagset_file), tagset_file

    if do_partition:
        script = scriptpath('partition-graph.py')
        args = [outfile]
        if stop_big_traverse:
            args.insert(0, '--no-big-traverse')
        utils.runscript(script, args)

        script = scriptpath('merge-partitions.py')
        args = [outfile, '-k', str(ksize)]
        utils.runscript(script, args)

        final_pmap_file = outfile + '.pmap.merged'
        assert os.path.exists(final_pmap_file)

        if annotate_partitions:
            script = scriptpath('annotate-partitions.py')
            args = ["-k", str(ksize), outfile, infilename]

            in_dir = os.path.dirname(outfile)
            utils.runscript(script, args, in_dir)

            baseinfile = os.path.basename(infilename)
            assert os.path.exists(os.path.join(in_dir, baseinfile + '.part'))

    return outfile


def _DEBUG_make_graph(infilename, min_hashsize=1e7, n_hashes=2, ksize=20,
                      do_partition=False,
                      annotate_partitions=False,
                      stop_big_traverse=False):
    script = scriptpath('load-graph.py')
    args = ['-x', str(min_hashsize), '-N', str(n_hashes), '-k', str(ksize)]

    outfile = utils.get_temp_filename('out')
    infile = utils.get_test_data(infilename)

    args.extend([outfile, infile])

    utils.runscript(script, args)

    ht_file = outfile + '.ct'
    assert os.path.exists(ht_file), ht_file

    tagset_file = outfile + '.tagset'
    assert os.path.exists(tagset_file), tagset_file

    if do_partition:
        print ">>>> DEBUG: Partitioning <<<"
        script = scriptpath('partition-graph.py')
        args = [outfile]
        if stop_big_traverse:
            args.insert(0, '--no-big-traverse')
        utils.runscript(script, args)

        print ">>>> DEBUG: Merging Partitions <<<"
        script = scriptpath('merge-partitions.py')
        args = [outfile, '-k', str(ksize)]
        utils.runscript(script, args)

        final_pmap_file = outfile + '.pmap.merged'
        assert os.path.exists(final_pmap_file)

        if annotate_partitions:
            print ">>>> DEBUG: Annotating Partitions <<<"
            script = scriptpath('annotate-partitions.py')
            args = ["-k", str(ksize), outfile, infilename]

            in_dir = os.path.dirname(outfile)
            utils.runscript(script, args, in_dir)

            baseinfile = os.path.basename(infilename)
            assert os.path.exists(os.path.join(in_dir, baseinfile + '.part'))

    return outfile


def test_partition_graph_1():
    graphbase = _make_graph(utils.get_test_data('random-20-a.fa'))

    script = scriptpath('partition-graph.py')
    args = [graphbase]

    utils.runscript(script, args)

    script = scriptpath('merge-partitions.py')
    args = [graphbase, '-k', str(20)]
    utils.runscript(script, args)

    final_pmap_file = graphbase + '.pmap.merged'
    assert os.path.exists(final_pmap_file)

    ht = khmer.load_hashbits(graphbase + '.pt')
    ht.load_tagset(graphbase + '.tagset')
    ht.load_partitionmap(final_pmap_file)

    x = ht.count_partitions()
    assert x == (1, 0), x          # should be exactly one partition.


def test_partition_graph_nojoin_k21():
    # test with K=21
    graphbase = _make_graph(utils.get_test_data('random-20-a.fa'), ksize=21)

    script = scriptpath('partition-graph.py')
    args = [graphbase]

    utils.runscript(script, args)

    script = scriptpath('merge-partitions.py')
    args = [graphbase, '-k', str(21)]
    utils.runscript(script, args)

    final_pmap_file = graphbase + '.pmap.merged'
    assert os.path.exists(final_pmap_file)

    ht = khmer.load_hashbits(graphbase + '.pt')
    ht.load_tagset(graphbase + '.tagset')
    ht.load_partitionmap(final_pmap_file)

    x = ht.count_partitions()
    assert x == (99, 0), x          # should be 99 partitions at K=21


def test_partition_graph_nojoin_stoptags():
    # test with stoptags
    graphbase = _make_graph(utils.get_test_data('random-20-a.fa'))

    # add in some stop tags
    ht = khmer.load_hashbits(graphbase + '.pt')
    ht.add_stop_tag('TTGCATACGTTGAGCCAGCG')
    stoptags_file = graphbase + '.stoptags'
    ht.save_stop_tags(stoptags_file)
    del ht

    # run script with stoptags option
    script = scriptpath('partition-graph.py')
    args = ['--stoptags', stoptags_file, graphbase]

    utils.runscript(script, args)

    script = scriptpath('merge-partitions.py')
    args = [graphbase, '-k', str(20)]
    utils.runscript(script, args)

    final_pmap_file = graphbase + '.pmap.merged'
    assert os.path.exists(final_pmap_file)

    ht = khmer.load_hashbits(graphbase + '.pt')
    ht.load_tagset(graphbase + '.tagset')
    ht.load_partitionmap(final_pmap_file)

    x = ht.count_partitions()
    assert x == (2, 0), x          # should be 2 partitions


def test_partition_graph_big_traverse():
    graphbase = _make_graph(utils.get_test_data('biglump-random-20-a.fa'),
                            do_partition=True, stop_big_traverse=False)

    final_pmap_file = graphbase + '.pmap.merged'
    assert os.path.exists(final_pmap_file)

    ht = khmer.load_hashbits(graphbase + '.pt')
    ht.load_tagset(graphbase + '.tagset')
    ht.load_partitionmap(final_pmap_file)

    x = ht.count_partitions()
    assert x == (1, 0), x          # should be exactly one partition.


def test_partition_graph_no_big_traverse():
    # do NOT exhaustively traverse
    graphbase = _make_graph(utils.get_test_data('biglump-random-20-a.fa'),
                            do_partition=True, stop_big_traverse=True)

    final_pmap_file = graphbase + '.pmap.merged'
    assert os.path.exists(final_pmap_file)

    ht = khmer.load_hashbits(graphbase + '.pt')
    ht.load_tagset(graphbase + '.tagset')
    ht.load_partitionmap(final_pmap_file)

    x = ht.count_partitions()
    assert x[0] == 4, x       # should be four partitions, broken at knot.


def test_annotate_partitions():
    seqfile = utils.get_test_data('random-20-a.fa')
    graphbase = _make_graph(seqfile, do_partition=True)
    in_dir = os.path.dirname(graphbase)

    # get the final pmap file
    final_pmap_file = graphbase + '.pmap.merged'
    assert os.path.exists(final_pmap_file)

    script = scriptpath('annotate-partitions.py')
    args = ["-k", "20", graphbase, seqfile]
    utils.runscript(script, args, in_dir)

    partfile = os.path.join(in_dir, 'random-20-a.fa.part')

    parts = [r.name.split('\t')[1] for r in screed.open(partfile)]
    parts = set(parts)
    assert '2' in parts
    assert len(parts) == 1


def test_annotate_partitions_2():
    # test with K=21 (no joining of sequences)
    seqfile = utils.get_test_data('random-20-a.fa')
    graphbase = _make_graph(seqfile, do_partition=True,
                            ksize=21)
    in_dir = os.path.dirname(graphbase)

    # get the final pmap file
    final_pmap_file = graphbase + '.pmap.merged'
    assert os.path.exists(final_pmap_file)

    script = scriptpath('annotate-partitions.py')
    args = ["-k", "21", graphbase, seqfile]
    utils.runscript(script, args, in_dir)

    partfile = os.path.join(in_dir, 'random-20-a.fa.part')

    parts = [r.name.split('\t')[1] for r in screed.open(partfile)]
    parts = set(parts)
    print parts
    assert len(parts) == 99, len(parts)


def test_extract_partitions():
    seqfile = utils.get_test_data('random-20-a.fa')
    graphbase = _make_graph(
        seqfile, do_partition=True, annotate_partitions=True)
    in_dir = os.path.dirname(graphbase)

    # get the final part file
    partfile = os.path.join(in_dir, 'random-20-a.fa.part')

    # ok, now run extract-partitions.
    script = scriptpath('extract-partitions.py')
    args = ['extracted', partfile]

    utils.runscript(script, args, in_dir)

    distfile = os.path.join(in_dir, 'extracted.dist')
    groupfile = os.path.join(in_dir, 'extracted.group0000.fa')
    assert os.path.exists(distfile)
    assert os.path.exists(groupfile)

    dist = open(distfile).readline()
    assert dist.strip() == '99 1 1 99'

    parts = [r.name.split('\t')[1] for r in screed.open(partfile)]
    assert len(parts) == 99, len(parts)
    parts = set(parts)
    assert len(parts) == 1, len(parts)


def test_extract_partitions_fq():
    seqfile = utils.get_test_data('random-20-a.fq')
    graphbase = _make_graph(
        seqfile, do_partition=True, annotate_partitions=True)
    in_dir = os.path.dirname(graphbase)

    # get the final part file
    partfile = os.path.join(in_dir, 'random-20-a.fq.part')

    # ok, now run extract-partitions.
    script = scriptpath('extract-partitions.py')
    args = ['extracted', partfile]

    utils.runscript(script, args, in_dir)

    distfile = os.path.join(in_dir, 'extracted.dist')
    groupfile = os.path.join(in_dir, 'extracted.group0000.fq')
    assert os.path.exists(distfile)
    assert os.path.exists(groupfile)

    dist = open(distfile).readline()
    assert dist.strip() == '99 1 1 99'

    parts = [r.name.split('\t')[1] for r in screed.open(partfile)]
    assert len(parts) == 99, len(parts)
    parts = set(parts)
    assert len(parts) == 1, len(parts)

    quals = set([r.accuracy for r in screed.open(partfile)])
    quals = list(quals)
    assert quals[0], quals


def test_extract_partitions_output_unassigned():
    seqfile = utils.get_test_data('random-20-a.fa')
    graphbase = _make_graph(
        seqfile, do_partition=True, annotate_partitions=True)
    in_dir = os.path.dirname(graphbase)

    # get the final part file
    partfile = os.path.join(in_dir, 'random-20-a.fa.part')

    # ok, now run extract-partitions.
    script = scriptpath('extract-partitions.py')
    args = ['-U', 'extracted', partfile]

    utils.runscript(script, args, in_dir)

    distfile = os.path.join(in_dir, 'extracted.dist')
    groupfile = os.path.join(in_dir, 'extracted.group0000.fa')
    unassigned_file = os.path.join(in_dir, 'extracted.unassigned.fa')
    assert os.path.exists(distfile)
    assert os.path.exists(groupfile)
    assert os.path.exists(unassigned_file)

    dist = open(distfile).readline()
    assert dist.strip() == '99 1 1 99'

    parts = [r.name.split('\t')[1] for r in screed.open(partfile)]
    assert len(parts) == 99, len(parts)
    parts = set(parts)
    assert len(parts) == 1, len(parts)


def test_extract_partitions_no_output_groups():
    seqfile = utils.get_test_data('random-20-a.fq')
    graphbase = _make_graph(
        seqfile, do_partition=True, annotate_partitions=True)
    in_dir = os.path.dirname(graphbase)

    # get the final part file
    partfile = os.path.join(in_dir, 'random-20-a.fa.part')

    # ok, now run extract-partitions.
    script = scriptpath('extract-partitions.py')
    args = ['-n', 'extracted', partfile]

    # We expect a sys.exit -> we need the test to be tolerant
    utils.runscript(script, args, in_dir, fail_ok=True)

    # Group files are created after output_groups is
    # checked. They should not exist in this scenario
    groupfile = os.path.join(in_dir, 'extracted.group0000.fa')
    assert not os.path.exists(groupfile)


def test_extract_partitions_pid_0():
    basefile = utils.get_test_data('random-20-a.fa.part')
    partfile = utils.get_temp_filename('random-20-a.fa.part')
    shutil.copyfile(basefile, partfile)

    in_dir = os.path.dirname(partfile)
    # ok, now run extract-partitions.
    script = scriptpath('extract-partitions.py')
    args = ['-U', 'extracted', partfile]

    utils.runscript(script, args, in_dir)

    distfile = os.path.join(in_dir, 'extracted.dist')
    groupfile = os.path.join(in_dir, 'extracted.group0000.fa')
    unassigned_file = os.path.join(in_dir, 'extracted.unassigned.fa')
    assert os.path.exists(distfile)
    assert os.path.exists(groupfile)
    assert os.path.exists(unassigned_file)

    # Assert unassigned file not empty
    unassigned_content = open(unassigned_file).readline()
    assert unassigned_content.strip().split('\t')[0] != ''


def test_extract_partitions_multi_groups():
    basefile = utils.get_test_data('random-20-a.fa.part')
    partfile = utils.get_temp_filename('random-20-a.fa.part')
    shutil.copyfile(basefile, partfile)

    in_dir = os.path.dirname(partfile)

    # ok, now run extract-partitions.
    script = scriptpath('extract-partitions.py')
    args = ['-m', '1', '-X', '1', 'extracted', partfile]

    utils.runscript(script, args, in_dir)

    # Multiple group files are created after should be created
    groupfile1 = os.path.join(in_dir, 'extracted.group0000.fa')
    groupfile2 = os.path.join(in_dir, 'extracted.group0001.fa')
    groupfile3 = os.path.join(in_dir, 'extracted.group0002.fa')
    assert os.path.exists(groupfile1)
    assert os.path.exists(groupfile2)
    assert os.path.exists(groupfile3)


def test_extract_partitions_no_groups():
    empty_file = utils.get_temp_filename('empty-file')
    basefile = utils.get_test_data('empty-file')

    shutil.copyfile(basefile, empty_file)
    in_dir = os.path.dirname(empty_file)

    # ok, now run extract-partitions.
    script = scriptpath('extract-partitions.py')
    args = ['extracted', empty_file]

    utils.runscript(script, args, in_dir, fail_ok=True)

    # No group files should be created
    groupfile = os.path.join(in_dir, 'extracted.group0000.fa')

    assert not os.path.exists(groupfile)


def test_abundance_dist():
    infile = utils.get_temp_filename('test.fa')
    outfile = utils.get_temp_filename('test.dist')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)

    htfile = _make_counting(infile, K=17)

    script = scriptpath('abundance-dist.py')
    args = ['-z', htfile, infile, outfile]
    utils.runscript(script, args, in_dir)

    fp = iter(open(outfile))
    line = fp.next().strip()
    assert line == '1 96 96 0.98', line
    line = fp.next().strip()
    assert line == '1001 2 98 1.0', line


def test_abundance_dist_nobigcount():
    infile = utils.get_temp_filename('test.fa')
    outfile = utils.get_temp_filename('test.dist')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)

    htfile = _make_counting(infile, K=17, BIGCOUNT=False)

    script = scriptpath('abundance-dist.py')
    args = ['-z', htfile, infile, outfile]
    utils.runscript(script, args, in_dir)

    fp = iter(open(outfile))
    line = fp.next().strip()
    assert line == '1 96 96 0.98', line
    line = fp.next().strip()
    assert line == '255 2 98 1.0', line


def test_abundance_dist_single():
    infile = utils.get_temp_filename('test.fa')
    outfile = utils.get_temp_filename('test.dist')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)

    script = scriptpath('abundance-dist-single.py')
    args = ['-x', '1e7', '-N', '2', '-k', '17', '-z', '-t', infile,
            outfile]
    (status, out, err) = utils.runscript(script, args, in_dir)

    assert 'Total number of k-mers: 98' in err, err

    fp = iter(open(outfile))
    line = fp.next().strip()
    assert line == '1 96 96 0.98', line
    line = fp.next().strip()
    assert line == '1001 2 98 1.0', line


def test_abundance_dist_single_nobigcount():
    infile = utils.get_temp_filename('test.fa')
    outfile = utils.get_temp_filename('test.dist')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)

    script = scriptpath('abundance-dist-single.py')
    args = ['-x', '1e7', '-N', '2', '-k', '17', '-z', '-b', infile, outfile]
    utils.runscript(script, args, in_dir)

    fp = iter(open(outfile))
    line = fp.next().strip()
    assert line == '1 96 96 0.98', line
    line = fp.next().strip()
    assert line == '255 2 98 1.0', line


def test_do_partition():
    seqfile = utils.get_test_data('random-20-a.fa')
    graphbase = utils.get_temp_filename('out')
    in_dir = os.path.dirname(graphbase)

    script = scriptpath('do-partition.py')
    args = ["-k", "20", graphbase, seqfile]

    utils.runscript(script, args, in_dir)

    partfile = os.path.join(in_dir, 'random-20-a.fa.part')

    parts = [r.name.split('\t')[1] for r in screed.open(partfile)]
    parts = set(parts)
    assert '2' in parts
    assert len(parts) == 1


def test_do_partition_2():
    # test with K=21 (no joining of sequences)
    seqfile = utils.get_test_data('random-20-a.fa')
    graphbase = utils.get_temp_filename('out')
    in_dir = os.path.dirname(graphbase)

    script = scriptpath('do-partition.py')
    args = ["-k", "21", graphbase, seqfile]

    utils.runscript(script, args, in_dir)

    partfile = os.path.join(in_dir, 'random-20-a.fa.part')

    parts = [r.name.split('\t')[1] for r in screed.open(partfile)]
    parts = set(parts)

    assert len(parts) == 99, len(parts)

#


def test_interleave_reads_1_fq():
    # test input files
    infile1 = utils.get_test_data('paired.fq.1')
    infile2 = utils.get_test_data('paired.fq.2')

    # correct output
    ex_outfile = utils.get_test_data('paired.fq')

    # actual output file
    outfile = utils.get_temp_filename('out.fq')

    script = scriptpath('interleave-reads.py')
    args = [infile1, infile2, '-o', outfile]

    utils.runscript(script, args)

    r = open(ex_outfile).read()
    q = open(outfile).read()

    assert r == q, (r, q)


def test_interleave_reads_2_fa():
    # test input files
    infile1 = utils.get_test_data('paired.fa.1')
    infile2 = utils.get_test_data('paired.fa.2')

    # correct output
    ex_outfile = utils.get_test_data('paired.fa')

    # actual output file
    outfile = utils.get_temp_filename('out.fa')

    script = scriptpath('interleave-reads.py')
    args = [infile1, infile2, '-o', outfile]

    utils.runscript(script, args)

    n = 0
    for r, q in zip(screed.open(ex_outfile), screed.open(outfile)):
        n += 1
        assert r.name == q.name
        assert r.sequence == q.sequence
    assert n > 0


def test_make_initial_stoptags():
    # gen input files using load-graph.py -t
    # should keep test_data directory size down
    # or something like that
    # this assumes (obv.) load-graph works properly
    bzinfile = utils.get_temp_filename('test-reads.fq.bz2')
    shutil.copyfile(utils.get_test_data('test-reads.fq.bz2'), bzinfile)
    in_dir = os.path.dirname(bzinfile)

    genscript = scriptpath('load-graph.py')
    genscriptargs = ['-t', 'test-reads', 'test-reads.fq.bz2']
    utils.runscript(genscript, genscriptargs, in_dir)

    # test input file gen'd by load-graphs
    infile = utils.get_temp_filename('test-reads.pt')
    infile2 = utils.get_temp_filename('test-reads.tagset', in_dir)

    # get file to compare against
    ex_outfile = utils.get_test_data('test-reads.stoptags')

    # actual output file
    outfile1 = utils.get_temp_filename('test-reads.stoptags', in_dir)

    script = scriptpath('make-initial-stoptags.py')
    # make-initial-stoptags has weird file argument syntax
    # read the code before modifying
    args = ['test-reads']

    utils.runscript(script, args, in_dir)
    assert os.path.exists(outfile1), outfile1


def test_extract_paired_reads_1_fa():
    # test input file
    infile = utils.get_test_data('paired-mixed.fa')

    ex_outfile1 = utils.get_test_data('paired-mixed.fa.pe')
    ex_outfile2 = utils.get_test_data('paired-mixed.fa.se')

    # actual output files...
    outfile1 = utils.get_temp_filename('paired-mixed.fa.pe')
    in_dir = os.path.dirname(outfile1)
    outfile2 = utils.get_temp_filename('paired-mixed.fa.se', in_dir)

    script = scriptpath('extract-paired-reads.py')
    args = [infile]

    utils.runscript(script, args, in_dir)

    assert os.path.exists(outfile1), outfile1
    assert os.path.exists(outfile2), outfile2

    n = 0
    for r, q in zip(screed.open(ex_outfile1), screed.open(outfile1)):
        n += 1
        assert r.name == q.name
        assert r.sequence == q.sequence
    assert n > 0

    n = 0
    for r, q in zip(screed.open(ex_outfile2), screed.open(outfile2)):
        n += 1
        assert r.name == q.name
        assert r.sequence == q.sequence
    assert n > 0


def test_extract_paired_reads_2_fq():
    # test input file
    infile = utils.get_test_data('paired-mixed.fq')

    ex_outfile1 = utils.get_test_data('paired-mixed.fq.pe')
    ex_outfile2 = utils.get_test_data('paired-mixed.fq.se')

    # actual output files...
    outfile1 = utils.get_temp_filename('paired-mixed.fq.pe')
    in_dir = os.path.dirname(outfile1)
    outfile2 = utils.get_temp_filename('paired-mixed.fq.se', in_dir)

    script = scriptpath('extract-paired-reads.py')
    args = [infile]

    utils.runscript(script, args, in_dir)

    assert os.path.exists(outfile1), outfile1
    assert os.path.exists(outfile2), outfile2

    n = 0
    for r, q in zip(screed.open(ex_outfile1), screed.open(outfile1)):
        n += 1
        assert r.name == q.name
        assert r.sequence == q.sequence
        assert r.accuracy == q.accuracy
    assert n > 0

    n = 0
    for r, q in zip(screed.open(ex_outfile2), screed.open(outfile2)):
        n += 1
        assert r.name == q.name
        assert r.sequence == q.sequence
        assert r.accuracy == q.accuracy
    assert n > 0


def test_split_paired_reads_1_fa():
    # test input file
    infile = utils.get_test_data('paired.fa')

    ex_outfile1 = utils.get_test_data('paired.fa.1')
    ex_outfile2 = utils.get_test_data('paired.fa.2')

    # actual output files...
    outfile1 = utils.get_temp_filename('paired.fa.1')
    in_dir = os.path.dirname(outfile1)
    outfile2 = utils.get_temp_filename('paired.fa.2', in_dir)

    script = scriptpath('split-paired-reads.py')
    args = [infile]

    utils.runscript(script, args, in_dir)

    assert os.path.exists(outfile1), outfile1
    assert os.path.exists(outfile2), outfile2

    n = 0
    for r, q in zip(screed.open(ex_outfile1), screed.open(outfile1)):
        n += 1
        assert r.name == q.name
        assert r.sequence == q.sequence
    assert n > 0

    n = 0
    for r, q in zip(screed.open(ex_outfile2), screed.open(outfile2)):
        n += 1
        assert r.name == q.name
        assert r.sequence == q.sequence
    assert n > 0


def test_split_paired_reads_2_fq():
    # test input file
    infile = utils.get_test_data('paired.fq')

    ex_outfile1 = utils.get_test_data('paired.fq.1')
    ex_outfile2 = utils.get_test_data('paired.fq.2')

    # actual output files...
    outfile1 = utils.get_temp_filename('paired.fq.1')
    in_dir = os.path.dirname(outfile1)
    outfile2 = utils.get_temp_filename('paired.fq.2', in_dir)

    script = scriptpath('split-paired-reads.py')
    args = [infile]

    utils.runscript(script, args, in_dir)

    assert os.path.exists(outfile1), outfile1
    assert os.path.exists(outfile2), outfile2

    n = 0
    for r, q in zip(screed.open(ex_outfile1), screed.open(outfile1)):
        n += 1
        assert r.name == q.name
        assert r.sequence == q.sequence
        assert r.accuracy == q.accuracy
    assert n > 0

    n = 0
    for r, q in zip(screed.open(ex_outfile2), screed.open(outfile2)):
        n += 1
        assert r.name == q.name
        assert r.sequence == q.sequence
        assert r.accuracy == q.accuracy
    assert n > 0


def test_sample_reads_randomly():
    infile = utils.get_temp_filename('test.fq')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-fastq-reads.fq'), infile)

    script = scriptpath('sample-reads-randomly.py')
    # fix random number seed for reproducibility
    args = ['-N', '10', '-R', '1']
    args.append(infile)
    utils.runscript(script, args, in_dir)

    outfile = infile + '.subset'
    assert os.path.exists(outfile), outfile

    seqs = set([r.name for r in screed.open(outfile)])
    assert seqs == set(['895:1:1:1326:7273', '895:1:1:1373:4848',
                        '895:1:1:1264:15854', '895:1:1:1338:15407',
                        '895:1:1:1327:15301', '895:1:1:1265:2265',
                        '895:1:1:1327:13028', '895:1:1:1368:4434',
                        '895:1:1:1335:19932', '895:1:1:1340:19387'])


def test_fastq_to_fasta():

    script = scriptpath('fastq-to-fasta.py')
    clean_infile = utils.get_temp_filename('test-clean.fq')
    n_infile = utils.get_temp_filename('test-n.fq')

    shutil.copyfile(utils.get_test_data('test-fastq-reads.fq'), clean_infile)
    shutil.copyfile(utils.get_test_data('test-fastq-n-reads.fq'), n_infile)

    clean_outfile = clean_infile + '.keep.fa'
    n_outfile = n_infile + '.keep.fa'

    in_dir = os.path.dirname(clean_infile)
    in_dir_n = os.path.dirname(n_infile)

    args = [clean_infile, '-n', '-o', clean_outfile]
    (status, out, err) = utils.runscript(script, args, in_dir)
    assert len(out.splitlines()) == 2, len(out.splitlines())
    assert "No lines dropped" in err

    args = [n_infile, '-n', '-o', n_outfile]
    (status, out, err) = utils.runscript(script, args, in_dir_n)
    assert len(out.splitlines()) == 2
    assert "No lines dropped" in err

    args = [clean_infile, '-o', clean_outfile]
    (status, out, err) = utils.runscript(script, args, in_dir)
    assert len(out.splitlines()) == 2
    assert "0 lines dropped" in err

    args = [n_infile, '-o', n_outfile]
    (status, out, err) = utils.runscript(script, args, in_dir_n)
    assert len(out.splitlines()) == 2, out
    assert "4 lines dropped" in err, err

    args = [clean_infile]
    (status, out, err) = utils.runscript(script, args, in_dir)
    assert len(out.splitlines()) > 2
    assert "0 lines dropped" in err

    args = [n_infile]
    (status, out, err) = utils.runscript(script, args, in_dir_n)
    assert len(out.splitlines()) > 2
    assert "4 lines dropped" in err


def test_extract_long_sequences():

    script = scriptpath('extract-long-sequences.py')
    fq_infile = utils.get_temp_filename('test.fq')
    fa_infile = utils.get_temp_filename('test.fa')

    shutil.copyfile(utils.get_test_data('paired-mixed.fq'), fq_infile)
    shutil.copyfile(utils.get_test_data('paired-mixed.fa'), fa_infile)

    fq_outfile = fq_infile + '.keep.fq'
    fa_outfile = fa_infile + '.keep.fa'

    in_dir_fq = os.path.dirname(fq_infile)
    in_dir_fa = os.path.dirname(fa_infile)

    args = [fq_infile, '-l', '10', '-o', 'fq_outfile']
    (status, out, err) = utils.runscript(script, args, in_dir_fa)

    countlines = sum(1 for line in open(fq_infile))
    assert countlines == 44, countlines

    args = [fa_infile, '-l', '10', '-o', 'fa_outfile']
    (status, out, err) = utils.runscript(script, args, in_dir_fa)

    countlines = sum(1 for line in open(fa_infile))
    assert countlines == 22, countlines


def test_sample_reads_randomly_S():
    infile = utils.get_temp_filename('test.fq')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-fastq-reads.fq'), infile)

    script = scriptpath('sample-reads-randomly.py')

    # fix random number seed for reproducibility
    args = ['-N', '10', '-R', '1', '-S', '3']

    badargs = list(args)
    badargs.extend(['-o', 'test', 'test.fq', 'test.fq'])
    (status, out, err) = utils.runscript(script, badargs, in_dir, fail_ok=True)
    assert status == -1, (status, out, err)

    args.append('test.fq')

    utils.runscript(script, args, in_dir)

    outfile = infile + '.subset.0'
    assert os.path.exists(outfile), outfile

    seqs = set([r.name for r in screed.open(outfile)])
    print seqs
    assert seqs == set(['895:1:1:1298:13380', '895:1:1:1347:3237',
                        '895:1:1:1295:6189', '895:1:1:1342:11001',
                        '895:1:1:1252:19493', '895:1:1:1318:10532',
                        '895:1:1:1314:10430', '895:1:1:1347:8723',
                        '895:1:1:1381:4958', '895:1:1:1338:6614'])

    outfile = infile + '.subset.1'
    assert os.path.exists(outfile), outfile

    seqs = set([r.name for r in screed.open(outfile)])
    print seqs
    assert seqs == set(['895:1:1:1384:20217', '895:1:1:1347:3237',
                        '895:1:1:1348:18672', '895:1:1:1290:11501',
                        '895:1:1:1386:7536', '895:1:1:1373:13994',
                        '895:1:1:1355:13535', '895:1:1:1303:6251',
                        '895:1:1:1381:4958', '895:1:1:1338:6614'])

    outfile = infile + '.subset.2'
    assert os.path.exists(outfile), outfile

    seqs = set([r.name for r in screed.open(outfile)])
    print seqs
    assert seqs == set(['895:1:1:1326:7273', '895:1:1:1384:20217',
                        '895:1:1:1347:3237', '895:1:1:1353:6642',
                        '895:1:1:1340:19387', '895:1:1:1252:19493',
                        '895:1:1:1381:7062', '895:1:1:1383:3089',
                        '895:1:1:1342:20695', '895:1:1:1303:6251'])


def test_count_overlap():
    seqfile1 = utils.get_temp_filename('test-overlap1.fa')
    in_dir = os.path.dirname(seqfile1)
    seqfile2 = utils.get_temp_filename('test-overlap2.fa', in_dir)
    outfile = utils.get_temp_filename('overlap.out', in_dir)
    curvefile = utils.get_temp_filename('overlap.out.curve', in_dir)
    shutil.copy(utils.get_test_data('test-overlap1.fa'), seqfile1)
    shutil.copy(utils.get_test_data('test-overlap2.fa'), seqfile2)
    htfile = _make_graph(seqfile1, ksize=20)
    script = scriptpath('count-overlap.py')
    args = ['--ksize', '20', '--n_tables', '2', '--min-tablesize', '10000000',
            htfile + '.pt', seqfile2, outfile]
    (status, out, err) = utils.runscript(script, args, in_dir)
    assert status == 0
    assert os.path.exists(outfile), outfile
    data = [x.strip() for x in open(outfile)]
    data = set(data)
    assert '# of unique k-mers in dataset2: 759047' in data
    assert '# of overlap unique k-mers: 245621' in data
    assert os.path.exists(curvefile), curvefile
    data = [x.strip() for x in open(curvefile)]
    data = set(data)
    assert '178633 1155' in data
    assert '496285 2970' in data
    assert '752053 238627' in data
