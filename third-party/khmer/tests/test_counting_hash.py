#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
# pylint: disable=missing-docstring,protected-access
import gzip

import khmer
import khmer_tst_utils as utils
from khmer import ReadParser

from nose.plugins.attrib import attr

MAX_COUNT = 255
MAX_BIGCOUNT = 65535

#

# from http://www.rsok.com/~jrm/printprimes.html
PRIMES_1m = [1000003, 1009837]
PRIMES_100m = [100009979, 100000007]
PRIMES_1b = [1000000007, 1000000919]
PRIMES_2b = [1999999973, 1999999943]
PRIMES_4b = [4000000007, 4000000009]
PRIMES_8b = [8000000011, 8000000051]

DNA = "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"


def teardown():
    utils.cleanup()


class Test_CountingHash(object):

    def setup(self):
        self.hi = khmer._new_counting_hash(12, PRIMES_1m)

    def test_collision_1(self):

        GG = 'G' * 12                   # forward_hash: 11184810
        assert khmer.forward_hash(GG, 12) == 11184810

        collision_1 = 'AAACGTATGACT'
        assert khmer.forward_hash(collision_1, 12) == 184777L

        collision_2 = 'AAATACCGAGCG'
        assert khmer.forward_hash(collision_2, 12) == 76603L

        # note, hash(GG) % 1000003 == hash(collision_1)
        # note, hash(GG) % 1009837 == hash(collision_2)

        hi = self.hi
        hi.consume(GG)
        hi.consume(collision_1)

        assert hi.get(GG) == 1

    def test_collision_2(self):

        GG = 'G' * 12                   # forward_hash: 11184810
        assert khmer.forward_hash(GG, 12) == 11184810

        collision_1 = 'AAACGTATGACT'
        assert khmer.forward_hash(collision_1, 12) == 184777L

        collision_2 = 'AAATACCGAGCG'
        assert khmer.forward_hash(collision_2, 12) == 76603L

        # hash(GG) % 1000003 == hash(collision_1)
        # hash(GG) % 1009837 == hash(collision_2)

        hi = self.hi
        hi.consume(GG)
        hi.consume(collision_2)

        assert hi.get(GG) == 1

    def test_collision_3(self):

        GG = 'G' * 12                   # forward_hash: 11184810
        assert khmer.forward_hash(GG, 12) == 11184810

        collision_1 = 'AAACGTATGACT'
        assert khmer.forward_hash(collision_1, 12) == 184777L

        collision_2 = 'AAATACCGAGCG'
        assert khmer.forward_hash(collision_2, 12) == 76603L

        # hash(GG) % 1000003 == hash(collision_1)
        # hash(GG) % 1009837 == hash(collision_2)

        hi = self.hi
        hi.consume(GG)
        hi.consume(collision_1)
        hi.consume(collision_2)

        assert hi.get(GG) == 2


def test_3_tables():
    x = list(PRIMES_1m)
    x.append(1000005)

    hi = khmer._new_counting_hash(12, x)

    GG = 'G' * 12                   # forward_hash: 11184810
    assert khmer.forward_hash(GG, 12) == 11184810

    collision_1 = 'AAACGTATGACT'
    assert khmer.forward_hash(collision_1, 12) == 184777L

    collision_2 = 'AAATACCGAGCG'
    assert khmer.forward_hash(collision_2, 12) == 76603L

    collision_3 = 'AAACGTATCGAG'
    assert khmer.forward_hash(collision_3, 12) == 184755L

    # hash(GG) % 1000003 == hash(collision_1)
    # hash(GG) % 1009837 == hash(collision_2)
    # hash(GG) % 1000005 == hash(collision_3)
    hi.consume(GG)
    assert hi.get(GG) == 1

    hi.consume(collision_1)
    assert hi.get(GG) == 1

    hi.consume(collision_2)
    assert hi.get(GG) == 1

    hi.consume(collision_3)
    assert hi.get(GG) == 2


def test_simple_median():
    hi = khmer.new_counting_hash(6, 1e6, 2)

    hi.consume("AAAAAA")
    (median, average, stddev) = hi.get_median_count("AAAAAA")
    print median, average, stddev
    assert median == 1
    assert average == 1.0
    assert stddev == 0.0

    hi.consume("AAAAAA")
    (median, average, stddev) = hi.get_median_count("AAAAAA")
    print median, average, stddev
    assert median == 2
    assert average == 2.0
    assert stddev == 0.0

    hi.consume("AAAAAT")
    (median, average, stddev) = hi.get_median_count("AAAAAAT")
    print median, average, stddev
    assert median == 2
    assert average == 1.5
    assert int(stddev * 100) == 50        # .5

    hi.consume("AAAAAT")
    (median, average, stddev) = hi.get_median_count("AAAAAAT")
    print median, average, stddev
    assert median == 2
    assert average == 2.0
    assert stddev == 0.0

    hi.consume("AAAAAT")
    (median, average, stddev) = hi.get_median_count("AAAAAAT")
    print median, average, stddev
    assert median == 3
    assert average == 2.5
    assert int(stddev * 100) == 50        # .5


def test_simple_kadian():
    hi = khmer.new_counting_hash(6, 1e6, 2)
    hi.consume("ACTGCTATCTCTAGAGCTATG")
    assert hi.get_kadian_count("ACTGCTATCTCTAGAGCTATG") == 1

    hi = khmer.new_counting_hash(6, 1e6, 2)
    hi.consume("ACTGCTATCTCTAGAGCTATG")
    hi.consume("ACTGCTATCTCTAGAcCTATG")
    #           ---------------^
    x = hi.get_kadian_count("ACTGCTATCTCTAGAGCTATG")
    assert x == 2, x

    hi = khmer.new_counting_hash(6, 1e6, 2)
    hi.consume("ACTGCTATCTCTAGAGCTATG")
    hi.consume("ACTGCTATCTCTAGAcCTATG")
    #           ---------------^---^
    x = hi.get_kadian_count("ACTGCTATCTCTAGAGCTATG")
    assert x == 2

    hi = khmer.new_counting_hash(6, 1e6, 2)
    hi.consume("ACTGCTATCTCTAGAGCTATG")
    hi.consume("ACTGCTATCTCTAGtGCTAcG")
    #           --------------^^---^
    x = hi.get_kadian_count("ACTGCTATCTCTAGAGCTATG")
    assert x == 1, x


def test_simple_kadian_2():
    hi = khmer.new_counting_hash(6, 1e6, 2)
    hi.consume("ACTGCTATCTCTAGAGCTATG")
    assert hi.get_kadian_count("ACTGCTATCTCTAGAGCTATG") == 1

    hi = khmer.new_counting_hash(6, 1e6, 2)
    hi.consume("ACTGCTATCTCTAGAGCTATG")
    # hi.consume("ACaGCTATCTCTAGAGCTATG")
    hi.consume("ACAGCTATCTCTAGAGCTATG")
    #           --^
    x = hi.get_kadian_count("ACTGCTATCTCTAGAGCTATG")
    assert x == 2, x

    hi = khmer.new_counting_hash(6, 1e6, 2)
    hi.consume("ACTGCTATCTCTAGAGCTATG")
    # hi.consume("ACaGCTATCTCTAGAcCTATG")
    hi.consume("ACAGCTATCTCTAGACCTATG")
    #           --^          --^
    x = hi.get_kadian_count("ACTGCTATCTCTAGAGCTATG")
    assert x == 1, x

    hi = khmer.new_counting_hash(6, 1e6, 2)
    hi.consume("ACTGCTATCTCTAGAGCTATG")
    # hi.consume("ACTGCTATCgCTAGAGCTATG")
    hi.consume("ACTGCTATCGCTAGAGCTATG")
    #                  --^
    x = hi.get_kadian_count("ACTGCTATCTCTAGAGCTATG")
    assert x == 2, x


def test_2_kadian():
    hi = khmer.new_counting_hash(6, 1e6, 2)
    hi.consume("ACTGCTATCTCTAGAGCTATG")
    assert hi.get_kadian_count("ACTGCTATCTCTAGAGCTATG", 2) == 1

    hi = khmer.new_counting_hash(6, 1e6, 2)
    hi.consume("ACTGCTATCTCTAGAGCTATG")
    # hi.consume("ACTGCTATCTCTAGAcCTATG")
    hi.consume("ACTGCTATCTCTAGACCTATG")
    #           ---------------^
    x = hi.get_kadian_count("ACTGCTATCTCTAGAGCTATG", 2)
    assert x == 2, x

    hi = khmer.new_counting_hash(6, 1e6, 2)
    hi.consume("ACTGCTATCTCTAGAGCTATG")
    # hi.consume("ACTGCTATCTCTAGAcCTAtG")
    hi.consume("ACTGCTATCTCTAGACCTATG")
    #           ---------------^---^
    assert hi.get_kadian_count("ACTGCTATCTCTAGAGCTATG", 2) == 2

    hi = khmer.new_counting_hash(6, 1e6, 2)
    hi.consume("ACTGCTATCTCTAGAGCTATG")
    # hi.consume("ACTGCTATCTCTACtcCTAtG")
    hi.consume("ACTGCTATCTCTACTCCTATG")
    #           --------------^^---^
    x = hi.get_kadian_count("ACTGCTATCTCTAGAGCTATG", 2)
    assert x == 2, x

    hi = khmer.new_counting_hash(6, 1e6, 2)
    hi.consume("ACTGCTATCTCTAGAGCTATG")
    # hi.consume("ACTGCTgTCTCTACtcCTAtG")
    hi.consume("ACTGCTGTCTCTACTCCTATG")
    #           ------^-------^^---^
    x = hi.get_kadian_count("ACTGCTATCTCTAGAGCTATG", 2)
    assert x == 1, x


@attr('highmem')
def test_save_load():
    inpath = utils.get_test_data('random-20-a.fa')
    savepath = utils.get_temp_filename('tempcountingsave0.ht')

    sizes = list(PRIMES_1m)
    sizes.append(1000005)

    hi = khmer._new_counting_hash(12, sizes)
    hi.consume_fasta(inpath)
    hi.save(savepath)

    ht = khmer._new_counting_hash(12, sizes)
    ht.load(savepath)

    tracking = khmer._new_hashbits(12, sizes)
    x = hi.abundance_distribution(inpath, tracking)

    tracking = khmer._new_hashbits(12, sizes)
    y = ht.abundance_distribution(inpath, tracking)

    assert sum(x) == 3966, sum(x)
    assert x == y, (x, y)


@attr('highmem')
def test_load_gz():
    inpath = utils.get_test_data('random-20-a.fa')

    savepath = utils.get_temp_filename('tempcountingsave1.ht')
    loadpath = utils.get_temp_filename('tempcountingsave1.ht.gz')

    sizes = list(PRIMES_1m)
    sizes.append(1000005)

    # save uncompressed hashtable.
    hi = khmer._new_counting_hash(12, sizes)
    hi.consume_fasta(inpath)
    hi.save(savepath)

    # compress.
    in_file = open(savepath, 'rb')
    out_file = gzip.open(loadpath, 'wb')
    out_file.writelines(in_file)
    out_file.close()
    in_file.close()

    # load compressed hashtable.
    ht = khmer._new_counting_hash(12, sizes)
    ht.load(loadpath)

    tracking = khmer._new_hashbits(12, sizes)
    x = hi.abundance_distribution(inpath, tracking)

    tracking = khmer._new_hashbits(12, sizes)
    y = ht.abundance_distribution(inpath, tracking)

    assert sum(x) == 3966, sum(x)
    assert x == y, (x, y)


@attr('highmem')
def test_save_load_gz():
    inpath = utils.get_test_data('random-20-a.fa')
    savepath = utils.get_temp_filename('tempcountingsave2.ht.gz')

    sizes = list(PRIMES_1m)
    sizes.append(1000005)

    hi = khmer._new_counting_hash(12, sizes)
    hi.consume_fasta(inpath)
    hi.save(savepath)

    ht = khmer._new_counting_hash(12, sizes)
    ht.load(savepath)

    tracking = khmer._new_hashbits(12, sizes)
    x = hi.abundance_distribution(inpath, tracking)

    tracking = khmer._new_hashbits(12, sizes)
    y = ht.abundance_distribution(inpath, tracking)

    assert sum(x) == 3966, sum(x)
    assert x == y, (x, y)


def test_trim_full():
    hi = khmer.new_counting_hash(6, 1e6, 2)

    hi.consume(DNA)
    hi.consume(DNA)

    seq, pos = hi.trim_on_abundance(DNA, 2)
    assert DNA == seq, seq


def test_trim_short():
    hi = khmer.new_counting_hash(6, 1e6, 2)

    hi.consume(DNA)
    hi.consume(DNA[:50])

    seq, pos = hi.trim_on_abundance(DNA, 2)
    assert DNA[:50] == seq, (seq, pos)
    assert hi.get(seq[-6:]) == 2
    assert hi.get(DNA[:51][-6:]) == 1


def test_maxcount():
    # hashtable should saturate at some point so as not to overflow counter
    kh = khmer.new_counting_hash(4, 4 ** 4, 4)
    kh.set_use_bigcount(False)

    last_count = None
    for i in range(0, 1000):
        kh.count('AAAA')
        c = kh.get('AAAA')

        if c == last_count:
            break
        last_count = c

    assert c != 1000, "should not be able to count to 1000: %d" % c
    assert c == MAX_COUNT, c       # this will depend on HashcountType...


def test_maxcount_with_bigcount():
    # hashtable should not saturate, if use_bigcount is set.
    kh = khmer.new_counting_hash(4, 4 ** 4, 4)
    kh.set_use_bigcount(True)

    last_count = None
    for i in range(0, 1000):
        kh.count('AAAA')
        c = kh.get('AAAA')

        if c == last_count:
            break
        last_count = c

    assert c == 1000, "should be able to count to 1000: %d" % c
    assert c != MAX_COUNT, c


def test_maxcount_with_bigcount_save():
    # hashtable should not saturate, if use_bigcount is set.
    kh = khmer.new_counting_hash(4, 4 ** 4, 4)
    kh.set_use_bigcount(True)

    for i in range(0, 1000):
        kh.count('AAAA')
        c = kh.get('AAAA')

    savepath = utils.get_temp_filename('tempcountingsave.ht')
    kh.save(savepath)

    kh = khmer.new_counting_hash(1, 1, 1)
    kh.load(savepath)

    c = kh.get('AAAA')
    assert c == 1000, "should be able to count to 1000: %d" % c
    assert c != MAX_COUNT, c


def test_bigcount_save():
    # hashtable should not saturate, if use_bigcount is set.
    kh = khmer.new_counting_hash(4, 4 ** 4, 4)
    kh.set_use_bigcount(True)

    savepath = utils.get_temp_filename('tempcountingsave.ht')
    kh.save(savepath)

    kh = khmer.new_counting_hash(1, 1, 1)
    kh.load(savepath)

    # set_use_bigcount should still be True after load (i.e. should be saved)

    assert kh.get('AAAA') == 0

    for i in range(0, 1000):
        kh.count('AAAA')
        kh.get('AAAA')

    assert kh.get('AAAA') == 1000


def test_nobigcount_save():
    kh = khmer.new_counting_hash(4, 4 ** 4, 4)
    # kh.set_use_bigcount(False) <-- this is the default

    savepath = utils.get_temp_filename('tempcountingsave.ht')
    kh.save(savepath)

    kh = khmer.new_counting_hash(1, 1, 1)
    kh.load(savepath)

    # set_use_bigcount should still be True after load (i.e. should be saved)

    assert kh.get('AAAA') == 0

    for i in range(0, 1000):
        kh.count('AAAA')
        kh.get('AAAA')

    assert kh.get('AAAA') == MAX_COUNT


@attr('highmem')
def test_bigcount_abund_dist():
    kh = khmer.new_counting_hash(18, 1e7, 4)
    tracking = khmer.new_hashbits(18, 1e7, 4)
    kh.set_use_bigcount(True)

    seqpath = utils.get_test_data('test-abund-read-2.fa')

    kh.consume_fasta(seqpath)

    dist = kh.abundance_distribution(seqpath, tracking)
    print kh.get('GGTTGACGGGGCTCAGGG')

    pdist = [(i, dist[i]) for i in range(len(dist)) if dist[i]]
    assert dist[1001] == 1, pdist


@attr('highmem')
def test_bigcount_abund_dist_2():
    kh = khmer.new_counting_hash(18, 1e7, 4)
    tracking = khmer.new_hashbits(18, 1e7, 4)
    kh.set_use_bigcount(True)

    seqpath = utils.get_test_data('test-abund-read.fa')

    kh.consume_fasta(seqpath)
    for i in range(1000):
        kh.count('GGTTGACGGGGCTCAGGG')

    dist = kh.abundance_distribution(seqpath, tracking)
    print kh.get('GGTTGACGGGGCTCAGGG')

    pdist = [(i, dist[i]) for i in range(len(dist)) if dist[i]]
    assert dist[1001] == 1, pdist


def test_bigcount_overflow():
    kh = khmer.new_counting_hash(18, 1e7, 4)
    kh.set_use_bigcount(True)

    for i in range(0, 70000):
        kh.count('GGTTGACGGGGCTCAGGG')

    assert kh.get('GGTTGACGGGGCTCAGGG') == MAX_BIGCOUNT


def test_get_ksize():
    kh = khmer.new_counting_hash(22, 1, 1)
    assert kh.ksize() == 22


def test_get_hashsizes():
    kh = khmer.new_counting_hash(22, 100, 4)
    assert kh.hashsizes() == [101, 103, 107, 109], kh.hashsizes()

# def test_collect_high_abundance_kmers():
#    seqpath = utils.get_test_data('test-abund-read-2.fa')
#
#    kh = khmer.new_counting_hash(18, 1e6, 4)
#    hb = kh.collect_high_abundance_kmers(seqpath, 2, 4)


####


def test_load_notexist_should_fail():
    savepath = utils.get_temp_filename('tempcountingsave0.ht')

    hi = khmer.new_counting_hash(12, 1000)
    try:
        hi.load(savepath)
        assert 0, "load should fail"
    except IOError, e:
        print str(e)


def test_load_truncated_should_fail():
    inpath = utils.get_test_data('random-20-a.fa')
    savepath = utils.get_temp_filename('tempcountingsave0.ht')

    hi = khmer.new_counting_hash(12, 1000)
    hi.consume_fasta(inpath)
    hi.save(savepath)

    fp = open(savepath, 'rb')
    data = fp.read()
    fp.close()

    fp = open(savepath, 'wb')
    fp.write(data[:1000])
    fp.close()

    hi = khmer.new_counting_hash(12, 1)
    try:
        hi.load(savepath)
        assert 0, "load should fail"
    except IOError, e:
        print str(e)


def test_load_gz_notexist_should_fail():
    savepath = utils.get_temp_filename('tempcountingsave0.ht.gz')

    hi = khmer.new_counting_hash(12, 1000)
    try:
        hi.load(savepath)
        assert 0, "load should fail"
    except IOError, e:
        print str(e)


def test_load_gz_truncated_should_fail():
    inpath = utils.get_test_data('random-20-a.fa')
    savepath = utils.get_temp_filename('tempcountingsave0.ht.gz')

    hi = khmer.new_counting_hash(12, 1000)
    hi.consume_fasta(inpath)
    hi.save(savepath)

    fp = open(savepath, 'rb')
    data = fp.read()
    fp.close()

    fp = open(savepath, 'wb')
    fp.write(data[:1000])
    fp.close()

    hi = khmer.new_counting_hash(12, 1)
    try:
        hi.load(savepath)
        assert 0, "load should fail"
    except IOError, e:
        print str(e)


def test_counting_file_version_check():
    ht = khmer.new_counting_hash(12, 1, 1)

    inpath = utils.get_test_data('badversion-k12.kh')

    try:
        ht.load(inpath)
        assert 0, "this should fail"
    except IOError, e:
        print str(e)


def test_counting_gz_file_version_check():
    ht = khmer.new_counting_hash(12, 1, 1)

    inpath = utils.get_test_data('badversion-k12.kh.gz')

    try:
        ht.load(inpath)
        assert 0, "this should fail"
    except IOError, e:
        print str(e)


def test_counting_file_type_check():
    inpath = utils.get_test_data('goodversion-k12.ht')

    kh = khmer.new_counting_hash(12, 1, 1)

    try:
        kh.load(inpath)
        assert 0, "this should fail"
    except IOError, e:
        print str(e)


def test_counting_gz_file_type_check():
    ht = khmer.new_hashbits(12, 1, 1)

    inpath = utils.get_test_data('goodversion-k12.ht.gz')

    kh = khmer.new_counting_hash(12, 1, 1)

    try:
        kh.load(inpath)
        assert 0, "this should fail"
    except IOError, e:
        print str(e)


def test_counting_bad_primes_list():
    try:
        ht = khmer._new_counting_hash(12, ["a", "b", "c"], 1)
        assert 0, "bad list of primes should fail"
    except TypeError, e:
        print str(e)


def test_bad_use_bigcount():
    countingtable = khmer.new_counting_hash(4, 4 ** 4, 4)
    countingtable.set_use_bigcount(True)
    assert countingtable.get_use_bigcount()
    try:
        countingtable.get_use_bigcount(True)
        assert 0, "this should fail"
    except TypeError, err:
        print str(err)


def test_consume_absentfasta():
    countingtable = khmer.new_counting_hash(4, 4 ** 4, 4)
    try:
        countingtable.consume_fasta("absent_file.fa")
        assert 0, "This should fail"
    except IOError, err:
        print str(err)


def test_consume_absentfasta_with_reads_parser():
    countingtable = khmer.new_counting_hash(4, 4 ** 4, 4)
    try:
        countingtable.consume_fasta_with_reads_parser()
        assert 0, "this should fail"
    except TypeError, err:
        print str(err)
    readparser = ReadParser(utils.get_test_data('empty-file'))
    try:
        countingtable.consume_fasta_with_reads_parser(readparser)
        assert 0, "this should fail"
    except IOError, err:
        print str(err)


def test_badconsume():
    countingtable = khmer.new_counting_hash(4, 4 ** 4, 4)
    try:
        countingtable.consume()
        assert 0, "this should fail"
    except TypeError, err:
        print str(err)
    try:
        countingtable.consume("AAA")
        assert 0, "this should fail"
    except ValueError, err:
        print str(err)


def test_get_badmin_count():
    countingtable = khmer.new_counting_hash(4, 4 ** 4, 4)
    try:
        countingtable.get_min_count()
        assert 0, "this should fail"
    except TypeError, err:
        print str(err)
    try:
        countingtable.get_min_count("AAA")
        assert 0, "this should fail"
    except ValueError, err:
        print str(err)


def test_get_badmax_count():
    countingtable = khmer.new_counting_hash(4, 4 ** 4, 4)
    try:
        countingtable.get_max_count()
        assert 0, "this should fail"
    except TypeError, err:
        print str(err)
    try:
        countingtable.get_max_count("AAA")
        assert 0, "this should fail"
    except ValueError, err:
        print str(err)


def test_get_badmedian_count():
    countingtable = khmer.new_counting_hash(4, 4 ** 4, 4)
    try:
        countingtable.get_median_count()
        assert 0, "this should fail"
    except TypeError, err:
        print str(err)
    try:
        countingtable.get_median_count("AAA")
        assert 0, "this should fail"
    except ValueError, err:
        print str(err)


def test_get_badkadian_count():
    countingtable = khmer.new_counting_hash(4, 4 ** 4, 4)
    try:
        countingtable.get_kadian_count()
        assert 0, "this should fail"
    except TypeError, err:
        print str(err)
    try:
        countingtable.get_kadian_count("AAA")
        assert 0, "this should fail"
    except ValueError, err:
        print str(err)


def test_badget():
    countingtable = khmer.new_counting_hash(4, 4 ** 4, 4)
    try:
        countingtable.get()
        assert 0, "this should fail"
    except TypeError, err:
        print str(err)


def test_badget_2():
    countingtable = khmer.new_counting_hash(6, 1e6)

    countingtable.consume(DNA)

    assert countingtable.get("AGCTTT") == 1

    assert countingtable.get("GATGAG") == 0

    try:
        countingtable.get("AGCTT")
        assert 0, "this should fail"
    except ValueError, err:
        print str(err)


def test_badtrim():
    countingtable = khmer.new_counting_hash(6, 1e6, 2)

    countingtable.consume(DNA)
    try:
        countingtable.trim_on_abundance()
        assert 0, "this should fail"
    except TypeError, err:
        print str(err)
    countingtable.trim_on_abundance("AAAAAA", 1)


def test_badfasta_count_kmers_by_position():
    countingtable = khmer.new_counting_hash(4, 4 ** 4, 4)
    try:
        countingtable.fasta_count_kmers_by_position()
    except TypeError, err:
        print str(err)

    filename = utils.get_test_data("test-short.fa")
    try:
        countingtable.fasta_count_kmers_by_position(filename, -1, 0)
        assert 0, "this should fail"
    except ValueError, err:
        print str(err)
    try:
        countingtable.fasta_count_kmers_by_position(filename, 0, -1)
        assert 0, "this should fail"
    except ValueError, err:
        print str(err)


def test_badload():
    countingtable = khmer.new_counting_hash(4, 4 ** 4, 4)
    try:
        countingtable.load()
        assert 0, "this should fail"
    except TypeError, err:
        print str(err)


def test_badsave():
    countingtable = khmer.new_counting_hash(4, 4 ** 4, 4)
    try:
        countingtable.save()
        assert 0, "this should fail"
    except TypeError, err:
        print str(err)


def test_badksize():
    countingtable = khmer.new_counting_hash(4, 4 ** 4, 4)
    try:
        countingtable.ksize(True)
        assert 0, "this should fail"
    except TypeError, err:
        print str(err)


def test_badhashsizes():
    countingtable = khmer.new_counting_hash(4, 4 ** 4, 4)
    try:
        countingtable.hashsizes(True)
        assert 0, "this should fail"
    except TypeError, err:
        print str(err)


def test_badconsume_and_tag():
    countingtable = khmer.new_counting_hash(4, 4 ** 4, 4)
    try:
        countingtable.consume_and_tag()
        assert 0, "this should fail"
    except TypeError, err:
        print str(err)


def test_consume_fasta_and_tag():
    countingtable = khmer.new_counting_hash(4, 4 ** 4, 4)
    try:
        countingtable.consume_fasta_and_tag()
        assert 0, "this should fail"
    except TypeError, err:
        print str(err)
    countingtable.consume_fasta_and_tag(utils.get_test_data("test-graph2.fa"))
