#! /usr/bin/env python2
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
import khmer
import sys
import screed
import graph_tool.all as gt

input_fasta = sys.argv[3]
K = int(sys.argv[1])
x = float(sys.argv[2])


ht = khmer.new_hashbits(K, x, 4)

sparse_graph = gt.Graph()
hashes = sparse_graph.new_vertex_property("long long")


for n, record in enumerate(screed.open(input_fasta)):
    if n % 1000 == 0:
        print >>sys.stderr, '...loaded and tagged {} sequences'.format(n)
    name = record.name
    sequence = record.sequence

    ht.consume_sequence_and_tag_with_labels(sequence, n)
    tags = ht.sweep_tag_neighborhood(sequence, 0)
    for i in xrange(len(tags) - 1):
        src = tags[i]
        dst = tags[i + 1]

        new = False

        srcv = gt.find_vertex(sparse_graph, hashes, src)
        if not srcv:
            srcv = sparse_graph.add_vertex()
            hashes[srcv] = src
            new = True
        else:
            srcv = srcv[0]

        dstv = gt.find_vertex(sparse_graph, hashes, dst)
        if not dstv:
            dstv = sparse_graph.add_vertex()
            hashes[dstv] = dst
            new = True
        else:
            dstv = dstv[0]

        if new:
            e = sparse_graph.add_edge(srcv, dstv)

print 'Sparse graph has {} nodes, {} edges'.format(sparse_graph.num_vertices(), sparse_graph.num_edges())
comp = gt.label_largest_component(sparse_graph, directed=False)
#pos = gt.radial_tree_layout(sparse_graph, sparse_graph.vertex(0))
gt.graph_draw(sparse_graph, output_size=(
    5000, 5000), output=input_fasta + '_sparse.png')
sparse_graph.set_vertex_filter(comp)
gt.graph_draw(sparse_graph, output_size=(
    5000, 5000), output=input_fasta + '_sparse_comp.png')
