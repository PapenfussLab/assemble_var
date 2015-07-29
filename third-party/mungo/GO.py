"""
GO module
"""

from GoTools import *


slimFilename = '/Users/papenfuss/databases/GO/goslim_generic.obo'

filename='/Users/papenfuss/databases/GO/gene_ontology_edit.obo'
onto = Ontology.load(filename)

biological_process = onto.findName('biological_process')
molecular_function = onto.findName('molecular_function')
cellular_component = onto.findName('cellular_component')

apoptosis = onto.findName('apoptosis')


def childrenOfBiologicalProcess():
    print "Children of biological_process:"
    for gid in biological_process.children:
        print onto[gid]
        print


def testTraversals(gid):
    print 'LR traversal'
    traversal = onto.traverse(gid)
    for level in traversal:
        print level
    
    print '-'*60
    print 'Traversal iterator'
    for level in onto.traverseIter(gid):
        print level
    
    print '-'*60
    print 'Traversal iterator'
    for level in onto.traverseIter2(gid):
        print level


def test():
    gid = apoptosis.id
    subset = onto.traverseToGivenNodes(gid, biological_process.children)
    subset.sort()
    print subset
    print
    for gid in subset:
        print onto[gid]
        print

if __name__=='__main__':
    gid = 'GO:0000001'
    testTraversals(gid)
    print
    print 'Rank:', onto.rank(gid)
    
