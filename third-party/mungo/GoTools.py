"""
GoTools module
"""

from mungoCore import *


class GoTerm:
    """Class representing a GoTerm"""
    
    def __init__(self, goAccession='', name='', namespace='', defn=''):
        self.goAccession = goAccession
        self.name = name
        self.namespace = namespace
        self.defn = defn
        self.is_a = []
        self.children = []
        self.alt_id = []
        self.relationships = []
        self.replaced_by = []
        self.is_obsolete = False
    
    def __repr__(self):
        lines = []
        lines.append('[Term]')
        lines.append('goAccession: %s' % self.goAccession)
        lines.append('name: %s' % self.name)
        lines.append('namespace: %s' % self.namespace)
        lines.append('defn: %s' % self.defn)
        lines.append('is_a: %s' % str(self.is_a))
        lines.append('alt_id: %s' % str(self.alt_id))
        lines.append('children: %s' % str(self.children))
        lines.append('relationships: %s' % (self.relationships))
        lines.append('is_obsolete: %s' % self.is_obsolete)
        lines.append('replaced_by: %s' % str(self.replaced_by))
        return '\n'.join(lines) + '\n'
    
    @staticmethod
    def fromKvList(kvList):
        g = GoTerm()
        trans = {'id':'goAccession', 'def':'defn'}
        for k,v in kvList:
            if k in ['name', 'namespace']:
                g.__dict__[k] = v
            elif k=='is_obsolete' and v=='true':
                g.__dict__[k] = True
            elif k in ['id', 'def']:
                g.__dict__[trans[k]] = v
            elif k in ['is_a', 'alt_id', 'replaced_by', 'relationship']:
                try:
                    g.__dict__[k].append(v)
                except KeyError:
                    g.__dict__[k] = [v]
            elif k in ['synonym', 'subset', 'comment', 'consider', 'xref', 'disjoint_from']:
                pass
            else:
                print k,v
        return g


class Ontology:
    """Ontology class representing complete ontology,
    usually read from an obo file using load."""
    
    def __init__(self):
        "Constructor"
        self.terms = {}
        self.obsolete = {}
        self.alt_id = {}
        self.parents = {}
        self.children = {}
        self.roots = []
        self.leaves = set()
    
    @staticmethod
    def load(filename):
        onto = Ontology()
        for i,g in enumerate(GoFile(filename)):
            onto.add(g)
        onto._findChildren()
        return onto
    
    def add(self, g):
        """Add a GO term to the ontology.
        
        @param g: GoTerm object
        """
        # Ignore obsolete terms
        if g.is_obsolete or 'OBSOLETE' in g.defn:
            self.obsolete[g.goAccession] = g
        else:
            self.terms[g.goAccession] = g
            self.parents[g.goAccession] = g.is_a
            if len(g.is_a)==0:
                self.roots.append(g.goAccession)
        
        if g.alt_id:
            for a in g.alt_id:
                self.alt_id[a] = g.goAccession
    
    def _findChildren(self):
        """Fill in the children attributes of GoTerms across the ontology.
        Called by load."""
        for goAccession,term in self.terms.iteritems():
            for parent in term.is_a:
                try:
                    self.terms[parent].children.append(goAccession)
                except:
                    print 'Trouble', goAccession, parent
        
        # Build copy of children & construct leaves
        for goAccession,term in self.terms.iteritems():
            self.children[goAccession] = term.children
            if len(term.children)==0:
                self.leaves.add(goAccession)
        self.leaves = list(self.leaves)
    
    def __getitem__(self, goAccession):
        """Mapping interface
        
        @param goAccession: GO id
        @return: Corresponding GO term
        """
        return self.terms[goAccession]
    
    def findName(self, name):
        """Find the GO description (perfect match only)"""
        for term in self.terms.values():
            if term.name==name:
                return term
    
    def traverse(self, goAccession):
        "Traverse ontology from goAccession to root"
        if type(goAccession)==str:
            traversal = [set([goAccession])]
        elif type(goAccession)==list:
            traversal = [set(goAccession)]
        else:
            raise Exception('goAccession type not supported')
        
        i = 0
        while True:
            j = 0
            traversal.append(set())
            for node in traversal[i]:
                for p in self.terms[node].is_a:
                    traversal[-1].add(p)
            if len(traversal[-1])==0:
                del traversal[-1]
                break
            i += 1
        return traversal
    
    def traverse_old(self, goAccession):
        "Traverse ontology from goAccession to root"
        traversal = []
        if type(goAccession)==str:
            traversal.append([goAccession])
        elif type(goAccession)==list:
            traversal.append(goAccession)
        else:
            raise Exception('goAccession type not supported')
        
        i = 0
        while True:
            j = 0
            traversal.append([])
            for node in traversal[i]:
                parents = self.terms[node].is_a
                if parents:
                    for p in parents:
                        traversal[-1].append(p)
                else:
                    traversal[-1].append(node)
            i += 1
            if traversal[-1]==traversal[-2]:
                del traversal[-1]
                break
        return traversal
    
    def rank(self, goAccession):
        """Minumum distance to root"""
        traversal = self.traverse(goAccession)
        root = traversal[-1][0]
        for i,level in enumerate(traversal):
            if root in level:
                return i
    
    def traverseIter(self, goAccession):
        "Traverse ontology from goAccession to root"
        if type(goAccession)==str:
            nodes = [goAccession]
        elif type(goAccession)==list:
            nodes = goAccession
        else:
            raise Exception('goAccession type not supported')
        yield nodes
        
        while True:
            parentNodes = []
            for node in nodes:
                parents = self.terms[node].is_a
                if parents:
                    for p in parents:
                        parentNodes.append(p)
                else:
                    parentNodes.append(node)
            if nodes==parentNodes:
                raise StopIteration
            yield parentNodes
            nodes = copy.copy(parentNodes)
    
    def traverseToGivenNodes(self, goAccession, givenNodes):
        "Traverse ontology from goAccession to given nodes"
        included = []
        for nodes in self.traverseIter(goAccession):
            for node in nodes:
                if node in givenNodes:
                    included.append(node)
        return included


def GoFile(iFileHandle):
    """Factory function interface to GoReader
    
    @param iFileHandle: filename or file object
    """
    return GoReader(iFileHandle)


class States:
    """Labels for stateful parser"""
    BEGIN = 0
    TERM = 1
    TYPEDEF = 2


class GoReader(AbstractDataReader):
    """GoReader class"""
    
    def _generator(self):
        kvList = []
        state = States.BEGIN
        for line in self.iFile:
            line = line.strip()
            if not line:
                continue
            elif state==States.BEGIN:
                if line=='[Term]':
                    state = States.TERM
                continue
            
            if state==States.TERM:
                if line=='[Term]':
                    if kvList:
                        yield GoTerm.fromKvList(kvList)
                    kvList = []
                elif line=='[Typedef]':
                    if kvList:
                        yield GoTerm.fromKvList(kvList)
                    state = States.TYPEDEF
                    kvList = []
                else:
                    try:
                        tokens = line.split(': ')
                        key = tokens[0]
                        value = ': '.join(tokens[1:])
                        value = value.split('!')[0].rstrip()
                    except Exception, e:
                        print >> sys.stderr, e
                        print >> sys.stderr, line
                        sys.exit(-1)
                    kvList.append((key, value))
        if kvList:
            print kvlist
            yield GoTerm.fromKvList(kvList)
