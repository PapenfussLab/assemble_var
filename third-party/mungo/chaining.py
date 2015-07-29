"""
chaining module
"""

import sys, copy


def hspSort(hsps):
    hsps.sort(key=lambda h: (h.subjectId, h.strand(), h.sStart))
    return hsps


def findAllChains(hsps):
    if len(hsps)==0:
        return None
    elif len(hsps)==1:
        return [hsps]
    else:
        hsps = hspSort(hsps)
        chains = [[hsps[0]]]
        i = 1
        for h in hsps[1:]:
            if h.subjectId!=hsps[i-1].subjectId \
             or h.strand()!=hsps[i-1].strand() \
             or h.sStart-hsps[i-1].sStart>50000 \
             or ((h.qStart<=hsps[i-1].qStart and h.strand()=='+') \
             or (h.qStart>=hsps[i-1].qStart and h.strand()=='-')):
                chains.append([h])
            else:
                chains[-1].append(h)
            i += 1
        return chains


def scoreChain(chain):
    return sum([hsp.bitScore for hsp in chain])


def findBestChain(hsps):
    chains = findAllChains(hsps)
    scores = [scoreChain(chain) for chain in chains]
    chains = zip(scores, chains)
    chains.sort()
    bestChain = chains[-1][1]
    return bestChain
