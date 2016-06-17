#To ensure floating division
from __future__ import division

from taxoTree import TaxoTree
from totalratio import inSample,takeNodesInTree
from parsingMatrix import parseMatrix
from parsingTree import parseTree

#Patterns in the taxonomic tree are divided in two categories: common patterns, which are shared by the two trees induced by the samples of each group, and specific patterns, that only appear in one of the two trees. 

#@enumerateCommonPattern gives the list of common patterns to both samples
def enumerateCommonPattern(tree,sampleName1,sampleName2):
    #@nodesList1 is the list of nodes (name,rank,sampleHitList) for sampleName1 (cf. calculus.py)
    nodesList1,_,_ = takeNodesInTree(tree,[sampleName1])
    #print "nodesList1\n"
    #for node in nodesList1:
        #print str(node) + "\n"
    #print "NODESLIST1"
    #print nodesList1
    commonPatternList = []
    specificPattern1List = []
    specificPattern2List = []
    #We consider every node of the tree as a potential root for a pattern (common or specific)
    for node in nodesList1:
        #print "POTENTIAL ROOT"
        #print node[0],node[1]
        patternList = []
        name,rank = node[0],node[1]
        root = tree.search(name,rank)
        #@candidateNodes is the list of TaxoTree nodes that can potentially be added to the pattern
        candidateNodes = [root]
        #@n is the size of the largest common subtree in the tree rooted at node
        n = 0
        while candidateNodes:
            child = candidateNodes.pop()
            #print child.name,child.rank
            sampleHitList = child.sampleHitList
            #Checking if node has been assigned in both samples
            isInSample = []
            for x in sampleHitList:
                if inSample(x,[sampleName1]) and inSample(x,[sampleName2]):
                    isInSample.append(x)
            #If list are not empty, then node has been assigned in both samples
            if len(isInSample):
                #print "is in both samples"
                assert (len(isInSample) == 1)
                #Adds to @n the sum of the two number of assignments to this node
                n += 1 #isInSample1[0][1] + isInSample2[0][1]
                patternList.append((child.name,child.rank))
                candidateNodes += child.children
        #print n
        if (maxPatternLength < n):
            maxPatternList = patternList
        maxPatternLength = max(maxPatternLength,n)
    return maxPatternLength,maxPatternList

def enumerateSpecificPatterns():
    ()

def patternRatio(tree,sampleName1,sampleName2):
    #@numberN1 and numberN2 are respectively the number of nodes in sample1 and sample2
    _,_,numberN1 = takeNodesInTree(tree,[sampleName1])
    _,_,numberN2 = takeNodesInTree(tree,[sampleName2])
    maxPatternLength,patternList = maxPattern(tree,sampleName1,sampleName2)
    return (numberN1 + numberN2 - 2*maxPatternLength),patternList

def patternRatioNormalized(tree,sampleName1,sampleName2):
    #@numberN1 and numberN2 are respectively the number of nodes in sample1 and sample2
    _,_,numberN1 = takeNodesInTree(tree,[sampleName1])
    _,_,numberN2 = takeNodesInTree(tree,[sampleName2])
    maxPatternLength,patternList = maxPattern(tree,sampleName1,sampleName2)
    #if numberN1 = numberN2 = 0
    if not numberN1 and not numberN2:
        return "+inf",[]
    return 2*maxPatternLength/(numberN1 + numberN2),patternList

