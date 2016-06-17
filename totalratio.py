from __future__ import division

from taxoTree import TaxoTree,printTree
from parsingMatrix import parseMatrix
from parsingTree import parseTree
from misc import mem

#Returns a list of (name,rank,sampleHitList) tuples associated with
#nodes of the taxonomic tree reduced to the sample sampleName
#and the number of assignments in this tree
#NB: We need to keep the whole sampleHitList (and not only the
#number associated with sampleName) in order to apply set operations (intersection, ...)
def inSample(element,sampleNameList):
    #element[1] (number associated to a sample) must be non-zero
    assert element[1]
    #If list is empty
    if not sampleNameList:
        return False
    for name in sampleNameList:
        if (element[0] == name):
            return True
    return False

#@sampleName is a list of names of samples
#takeNodesInTree(tree,sampleName) should return the taxonomic tree reduced to the nodes assigned in sample sampleName = the list of assigned node + their sampleHitList, and the number of assignments in the reduced tree, and the number of nodes in the reduced tree
def takeNodesInTree(tree,sampleNameList):
    #@sampleHitList (see TaxoTree) is a list of (name of sample, number) pairs attached to a node of a TaxoTree
    sample = []
    numberTotalAssignments = 0
    numberNodes = 0
    queue = [ tree ]
    while queue:
        node = queue.pop()
        isInSample = []
        for x in node.sampleHitList:
            if inSample(x,sampleNameList):
                isInSample.append(x)
        #if node is in sample, ie isInSample is not empty
        if isInSample:
            sample.append((node.name,node.rank,node.sampleHitList))
            numberTotalAssignments += isInSample[0][1]
            numberNodes += 1
        queue += node.children
    return sample,numberTotalAssignments,numberNodes

#Returns boolean and sampleHitList if true
def memAndSampleHitList(x,nodeList):
    sampleHitList = []
    nodeListCopy = []
    for nd in nodeList:
        nodeListCopy.append(nd)
    #While @nodeList is not empty and @sampleHitList is empty
    while nodeListCopy:
        node = nodeListCopy.pop()
        if (x[0] == node[0] and x[1] == node[1]):
            return True,node[2]
    return False,[]

#Returns @common,@in1,@in2,@numberA1,@numberA2
#@numberA1 is the number of assignments in the tree reduced to nodes assigned in sampleList #1 (same goes for @numberA2)
#@numberN1 is the number of nodes in the reduced tree
#@common is the set of common nodes (name,rank,sampleHitList (from 1),sampleHitList (from 2)), @in1 the set of nodes only assigned in 1, @in2 the set of nodes only assigned in 2
def compute(tree,sampleNameList1,sampleNameList2):
    nodeList1,numberAT1,numberN1 = takeNodesInTree(tree,sampleNameList1)
    nodeList2,numberAT2,numberN2 = takeNodesInTree(tree,sampleNameList2)
    common, in1, in2 = [],[],[]
    #@numberC is the number of nodes in common
    numberC = 0
    for node in nodeList1:
        boolean,sampleHitList = memAndSampleHitList(node,nodeList2)
        #If @sampleHitList is not empty
        #That is if node also belongs to the second list of samples
        if boolean:
            common += [(node[0],node[1],node[2],sampleHitList)]
            numberC += 1
        else:
            in1 += [node]
    for node in nodeList1:
        #@boolean answers: Is it a node that does not belong to the first list of samples?
        boolean = True
        for cNode in common:
            if (node[0] == cNode[0] and node[1] == cNode[1]):
                boolean = False
        if boolean:
            in2 += [node]
    return common,in1,in2,numberAT1,numberAT2,numberN1,numberN2,numberC

def countAssignments(in12,sampleNameList):
    s = 0
    for node in in12:
        sampleHitList = node[2]
        for sampleHit in sampleHitList:
            if mem(sampleHit[0],sampleNameList):
                s += sampleHit[1]
    return s

def countAssignmentsInCommon(common,sampleNameList1,sampleNameList2):
    s = 0
    for node in common:
        sampleHitList1 = node[2]
        sampleHitList2 = node[3]
        for sampleHit in sampleHitList1:
            if mem(sampleHit[0],sampleNameList1):
                s += sampleHit[1]
        for sampleHit in sampleHitList2:
            if mem(sampleHit[0],sampleNameList2):
                s += sampleHit[1]
    return s

#When totalRatio is close to 1, it means the two samples are quite alike
def totalRatioNormalized(commonA,numberA1,numberA2):
    #if @numberA1 = @numberA2 = @commonA = 0
    if not numberA1 and not numberA2 and not commonA:
        return "+inf"
    return (commonA/(numberA1 + commonA + numberA2))

def totalRatio(commonA,numberA1,numberA2):
    return (numberA1 + numberA2)

def diffRatioNormalized(commonA,numberA1,numberA2):
    #if @numberA1 = @common = @numberA2 = 0
    if not numberA1 and not numberA2 and not commonA:
        return "+inf"
    return ((numberA1 + numberA2)/(numberA1 + commonA + numberA2))

def diffRatio(commonA):
    return commonA
