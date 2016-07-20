from __future__ import division

from taxoTree import TaxoTree,printTree
from parsingMatrix import parseMatrix
from parsingTree import parseTree
from misc import inSample,takeNodesInTree,memAndSampleHitList

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
            if (sampleHit[0] in sampleNameList):
                s += sampleHit[1]
    return s

def countAssignmentsInCommon(common,sampleNameList1,sampleNameList2):
    s = 0
    for node in common:
        sampleHitList1 = node[2]
        sampleHitList2 = node[3]
        for sampleHit in sampleHitList1:
            if (sampleHit[0] in sampleNameList1):
                s += sampleHit[1]
        for sampleHit in sampleHitList2:
            if (sampleHit[0] in sampleNameList2):
                s += sampleHit[1]
    return s

#When totalRatio is close to 1, it means the two samples are quite alike. [ Please read README to see in which way ]
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
        return "-inf"
    return ((numberA1 + numberA2)/(numberA1 + commonA + numberA2))

def diffRatio(commonA):
    return commonA
