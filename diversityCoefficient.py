from __future__ import division
from misc import takeNodesInTree

#Returns the diversity coefficient and the list of (node (name,rank),number of assignations in this group of samples) pairs
def computeDiversityCoefficient(numberNodesInTree,sampleNameList,dataArray):
    #@sample is a list of (name of node,rank of node,sampleHitList)
    sample,numberTotalAssignments,numberNodes = takeNodesInTree(dataArray[7],sampleNameList)
    resultNodes = []
    for node in sample:
        assignments = 0
        if not (len(node) == 3):
            print " \n/!\ ERROR: List of samples wrong:",len(node),"."
            raise ValueError
        name,rank,sampleList = node[0],node[1],node[2]
        for hit in sampleList:
            if not (len(hit) == 2):
                print "\n/!\ ERROR: Sample Hit List wrong:",len(hit),"."
                raise ValueError
            assignments += hit[1]
        resultNodes.append(((name,rank),assignments))
    if not numberNodes:
        print "\n/!\ ERROR: numberNodes is null."
        raise ValueError
    assignmentsMean = numberTotalAssignments/numberNodes
    if not numberNodesInTree or not numberNodes or not assignmentsMean:
        print "\n/!\ ERROR: Taxonomic Tree is empty: whole tree:",numberNodesInTree,"reduced tree:",numberNodes,"assignmentsMean:",assignmentsMean,"."
        raise ValueError
    return numberTotalAssignments/(numberNodesInTree*assignmentsMean),resultNodes
