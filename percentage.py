from __future__ import division
import numpy as np
import re

from misc import mem
from taxoTree import TaxoTree
from parsingInfo import parseInfo
from parsingTree import parseTree

integer = re.compile("[0-9]+")

#memPositions is the list of positions in the infoList of information that belongs to metadataList
#In the order of the elements of metadataList
def memPositionsInfo(infoList,metadataList):
    positions = []
    for datum in metadataList:
        i = 0
        n = len(infoList)
        while i < n and not (datum == infoList[i]):
            i += 1
        if not (i == n):
            positions.append(i)
    return positions

def isInInterval(value,interval1,interval2):
    if interval2 == "+inf" and interval1 == "-inf":
        return True
    elif interval1 == "+inf" or interval2 == "-inf":
        return False
    elif interval1 == "-inf":
        return value <= interval2
    elif interval2 == "+inf":
        return value >= interval1
    else:
        return (value <= interval2 and value >= interval1)

#Computes the samples that match the requirements for the metadata group
#@samplesListInGroup is the list of lists of sample ID that match the requirements for a certain metadatum
#@interval1List and @interval2List are the lists of lower and upper interval bounds for the metadata
def computeSamplesInGroup(samplesInfoList,infoList,metadataList,interval1List,interval2List):
    samplesListInGroup = []
    infoPositions = memPositionsInfo(infoList,metadataList)[::-1]
    #@datumPos matches the position of datum corresponding to the information at @infoPos
    datumPos = 0
    while infoPositions:
        #@samplesInGroup contains the samples' ID that match the requirements for the metadatum at infoPos position
        samplesInGroup = []
        infoPos = infoPositions.pop()
        assert (datumPos < len(interval1List))
        for sample in samplesInfoList:
            assert (infoPos < len(sample))
            #asserted in parseInfo that len(sample) == len(infoList)
            #if strict equality with the default value is required
            if (interval1List[datumPos] == interval2List[datumPos]):
                #If the requirements are fulfilled
                if (sample[infoPos] == "N" and interval1List[datumPos] == "N" and interval2List[datumPos] == "N"):
                    #Adds the sampleID to the list
                    samplesInGroup += [sample[0]]
                #if sample[infoPos] == interval1List[datumPos] == interval2List[datumPos]
                elif not (sample[infoPos] == "N") and (int(sample[infoPos]) == interval1List[datumPos]):
                    #Adds the sampleID to the list
                    samplesInGroup += [sample[0]]
                #In any other case, the sample is not added
            #If the length of interval of acceptable values is non-zero
            else:
                if integer.match(sample[infoPos]) and isInInterval(int(sample[infoPos]),interval1List[datumPos],interval2List[datumPos]):
                    samplesInGroup += [sample[0]]
        datumPos += 1
        samplesListInGroup.append(samplesInGroup)
    return samplesListInGroup

#___________________________________________________________

#Symmetrically, we need to compute the positions in the sample of the assignments to the group of bacterias chosen
#@nodesGroup contains (name,rank) of bacterias chosen
def memPositionsNodes(speciesList,nodesGroup):
    positions = []
    for node in nodesGroup:
        n = len(speciesList)
        i = 0
        while (i < n) and not (speciesList[i][0] == node[0] and speciesList[i][1] == node[1]):
            i += 1
        positions.append(i)
    return positions

def memReturnSampleHitList(sampleID,samplesList):
    for sample in samplesList:
        if (sample[0] == sampleID):
            return sample
    return []

def sumList(ls):
    s = 0
    #Deletes sample Name
    for i in ls[1:]:
        s += i
    return s

#Returns the list of percentages of assignment to the group of bacterias in the samples corresponding to a certain metadatum (@samplesInGroup)
def computePercentageAssignmentNodes(samplesListInGroup,nodesGroup,samplesList,speciesList):
    percentageList = []
    nodePositions = memPositionsNodes(speciesList,nodesGroup)
    samplesListCopy = []
    for s in samplesListInGroup:
        samplesListCopy.append(s)
    while samplesListCopy:
        #@samplesInGroup is a list of sample ID corresponding to a certain metadatum
        samplesInGroup = samplesListCopy.pop()
        #@assignmentsInGroup contains the number of assignment to the nodes in nodesPositions in the samples in samplesInGroup
        assignmentsInGroup = 0
        #@totalAssignments contains the total number of assignments in the samples in samplesInGroup
        totalAssignments = 0
        #Computes the total number of assignments
        for sample in samplesInGroup:
            totalAssignments += sumList(memReturnSampleHitList(sample,samplesList))
        #Computes the number of assignments for the chosen group of bacterias
        for nodePos in nodePositions:
            for sample in samplesInGroup:
                sampleHitList = memReturnSampleHitList(sample,samplesList)
                #Adds for every sample the number corresponding the node in nodePositions
                #if @sampleHitList is non-empty
                if nodePos <= len(sampleHitList):
                    assignmentsInGroup += sampleHitList[nodePos]
                else:
                    print "\n/!\ ERROR: [BUG] [percentage/computePercentageAssignmentNodes] List out of range: nodePos:",nodePos,"length of sampleHitList:",len(sampleHitList),"sampleHitList:",sampleHitList
                    raise ValueError
        if totalAssignments:
            percentageList.append(100*assignmentsInGroup/totalAssignments)
        else:
            percentageList.append("+inf")
    return percentageList

#___________________________________________________________
    
#A second version of this algorithm is only giving the root of the subtrees desired, and computing the number of assignments to the whole subtree.
#@memPositionsSubtree computes the list of subtrees that one needs to check out to get the number of assignments.
#@nodesGroup contains (name,rank) of bacterias chosen
def memPositionsSubtree(nodesGroup,tree):
    subtreeList = []
    for node in nodesGroup:
        subtreeList.append(tree.search(node[0],node[1]))
    return subtreeList

#Returns the total number of assignments in tree + the number of assignments associated to samples contained in samplesInGroup
def countAssignmentsInTree(tree,samplesInGroup):
    totalNumber = 0
    assignmentsNumber = 0
    currNode = TaxoTree("None")
    queue = [tree]
    while queue:
        node = queue.pop()
        for sample in node.sampleHitList:
            totalNumber += sample[1]
            if mem(sample[0],samplesInGroup):
                assignmentsNumber += sample[1]
        queue = []
        for ch in node.children:
            queue.append(ch)
    return totalNumber,assignmentsNumber
    
#Returns the list of percentages of assignment to the group of bacterias with second version
def computePercentageAssignmentTree(samplesListInGroup,nodesGroup,tree):
    subtreeList = memPositionsSubtree(nodesGroup,tree)
    percentageList = []
    for samplesInGroup in samplesListInGroup:
        totalAssign = 0
        assignInGroup = 0
        for subtree in subtreeList:
            totalAssignments,assignmentsInGroup = countAssignmentsInTree(subtree,samplesInGroup)
            totalAssign += totalAssignments
            assignInGroup += assignmentsInGroup
        if totalAssign:
            percentageList.append(100*assignInGroup/totalAssign)
        else:
            percentageList.append("+inf")
    return percentageList

#__________________________________________________________

#Returns an array giving the percentage of assignment to a certain group of bacterias depending on the metadata
#if samplesListInGroup is of length N, then the function returns a matrix of dimensions 1xN
def percentageAssign(samplesInfoList,infoList,samplesListInGroup,tree,nodesGroup,samplesList,speciesList,usingTree):
    #the elements in @percentageList are in the order of elements in @metadataList, see computeSamplesInGroup
    if usingTree:
        percentageList = computePercentageAssignmentTree(samplesListInGroup,nodesGroup,tree)[::-1]
    else:
        percentageList = computePercentageAssignmentNodes(samplesListInGroup,nodesGroup,samplesList,speciesList)[::-1]
    if not (len(percentageList) == len(samplesListInGroup)):
        print "\n/!\ ERROR: [BUG] [percentage/percentageAssign] The number of rows does not match the number of groups of samples."
        raise ValueError
    n = len(samplesListInGroup)
    result = np.zeros(n)
    for i in range(n):
        percent = percentageList.pop()
        result[i] = percent
    return result
