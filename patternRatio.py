#To ensure floating division
from __future__ import division

from taxoTree import TaxoTree
from misc import inSample,takeNodesInTree,mergeList,trimList
from parsingMatrix import parseMatrix
from parsingTree import parseTree

#Patterns in the taxonomic tree are divided in two categories: common patterns, which are shared by the two trees induced by the samples of each group, and specific patterns, that only appear in one of the two trees. 

#@enumerateCommonPatterns gives the list of common patterns to both samples
#@commonPatternsList contains (pattern (list of (name,rank) pairs),number of assignments in nodes of this pattern in both samples,number of nodes in this pattern) 
def enumerateCommonPatterns(tree,sampleNameList1,sampleNameList2):
    commonPatternsList = []
    #@nodesList1 is the list of nodes (name,rank,sampleHitList) for sampleNameList1 (cf. misc.py)
    nodesList1,_,_ = takeNodesInTree(tree,sampleNameList1)
    #We consider every node of the tree as a potential root for a pattern
    for node in nodesList1:
        pattern = []
        numberAssignments = 0
        numberNodes = 0
        name,rank = node[0],node[1]
        #Gets the subtree (of the whole taxonomic tree) rooted at node
        root = tree.search(name,rank)
        #@candidateNodes is the list of TaxoTree nodes that can potentially be added to the pattern
        candidateNodes = [root]
        while candidateNodes:
            child = candidateNodes.pop()
            sampleHitList = child.sampleHitList
            #Checking if child has been assigned in both samples
            isInSampleList1 = []
            isInSampleList2 = []
            for x in sampleHitList:
                if inSample(x,sampleNameList1):
                    isInSampleList1.append(x)
                if inSample(x,sampleNameList2):
                    isInSampleList2.append(x)
            #If both lists are not empty, then node has been assigned in both samples
            if isInSampleList1 and isInSampleList2:
                #Merge the elements of both lists, deleting duplicates
                #e.g. if OPNA-J90 belongs to both sampleNameLists it would corrupt the result
                #as the assignments in this sample to child would be duplicated
                #(assuming there is no duplicate in each list)
                isInSample = mergeList(isInSampleList1,isInSampleList2)
                pattern.append((child.name,child.rank))
                numberNodes += 1
                for x in isInSample:
                    numberAssignments += x[1]
                candidateNodes += child.children
        #if the pattern is non-empty
        if pattern:
            commonPatternsList.append((pattern,numberAssignments,numberNodes))
    return commonPatternsList

#___________________________________________________________________________________________________________

#@enumerateSpecificPatterns gives the list of specific patterns in sampleNameListPattern (that is, patterns that belongs to the tree induced by the samples in sampleNameListPattern and not by the samples in sampleNameListOther. In the case one sample belongs to both sampleNameListPattern and sampleNameListOther, it will not taken into account)
#@specificPatternsList contains (pattern (list of (name,rank) pairs),number of assignments in nodes of this pattern in samples of sampleNameListPattern,number of nodes in this pattern) 
def enumerateSpecificPatterns(tree,sampleNameListPattern,sampleNameListOther):
    specificPatternsList = []
    #List from where samples in both lists are deleted and only elements from sampleNameListPattern remain
    sampleNameListPatternTrimmed = trimList(sampleNameListPattern,sampleNameListOther)
    #@nodesList is the list of nodes (name,rank,sampleHitList) for sampleNameListPatternTrimmed
    nodesList,_,_ = takeNodesInTree(tree,sampleNameListPatternTrimmed)
    #Pretty much the same procedure than for @enumerateCommonPatterns
    for node in nodesList:
        pattern = []
        numberAssignments = 0
        numberNodes = 0
        name,rank = node[0],node[1]
        root = tree.search(name,rank)
        candidateNodes = [root]
        while candidateNodes:
            child = candidateNodes.pop()
            sampleHitList = child.sampleHitList
            isInSampleListPattern = []
            for x in sampleHitList:
                if inSample(x,sampleNameListPatternTrimmed):
                    isInSampleListPattern.append(x)
                if inSample(x,sampleNameListOther):
                    #This node is assigned in the samples of sampleNameListOther
                    #so it is discarded from the pattern
                    isInSamplePattern = []
                    break
            if isInSampleListPattern:
                pattern.append((child.name,child.rank))
                numberNodes += 1
                for x in isInSampleListPattern:
                    numberAssignments += x[1]
                candidateNodes += child.children
        #if the pattern is non-empty
        if pattern:
            specificPatternsList.append((pattern,numberAssignments,numberNodes))
    return specificPatternsList

#____________________________________________________________________________________________________________

def patternRatio(commonPatternsList,specificPatternsList1,specificPatternsList2):
    #There will be no assignment counted more than once, since patterns are disjoint
    commonAssignments = 0
    for x in commonPatternsList:
        if (len(x[0]) > 1):
            commonAssignments += x[1]
    specificAssignments1 = 0
    for x in specificPatternsList1:
        specificAssignments1 += x[1]
    specificAssignments2 = 0
    for x in specificPatternsList2:
        specificAssignments2 += x[1]
    if not specificAssignments1 and not specificAssignments2:
        return "+inf"
    return commonAssignments/(specificAssignments1 + specificAssignments2)
