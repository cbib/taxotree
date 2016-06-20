from __future__ import division
import numpy as np
import re

from writeOnFiles import writeFile

integer = re.compile("[0-9]+")

#Returns a list containing elements of list1 that are not elements of list2
#If the elements in lists are tuples, it will sorted using lexigraphical order with suitable order for each projection
#Sorting allows to run in worst case complexity of O(nlog(n) + mlog(m)) where n,m are the length of the two lists
def trimList(list1,list2):
    lst1 = sorted(list1,key=lambda x:x)
    lst2 = sorted(list2,key=lambda x:x)
    trimList = []
    while lst1 and lst2:
        x1 = lst1.pop()
        x2 = lst2.pop()
        if (x1 < x2):
            lst1.append(x1)
        elif (x1 > x2):
            trimList.append(x1)
            lst2.append(x2)
    #At the end of the loop, at least one of the two lists is empty
    if lst1:
        #We do not care of the order at the end of the merge
        return trimList + lst1
    return trimList

#Merge the elements of both lists, deleting multiple occurrences
def mergeList(list1,list2):
    lst1 = sorted(list1,key=lambda x:x)
    lst2 = sorted(list1,key=lambda x:x)
    union = []
    while lst1 and lst2:
        x1 = lst1.pop()
        x2 = lst2.pop()
        if (x1 == x2):
            union.append(x1)
        elif (x1 < x2):
            union.append(x2)
            lst1.append(x1)
        else:
            union.append(x1)
            lst2.append(x2)
    #At the end of the loop, at least one of the two lists is empty
    if lst1:
        #We do not care of the order at the end of the merge
        return union + lst1
    else:
        #@ls2 may be empty
        return union + lst2

#Returns a list of (name,rank,sampleHitList) tuples associated with
#nodes of the taxonomic tree reduced to the sample sampleName
#and the number of assignments in this tree
#NB: We need to keep the whole sampleHitList (and not only the
#number associated with sampleName) in order to apply set operations (intersection, ...)
def inSample(element,sampleNameList):
    #element[1] (number associated to a sample) must be non-zero
    if not element[1]:
        print "\n/!\ ERROR: [BUG] [misc/inSample] Null element."
        raise ValueError
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

#Gets sample IDs
def getSampleIDList(samplesList):
    sampleIDList = []
    for sample in samplesList:
        if not mem(sample[0],sampleIDList):
            sampleIDList.append(sample[0])
    #Sorts sample IDs in alphabetical order
    return sorted(sampleIDList,key=lambda x:x)

#For sorting lists of nodes (name,rank) by decreasing order S > G > F > O > C > P > K > R. Returns 1 if rank1 => rank2, -1 if rank1 < rank2
def compare(rank1,rank2):
    if (rank2 == "S"):
        return 1
    elif (rank1 == "S"):
        return -1
    elif (rank2 == "G"):
        return 1
    elif (rank1 == "G"):
        return -1
    elif (rank2 == "F"):
        return 1
    elif (rank1 == "F"):
        return -1
    elif (rank2 == "O"):
        return 1
    elif (rank1 == "O"):
        return -1
    elif (rank2 == "C"):
        return 1
    elif (rank1 == "C"):
        return -1
    elif (rank2 == "P"):
        return 1
    elif (rank1 == "P"):
        return -1
    elif (rank2 == "K"):
        return 1
    elif (rank1 == "K"):
        return -1
    elif (rank2 == "R"):
        return 1
    elif (rank1 == "R"):
        return -1

def sanitize(name):
    ls = name.split(" ")
    if (len(ls) == 1):
        return ls[0]
    sName = ""
    sLs = []
    for l in ls:
        if not (l == "" or l == "(class)" or l == "\n" or l == "#"):
            sLs.append(l)
    for l in sLs[:-1]:
        sName = sName + l + " "
    sName = sName + sLs[-1]
    return sName.split("\n")[0]

#is member function
def mem(x,ls):
    n = len(ls)
    for i in range(n):
        if (x == ls[i]):
            return True
    return False

def containsSpecie(path,name,rank):
    for x in path:
        if (x[0] == name) and (x[1] == rank):
            return True
    return False

#@paths is the paths list of a TaxoTree
#@n = len(@paths)
def selectPath(paths,name,rank,n):
    i = 0
    while i < n and not containsSpecie(paths[i],name,rank):
	i += 1
    if (i == n):
        print "\n/!\ ERROR: [BUG] [misc/selectPath] No path for (%s,%s)."%(name,rank)
        raise ValueError
    else:
        path = []
        x = paths[i]
        #Memorizes only the part of the path that leads to (name,rank)
        n = len(x)
        i = 0
        while (i < n) and not (x[i][0] == name and x[i][1] == rank):
            path.append(x[i])
            i += 1
        return path

def isLeaf(paths,name,rank,allNodes):
    l = []
    for ls in paths:
        if containsSpecie(ls,name,rank):
            l.append(ls)
    if allNodes:
        #deletes path ending with this node
        l2 = []
        for ls in l:
            if not (ls[-1][0] == name and ls[-1][1] == rank):
                l2.append(ls)
        return (len(l2) == 0)
    else:
        #Searches which paths may lead to (name,rank)
        return (len(l) == 1)

def setOperations(paths,name1,rank1,name2,rank2,allNodes=False):
    path1 = selectPath(paths,name1,rank1)
    path2 = selectPath(paths,name2,rank2)
    n = min(len(path1),len(path2))
    #if there is more than one path to the nodes, or no path
    if (n < 1):
        print "\n/!\ ERROR: [BUG] [misc/setOperations] It is not a tree."
        raise ValueError
    else:
        commonPath = []
        i = 0
        #As long as path1 and path2 are not empty
        while (i < n) and (path1[i] == path2[i]):
            commonPath.append(path1[i])
            i += 1
        return commonPath,path1[i+1:],path2[i+1:]

#Computes LCA from the list paths of a TaxoTree
def taxoLCA(paths,name1,rank1,name2,rank2,allNodes=False):
    commonPath,_,_ = setOperations(paths,name1,rank1,name2,rank2,allNodes=False)
    return commonPath[-1]

#Computes the positions of number of assignments to each bacteria of the group in the occurrences matrix
def getPositionBacteria(bacteriaList,speciesList):
    i = 0
    positions = []
    n = len(speciesList)
    while i < n:
        if mem(speciesList[i],bacteriaList):
            positions.append(i)
        i += 1
    return positions

#Computes the set of number of assignments to a certain bacteria (name,rank)
#Returns a (sample,number of assignments to the group of bacterias in sample) pair list
def getValueBacteria(samplesOccList,speciesList,bacteriaPos,sampleList):
    resultList = []
    assignments = 0
    for sample in sampleList:
        i = 0
        n = len(samplesOccList)
        #gets the line associated to the sample in samplesOccList
        while i < n and not (samplesOccList[i][0] == sample):
            i += 1
        for pos in bacteriaPos:
            assignments += samplesOccList[i][pos]
        resultList.append((sample,assignments))
    return resultList
    
#Returns xArray, yArray where xArray and yArray are both arrays containing the number of assignments to a chosen group of bacterias depending on the sample considered
def getValueBacteriaBacteria(samplesOccList,speciesList,sampleIDList,bacterias1,bacterias2):
    xArray,yArray = [],[]
    #Stores the positions of number of assignments to each bacteria of the group in the occurrences matrix
    bacteriaPos1 = getPositionBacteria(bacterias1,speciesList)
    bacteriaPos2 = getPositionBacteria(bacterias2,speciesList)
    for sample in sampleIDList:
        #gets the number of assignments to the group of bacterias in the sample 
        xArray.append(getValueBacteria(samplesOccList,speciesList,bacteriaPos1,[sample])[0])
        yArray.append(getValueBacteria(samplesOccList,speciesList,bacteriaPos2,[sample])[0])
    answer = raw_input("Write both Bacteria files? Y/N\n")
    if (answer == "Y"):
        print "Saving first file"
        writeFile(xArray,"**** Values of assignments in all samples of nodes: " + str(bacteria1) + "\n\n","array")
        print "Saving second file"
        writeFile(yArray,"**** Values of assignments in all samples of nodes: " + str(bacteria2) + "\n\n","array")
    return xArray,yArray

#Returns xArray and yArray, where yArray contains the values of selected metadatum and xArray contains the number of assignments to a chosen group of bacterias depending on the value of the metadatum 
def getValueBacteriaMetadata(samplesInfoList,infoList,bacterias,metadatum):
    xArray = []
    #Stores the positions of number of assignments to each bacteria of the group in the occurrences matrix
    bacteriaPos = getPositionBacteria(bacterias,speciesList)
    #computes the number of colum which matches the metadatum in infoList
    i = 0
    n = len(infoList)
    while i < n and not (infoList[i] == metadatum):
        i += 1
    #Getting the set of values of the metadatum
    #Sorting samples according to the values of the metadatum
    sampleSorted = sorted(samplesInfoList,key=lambda x: x[i])
    #List of list of samples: one sublist matches a value of the metadatum
    valueSampleMetadatum = []
    #The set of values of the metadatum
    valueSet = []
    if not len(sampleSorted):
        print "\n/!\ ERROR: You have selected no sample."
        raise ValueError
    sample = sampleSorted.pop()
    while not integer.match(sample[i]):
        sample = sampleSorted.pop()
    currValue = sample[i]
    valueSet.append(currValue)
    while sampleSorted:
        valueSample = []
        while sample[i] == currValue:
            valueSample.append(sample)
            sample = sampleSorted.pop()
            while not integer.match(sample[i]):
                sample = sampleSorted.pop()
        valueSampleMetadatum.append(valueSample)
        currValue = sample[i]
        valueSet.append(currValue)
    #Integer values of metadatum are sorted
    yArray = sorted(valueSet,key=lambda x:x)
    for sampleValueList in valueSampleMetadatum:
        xArray.append(getValueBacteria(samplesOccList,speciesList,bacteriaPos,sampleValueList))
    answer = raw_input("Write both Bacteria and Metadatum files? Y/N\n")
    if (answer == "Y"):
        print "Saving first file"
        writeFile(xArray,"**** Values of assignments in all samples of nodes: " + str(bacteria) + "\n\n","array")
        print "Saving second file"
        writeFile(yArray,"**** Values of metadatum: " + str(metadatum) + "\n\n","array")
    return xArray,yArray
