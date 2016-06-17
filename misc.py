from __future__ import division
import numpy as np

from writeOnFiles import writeFile
from time import time

#Gets sample ID 
def getSampleIDList(samplesList):
    sampleIDList = []
    for sample in samplesList:
        if not mem(sample[0],sampleIDList):
            sampleIDList.append(sample[0])
    return sampleIDList

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

def containsSpecieCutPath(path,name,rank,includeEnd):
    i = 0
    n = len(path)
    while i < n and not (path[i] == (name,rank)):
        i += 1
    if i == n:
        return False,[]
    else:
        if includeEnd:
            return True,path[:i+1]
        else:
            return True,path[:i]

#@paths is the paths list of a TaxoTree
#The queries in @addNode (see TaxoTree) start with S-ranked nodes, then G-ranked nodes, and so on.
#@n = len(@paths)
def selectPath(paths,name,rank,n,includeEnd=False):
    i = 0
    while i < n:
        boolean,path = containsSpecieCutPath(paths[i],name,rank,includeEnd)
        if boolean:
            return path
        else:
            i += 1
    print "BUG: No path for %s,%s"%(name,rank)
    return []

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

def sumOp(x,y):
    if (x == "N") or (y == "N"):
        return "N"
    else:
        return (x+y)

#@xArray and @yArray have the same length and are tuples lists
def arraySum(xArray,yArray):
    return [(xArray[i][0],xArray[i][1],sumOp(yArray[i][1],xArray[i][2])) for i in range(len(xArray))]
            

#Computes the set of number of assignments to a certain bacteria (name,rank)
#The set contains (sampleID,numberAssignments) pairs
def getValueOneBacteria(samplesOccList,speciesList,name,rank):
    i = 0
    n = len(speciesList)
    while i < n and not (speciesList[i][0] == name and speciesList[i][1] == rank):
        i += 1
    values = []
    #len(sample) == len(speciesList)
    for sample in samplesOccList:
        values.append((sample[0],sample[i]))
    return values
    
#@bacteriaList is a list of (name,rank) pairs
def getValueBacteria(samplesOccList,speciesList,bacteriaList):
    result = [(sample[0],bacteria,0) for bacteria in bacteriaList for sample in samplesOccList]
    for bacteria in bacteriaList:
        #print bacteria[0],bacteria[1]
        r = getValueOneBacteria(samplesOccList,speciesList,bacteria[0],bacteria[1])
        result = arraySum(result,r)
    #print result
    answer = raw_input("Write Bacteria file? [Choose for formatting (string,(string,string),integer)]  Y/N\n")
    if (answer == "Y"):
        writeFile(result,"array")
    return result

#Computes the set of number of assignments to a sample verifying some requirements for a given metadatum
#See percentage to understand the role of defaultValue and intervalLength
def getValueOneMetadatumSelection(samplesInfoList,infoList,metadatum,defaultValue,intervalLength):
    i = 0
    n = len(infoList)
    while i < n and not (metadatum == infoList[i]):
        i += 1
    values = []
    for sample in samplesInfoList:
        if (sample[i] == defaultValue) or (sample[i] >= defaultValue - intervalLength/2 and sample[i] <= defaultValue + intervalLength/2):
            if not (sample[i] == "N"):
                values.append((sample[0],int(sample[i])))
            else:
                values.append((sample[0],sample[i]))
    return values

def getValueOneMetadatum(samplesInfoList,infoList,metadatum):
    i = 0
    n = len(infoList)
    while i < n and not (metadatum == infoList[i]):
        i += 1
    values = []
    for sample in samplesInfoList:
        if not (sample[i] == "N"):
            values.append((sample[0],int(sample[i])))
        else:
            values.append((sample[0],sample[i]))
    return values

#Selection of samples provided a default value for the metadata
def getValueMetadataSelection(samplesInfoList,infoList,metadata,defaultValues,intervalLengths):
    result = [(sample[0],metadatum,0) for metadatum in metadata for sample in samplesInfoList]
    print result
    assert (len(metadata) == len(defaultValues) and len(metadata) == len(intervalLengths))
    n = len(metadata)
    for i in range(n):
        r = getValueOneMetadatumSelection(samplesInfoList,infoList,metadata[i],defaultValues[i],intervalLengths[i])
        print r
        result = arraySum(result,r)
    print result
    answer = raw_input("Write Metadata file? [Choose for formatting (string,string,integer)] Y/N\n")
    if (answer == "Y"):
        writeFile(result,"array")
    return result

#@metadataList is a list of metadata
def getValueMetadata(samplesInfoList,infoList,metadata):
    result = [(sample[0],metadatum,0) for metadatum in metadata for sample in samplesInfoList]
    print result
    for metadatum in metadata:
        r = getValueOneMetadatum(samplesInfoList,infoList,metadatum)
        print r
        result = arraySum(result,r)
    print result
    answer = raw_input("Write Metadata file? [Choose for formatting (string,string,integer)] Y/N\n")
    if (answer == "Y"):
        writeFile(result,"array")
    return result
