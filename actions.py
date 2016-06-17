import numpy as np

from writeOnFiles import writeFile
from parsingMatrix import parseMatrix
from parsingInfo import parseInfo
from parsingTree import parseTree
from taxoTree import TaxoTree
from misc import getValueBacteria,getValueMetadata,getValueMetadataSelection,mem

from totalratio import compute,countAssignmentsInCommon,countAssignments,totalRatio,totalRatioNormalized,diffRatio,diffRatioNormalized
from patternRatio import patternRatio,patternRatioNormalized
from pearsonCorrelation import samplePearson,populationPearson,printProbabilityLawsList
from percentage import percentageAssign
from similarityCoefficient import similarity

#@dataArray = [samplesInfoList,infoList,samplesOccList,speciesList,paths,n,nodesList,taxoTree,sampleIDList]

#Parsing functions
def parseList(string):
    if not (len(string.split(",")) == 1):
        raise ValueError
    elif not (len(string.split(":")) == 1):
        raise ValueError
    return string.split(";")

def parseListNode(string):
    ls = string.split(";")
    res = []
    for node in ls:
        name,rank = node.split(",")[0].split("(")[-1],node.split(",")[-1].split(")")[0]
        res.append((name,rank))
    return res

def parseIntList(string):
    if not (len(string.split(",")) == 1):
        raise ValueError
    elif not (len(string.split(":")) == 1):
        raise ValueError
    l = string.split(";")
    resultList = []
    for s in l:
        if not (s == "+inf" or s == "-inf"):
            resultList.append(int(s))
        else:
            resultList.append(s)
    return resultList

#___________________________________________________________________________

#Macros
def listNodes(nodeList):
    string = ""
    for l in nodeList[:-1]:
        string += str(l) + ", "
    string += str(nodeList[-1])
    return string

def listMetadata(metadataList,interval1List,interval2List):
    string = ""
    n = len(metadataList)
    for i in range(n-1):
        if (interval1List[i] == interval2List[i]):
            string += metadataList[i] + " (value equal to " + str(interval1List[i]) + "), "
        else:
            string += metadataList[i] + " (value between " + str(interval1List[i]) + " and " + str(interval2List[i]) + "), "
    if (interval1List[-1] == interval2List[-1]):
        string += metadataList[-1] + " (value equal to " + str(interval1List[-1]) + ")"
    else:
        string += metadataList[-1] + " (value between " + str(interval1List[-1]) + " and " + str(interval2List[-1]) + ")"
    return string

def createSampleNameList(dataArray):
    sampleIDList = dataArray[8]
    sampleType = parseList(raw_input("Do you want to select one by one samples or to select samples matching requirements on metadata? Y/N \n"))
    if (sampleType == "Y"):
        print sampleIDList
        sampleNameList = parseList(raw_input("Input the first list of samples using the ID printed above. [ e.g. OPNA-J90;DUGA-J0;GATE-J0 ]\n"))
        isInDatabase(sampleNameList,sampleIDList)
    elif (sampleType == "N"):
        print dataArray[1]
        metadataList = parseList(raw_input("Input the list of metadata you want to consider among those written above. [ e.g. ATB_IV;ATB_Os ]\n"))
        isInDatabase(metadataList,dataArray[1])
        interval1List = parseIntList(raw_input("Input the list of lower interval bounds corresponding to metadata above. [ Please refer to README for more details. e.g. 1;2;1;0;0 ]\n"))
        assert (len(interval1List) == len(metadataList))
        interval2List = parseIntList(raw_input("Input the list of upper interval bounds corresponding to metadata above. [ Please refer to README for more details. e.g. 3;2;1;2;0 ]\n"))
        assert (len(interval2List) == len(metadataList))
        sampleNameList = getValueMetadataSelection(dataArray[0],dataArray[1],metadataList,interval1List,interval2List)
    else:
        raise ValueError
    return sampleNameList

#Checks if the elements in @parselist belong to @datalist else returns an error
def isInDatabase(parseList,dataList):
    for pl in parseList:
        if not mem(pl,dataList):
            raise ValueError

#____________________________________________________________________________

#Actions
def totalDiffRatioAct(dataArray):
    sampleIDList = dataArray[8]
    print "[First set of samples]"
    #sampleNameList1 = createSampleNameList(dataArray)
    print sampleIDList
    sampleNameList1 = parseList(raw_input("Input the first list of samples using the ID printed above. [e.g. OPNA-J90;DUGA-J0;GATE-J0 ]\n"))
    isInDatabase(sampleNameList1,sampleIDList)
    print "[Second set of samples]"
    #sampleNameList2 = createSampleNameList(dataArray)
    sampleNameList2 = parseList(raw_input("Input the second list of samples using the ID printed above. [e.g. OPNA-J90;DUGA-J0;GATE-J0 ]\n"))
    isInDatabase(sampleNameList2,sampleIDList)
    common,in1,in2,_,_,_,_,_ = compute(dataArray[7],sampleNameList1,sampleNameList2)
    commonA = countAssignmentsInCommon(common,sampleNameList1,sampleNameList2)
    numberA1 = countAssignments(in1,sampleNameList1)
    numberA2 = countAssignments(in2,sampleNameList2)
    tratio = totalRatio(commonA,numberA1,numberA2)
    ntRatio = totalRatioNormalized(commonA,numberA1,numberA2)
    dratio = diffRatio(commonA)
    ndRatio = diffRatioNormalized(commonA,numberA1,numberA2)
    print "Total Ratio Distance is: " + str(tratio)
    print "normalized Total Ratio is: " + str(ntRatio) + "\n[The more it is close to 1, the more the two groups are alike]"
    print "Diff Ratio Distance is: " + str(dratio)
    print "normalized Diff Ratio is: " + str(ndRatio) + "\n[The more it is close to 0, the more the two groups are alike]"
    answer = raw_input("Save the results? Y/N\n")
    if (answer == "Y"):
        data = "Total Ratio Results ****\n for lists " + str(sampleNameList1) + "\n and " + str(sampleNameList2) + "\n\nTotal Ratio Distance is: " + str(tratio) + "\n normalized Total Ratio is: " + str(ntRatio) + "\nDiff Ratio Distance is: " + str(dratio) + "\n normalized Diff Ratio is: " + str(ndRatio) +"\n\n END ****"  
        writeFile(data,"","text")

#____________________________________________________________________________

def patternRatioAct(dataArray):
    #Prints the samples ID
    sampleIDList = dataArray[8]
    print sampleIDList
    sampleName1 = raw_input("Input the first sample name using the ID printed above. [e.g OPNA-J90 ]\n")
    isInDatabase([sampleName1],sampleIDList)
    sampleName2 = raw_input("Input the second sample name using the ID printed above. [e.g. DUGA-J0 ]\n")
    isInDatabase([sampleName2],sampleIDList)
    ratio,patternList = patternRatio(dataArray[7],sampleName1,sampleName2)
    nRatio,_ = patternRatioNormalized(dataArray[7],sampleName1,sampleName2)
    print "(non normalized) Pattern Ratio is: ",ratio
    print "normalized Pattern Ratio is: ",nRatio
    print "Pattern is:"
    print patternList
    answer = raw_input("Save the results? Y/N\n")
    if (answer == "Y"):
        data = "Pattern Ratio Results ****\n for samples " + sampleName1 + "\n and " + sampleName2 + "\n\n (non normalized) Pattern Ratio is: " + str(ratio) + "\n normalized Pattern Ratio is: " + str(nRatio) + "\n Pattern is: " + str(patternList) + "\n\n END ****"  
        writeFile(data,"","text")

#____________________________________________________________________________
        
#@dataArray = [samplesInfoList,infoList,samplesOccList,speciesList,paths,n,nodesList,taxoTree]
def percentageAct(dataArray):
    uTree = raw_input("Do you to get percentage of assignments to subtrees or to bacterias themselves? subtree/bacteria \n")
    usingTree = (uTree == "subtree")
    nodesGroup = parseListNode(raw_input("Input the list of nodes/roots of subtrees you want to consider. [Please look at the taxonomic tree file to help you: e.g. (Clostridium,G);(Elusimicrobia,C);(Bacteria,K) ]\n"))
    isInDatabase(nodesGroup,dataArray[6])
    print dataArray[1]
    metadataList = parseList(raw_input("Input the list of metadata you want to consider among those written above. [e.g. ATB_IV;ATB_Os ]\n"))
    isInDatabase(metadataList,dataArray[1])
    interval1List = parseIntList(raw_input("Input the list of lower interval bounds corresponding to metadata above. [Please refer to README for more details. e.g. 1;2;1;0;0 ]\n"))
    assert (len(interval1List) == len(metadataList))
    interval2List = parseIntList(raw_input("Input the list of upper interval bounds corresponding to metadata above. [Please refer to README for more details. e.g. 3;2;1;2;0 ]\n"))
    assert (len(interval2List) == len(metadataList))
    result = percentageAssign(dataArray[0],dataArray[1],metadataList,interval1List,interval2List,dataArray[7],nodesGroup,dataArray[2],dataArray[3],usingTree)
    print "[Preview.]"
    print result
    l = len(result)
    data = np.zeros(l)
    for i in range(l):
        data[i] = result[i]
    answer = raw_input("Save the results? Y/N\n")
    if (answer == "Y"):  
        writeFile(data,"Percentage of assignments in the group of nodes: " + listNodes(nodesGroup) + "\ndepending on metadata (for each line): " + listMetadata(metadataList,interval1List,interval2List),"array")

#_____________________________________________________________________________
        
def creatingArray(typeArray,valueArrayList,samplesOccList,speciesList,samplesInfoList,infoList):
    if (typeArray == "bacteria"):
        #then valueArrayList = (name,rank) list
        return getValueBacteria(samplesOccList,speciesList,valueArrayList)
    #typeArray == metadatum
    else:
        return getValueMetadata(samplesInfoList,infoList,valueArrayList)

#_____________________________________________________________________________

#@dataArray = [samplesInfoList,infoList,samplesOccList,speciesList,paths,n,nodesList,taxoTree]
def pearsonAct(dataArray):
    pearsonType = raw_input("Do you want to compute sample Pearson coefficient or population Pearson coefficient? sample/population \n")
    if (pearsonType == "sample"):
        print "\nFirst set of values\n"
        type1 = raw_input("Value are of type bacteria or metadatum? bacteria/metadatum\n")
        if (type1 == "bacteria"):
            value1 = parseListNode(raw_input("What value? [(name,rank) list for bacteria, e.g. (Clostridium,G);(Elusimicrobia,C) ]\n"))
            isInDatabase(value1,dataArray[6])
        elif (type1 == "metadatum"):
            value1 = parseList(raw_input("What value? [metadatum list for metadatum, e.g. ATB_IV;ATB_Os ]\n"))
            isInDatabase(value1,dataArray[1])
        else:
            raise ValueError
        xArray = creatingArray(type1,value1,dataArray[2],dataArray[3],dataArray[0],dataArray[1])
        print "\nSecond set of values\n"
        type2 = raw_input("Value are of type bacteria or metadatum? bacteria/metadatum\n")
        if (type2 == "bacteria"):
            value2 = parseListNode(raw_input("What value? [(name,rank) list for bacteria, e.g. (Clostridium,G);(Elusimicrobia,C) ]\n"))
            isInDatabase(value2,dataArray[6])
        elif (type2 == "metadatum"):
            value2 = parseList(raw_input("What value? [metadatum list for metadatum, e.g. ATB_IV;ATB_Os ]\n"))
            isInDatabase(value2,dataArray[1])
        else:
            raise ValueError
        yArray = creatingArray(type2,value2,dataArray[2],dataArray[3],dataArray[0],dataArray[1])
        pearson = samplePearson(xArray,yArray)
        print "Pearson Sample coefficient is: " + str(pearson)
    else:
        print "First set of values\n"
        type1 = raw_input("Value are of type bacteria or metadatum? bacteria/metadatum\n")
        if (type1 == "bacteria"):
            value1 = parseListNode(raw_input("What value? [(name,rank) list for bacteria, e.g. (Clostridium,G);(Elusimicrobia,C) ]\n"))
            isInDatabase(value1,dataArray[6])
        elif (type1 == "metadatum"):
            value1 = parseList(raw_input("What value? [metadatum list for metadatum, e.g. ATB_IV;ATB_Os ]\n"))
            isInDatabase(value1,dataArray[1])
        else:
            raise ValueError
        xArray = creatingArray(type1,value1,dataArray[2],dataArray[3],dataArray[0],dataArray[1])
        print "Second set of values\n"
        type2 = raw_input("Value are of type bacteria or metadatum? bacteria/metadatum\n")
        if (type2 == "bacteria"):
            value2 = parseListNode(raw_input("What value? [(name,rank) list for bacteria, e.g. (Clostridium,G);(Elusimicrobia,C) ]\n"))
            isInDatabase(value2,dataArray[6])
        elif (type2 == "metadatum"):
            value2 = parseList(raw_input("What value? [metadatum list for metadatum, e.g. ATB_IV;ATB_Os ]\n"))
            isInDatabase(value2,dataArray[1])
        else:
            raise ValueError
        yArray = creatingArray(type2,value2,dataArray[2],dataArray[3],dataArray[0],dataArray[1])
        printProbabilityLawsList()
        probList = ["UniformProbability","UniformProbabilityProduct"]
        p1 = raw_input("Enter the law of probability for first values [among the ones above]\n")
        isInDatabase(p1,probList)
        p2 = raw_input("Enter the law of probability for second values [among the ones above]\n")
        isInDatabase(p2,probList)
        p3 = raw_input("Enter the law of probability for the product of first values with second values [among the ones above]\n")
        isInDatabase(p3,probList)
        pearson = populationPearson(xArray,yArray,p1,p2,p3)
        print "Pearson Population coefficient is: " + str(pearson) 
    answer = raw_input("Save the results? Y/N\n")
    if (answer == "Y"):
        data = "The " + pearsonType + " Pearson coefficient for values: **** \n\n" + str(xArray) + "\n and " + str(yArray) + "\n is : " + str(pearson) + "\n\n END ****"
        writeFile(data,"","text")

#_____________________________________________________________________________

def similarityAct(dataArray,iMatrix):
    print "/!\ Computing similarity matrix..."
    m = similarity(dataArray[0],dataArray[1])
    print "[Preview.]"
    print m
    answer = raw_input("Save the results? Y/N\n")
    if (answer == "Y"):  
        writeFile(m,"Similarity coefficients between patients for file meta/" + iMatrix + ".csv:\n" + listNodes(dataArray[8]),"array")
    return m
