from __future__ import division
import numpy as np
import re

from writeOnFiles import writeFile
from parsingMatrix import parseMatrix
from parsingInfo import parseInfo
from parsingTree import parseTree
from taxoTree import TaxoTree,printTree
from misc import getValueBacteriaBacteria,getValueBacteriaMetadata,mem

from totalratio import compute,countAssignmentsInCommon,countAssignments,totalRatio,totalRatioNormalized,diffRatio,diffRatioNormalized
from patternRatio import patternRatio,enumerateCommonPatterns,enumerateSpecificPatterns
from pearsonCorrelation import samplePearson,populationPearson,printProbabilityLawsList
from percentage import percentageAssign,computeSamplesInGroup
from similarityCoefficient import similarity

#@dataArray = [samplesInfoList,infoList,samplesOccList,speciesList,paths,n,nodesList,taxoTree,sampleIDList,#similarityMatrix]

integer = re.compile("[0-9]+")

#Parsing functions
def parseList(string):
    if not (len(string.split(",")) == 1):
        print "\n/!\ ERROR: Do not use ',' as a separator: rather use ';'."
        raise ValueError
    elif not (len(string.split(":")) == 1):
        print "\n/!\ ERROR: Do not use ':' as a separator: rather use ';'."
        raise ValueError
    return string.split(";")

def parseListNode(string):
    if not (len(string.split(":")) == 1):
        print "\n/!\ ERROR: Do not use ':' as a separator: rather use ';'."
        raise ValueError
    ls = string.split(";")
    res = []
    for node in ls:
        nodeSplit = node.split(",")
        if not (len(nodeSplit) == 2):
            print "\n/!\ ERROR: Please use ',' as a separator for name,rank of a bacteria."
            raise ValueError
        nodeSplitName = nodeSplit[0].split("(")
        if not (len(nodeSplitName) == 2):
            print "\n/!\ ERROR: Please use the syntax '([name],[rank])' for each bacteria."
            raise ValueError
        nodeSplitRank = nodeSplit[-1].split(")")
        if not (len(nodeSplitRank) == 2):
            print "\n/!\ ERROR: Please use the syntax '([name],[rank])' for each bacteria."
            raise ValueError
        name,rank = nodeSplitName[-1],nodeSplitRank[0]
        res.append((name,rank))
    return res

def parseIntList(string):
    if not (len(string.split(",")) == 1):
        print "\n/!\ ERROR: Do not use ',' as a separator: rather use ';'."
        raise ValueError
    elif not (len(string.split(":")) == 1):
        print "\n/!\ ERROR: Do not use ':' as a separator: rather use ';'."
        raise ValueError
    l = string.split(";")
    resultList = []
    for s in l:
        if integer.match(s):
            resultList.append(int(s))
        elif s == "+inf" or s == "-inf":
            resultList.append(s)
        else:
            print "\n/!\ ERROR: Here you can only use integers or '+inf' or '-inf'."
            raise ValueError
    return resultList

#___________________________________________________________________________

#Macros for formatting
#Printing pretty lists of nodes
def listNodes(nodeList):
    string = ""
    for l in nodeList[:-1]:
        string += str(l) + ", "
    string += str(nodeList[-1])
    return string

#@stringNode is assumed to be a (name,rank) pair, with name and rank being strings
#@sanitizeNode allows it to be printed "(name,rank)" and not "('name','rank')"
def sanitizeNode(stringNode):
    return "(" + stringNode[0] + "," + stringNode[1] + ")"

#Printing pretty lists of metadata with their default values
def listSampleInvolved(metadataList,interval1List,interval2List,sampleNameList):
    string = ""
    if not metadataList and not interval1List and not interval2List and not sampleNameList:
      print "\n/!\ ERROR: You have selected no sample."
      raise ValueError
    #If samples were selected one by one
    elif sampleNameList:
        string += "\ndepending on the group of samples: "
        for sl in sampleNameList[:-1]:
            string += str(sl) + ", "
        string += str(sampleNameList[-1])
    #If samples were selected according to metadata values (len(metadataList) = len(interval1List) = len(interval2List))
    if metadataList:
        string += "\nselected on metadata (for each line): "
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

#Selecting samples in two ways: either choose each of them one by one, or selecting according to default values of certain metadatum
def createSampleNameList(dataArray,percentage=False):
    metadataList = []
    interval1List = []
    interval2List = []
    sampleIDList = dataArray[8]
    if percentage:
        i = raw_input("/!\ How many different lists of samples do you want?\n")
        if not integer.match(i):
            print "\n/!\ ERROR: You need to enter a integer here!"
            raise ValueError
        numberList = int(i)
        sampleNameList = []
        if (numberList < 1):
            print "\n/!\ ERROR: Empty set of lists of samples!"
            raise ValueError
        while numberList:
            answer = raw_input("Do you want to select samples one by one, or to select samples matching requirements on metadata? one/matching \n")
            if (answer == "one"):
                if (len(sampleIDList) < 2):
                    print "\n/!\ ERROR: List of samples is empty or only of length one!..."
                    raise ValueError
                print sampleIDList
                sampleNameList11 = parseList(raw_input("Input the list of samples using the ID printed above. [e.g. " + sampleIDList[0] + ";"+ sampleIDList[1] + " ]\n"))
            elif (answer == "matching"):
                print dataArray[1]
                metadataList = parseList(raw_input("Input the list of metadata you want to consider among those written above. [ e.g. " + dataArray[1][0] + ";" + dataArray[1][-1] + " ]\n"))
                isInDatabase(metadataList,dataArray[1])
                interval1List = parseIntList(raw_input("Input the list of lower interval bounds corresponding to metadatum/metadata above. [ Please refer to README for more details. e.g. 1;2 ]\n"))
                if not (len(interval1List) == len(metadataList)):
                    print "\n/!\ ERROR: You need to enter the same number of lower bounds than of metadata!"
                    raise ValueError
                interval2List = parseIntList(raw_input("Input the list of upper interval bounds corresponding to metadatum/metadata above. [ Please refer to README for more details. e.g. 3;2 ]\n"))
                if not (len(interval2List) == len(metadataList)):
                    print "\n/!\ ERROR: You need to enter the same number of upper bounds than of metadata!"
                    raise ValueError
                sampleNameList11 = computeSamplesInGroup(dataArray[0],dataArray[1],metadataList,interval1List,interval2List)[0]
            else:
                print "\n/!\ ERROR: You need to answer either 'one' or 'matching' and not: \"",answer,"\"."
                raise ValueError
            isInDatabase(sampleNameList11,sampleIDList)
            sampleNameList.append(sampleNameList11)
            numberList -= 1
    else:
            answer = raw_input("Do you want to select samples one by one, or to select samples matching requirements on metadata? one/matching \n")
            if (answer == "one"):
                if (len(sampleIDList) < 2):
                    print "\n/!\ ERROR: List of samples is empty or only of length one!..."
                    raise ValueError
                print sampleIDList
                sampleNameList = parseList(raw_input("Input the list of samples using the ID printed above. [e.g. " + sampleIDList[0] + ";"+ sampleIDList[1] + " ]\n"))
            elif (answer == "matching"):
                print dataArray[1]
                metadataList = parseList(raw_input("Input the list of metadata you want to consider among those written above. [ e.g. " + dataArray[1][0] + ";" + dataArray[1][-1] + " ]\n"))
                isInDatabase(metadataList,dataArray[1])
                interval1List = parseIntList(raw_input("Input the list of lower interval bounds corresponding to metadatum/metadata above. [ Please refer to README for more details. e.g. 1;-inf;3 ]\n"))
                if not (len(interval1List) == len(metadataList)):
                    print "\n/!\ ERROR: You need to enter the same number of lower bounds than of metadata!"
                    raise ValueError
                interval2List = parseIntList(raw_input("Input the list of upper interval bounds corresponding to metadatum/metadata above. [ Please refer to README for more details. e.g. +inf;2;1 ]\n"))
                if not (len(interval2List) == len(metadataList)):
                    print "\n/!\ ERROR: You need to enter the same number of upper bounds than of metadata!"
                    raise ValueError
                sampleNameList = computeSamplesInGroup(dataArray[0],dataArray[1],metadataList,interval1List,interval2List)[0]
            else:
                print "\n/!\ ERROR: You need to answer either 'one' or 'matching' and not: \"",answer,"\"."
                raise ValueError
            isInDatabase(sampleNameList,sampleIDList)
    return sampleNameList,metadataList,interval1List,interval2List

#Checks if the elements in @parselist belong to @datalist else returns an error
def isInDatabase(parseList,dataList):
    for pl in parseList:
        if not mem(pl,dataList):
            n = len(dataList)
            if not n:
                print "\n/!\ ERROR: [BUG] [actions/isInDatabase] Empty list."
            else:
                print "\n/!\ ERROR: '" + str(pl) + "' is not in the database beginning with: " + str(dataList[:min(n-1,3)]) + "."
            raise ValueError

#____________________________________________________________________________

#Actions
def totalDiffRatioAct(dataArray):
    print "First list of samples."
    sampleNameList1,metadataList1,interval1List1,interval2List1 = createSampleNameList(dataArray)
    print "Second list of samples."
    sampleNameList2,metadataList2,interval1List2,interval2List2 = createSampleNameList(dataArray)
    common,in1,in2,_,_,_,_,_ = compute(dataArray[7],sampleNameList1,sampleNameList2)
    commonA = countAssignmentsInCommon(common,sampleNameList1,sampleNameList2)
    numberA1 = countAssignments(in1,sampleNameList1)
    numberA2 = countAssignments(in2,sampleNameList2)
    tratio = totalRatio(commonA,numberA1,numberA2)
    ntRatio = totalRatioNormalized(commonA,numberA1,numberA2)
    dratio = diffRatio(commonA)
    ndRatio = diffRatioNormalized(commonA,numberA1,numberA2)
    print "\nTotal Ratio Distance is: " + str(tratio)
    print "normalized Total Ratio is: " + str(ntRatio) + "\n[The more it is close to 1, the more the two groups are alike]\n"
    print "Diff Ratio Distance is: " + str(dratio)
    print "normalized Diff Ratio is: " + str(ndRatio) + "\n[The more it is close to 0, the more the two groups are alike]\n"
    print "[If you have obtained +inf (resp. -inf), it could mean you have selected no sample.]\n"
    answer = raw_input("Save the results? Y/N\n")
    if (answer == "Y"):
        data = "Total Ratio Results ****\n for lists " + str(sampleNameList1) + "\n"
        if metadataList1:
            data += "selected on metadata: " + str(metadataList1) + "with extreme values: " + str(interval1List1) + " (lower bounds) and " + str(interval2List1) + " (upper bounds) \n"
        data += " and " + str(sampleNameList2) + "\n"
        if metadataList2:
            data += "selected on metadata: " + str(metadataList2) + "with extreme values: " + str(interval1List2) + " (lower bounds) and " + str(interval2List2) + " (upper bounds) \n"
        data += "\nTotal Ratio Distance is: " + str(tratio) + "\n normalized Total Ratio is: " + str(ntRatio) + "\nDiff Ratio Distance is: " + str(dratio) + "\n normalized Diff Ratio is: " + str(ndRatio) +"\n\nEND OF FILE ****"  
        writeFile(data,"","text")

#____________________________________________________________________________

def patternRatioAct(dataArray):
    print "First list of samples."
    sampleNameList1,metadata1,interval11,interval21 = createSampleNameList(dataArray)
    print "Second list of samples."
    sampleNameList2,metadata2,interval12,interval22 = createSampleNameList(dataArray)
    commonPatternsList = enumerateCommonPatterns(dataArray[7],sampleNameList1,sampleNameList2)
    specificPatternsList1 = enumerateSpecificPatterns(dataArray[7],sampleNameList1,sampleNameList2)
    specificPatternsList2 = enumerateSpecificPatterns(dataArray[7],sampleNameList2,sampleNameList1)
    pRatio = patternRatio(commonPatternsList,specificPatternsList1,specificPatternsList2)
    #Only printing patterns of length > 1
    print "\n--- Total number of common patterns: ",len(commonPatternsList)
    print "--- Common patterns of length > 1 ---"
    if commonPatternsList:
        for x in commonPatternsList:
            if len(x[0]) > 1:
                print x[0]
    else:
        print "No pattern of length > 1."
    print "\n--- Total number of specific patterns in",sampleNameList1
    if metadata1:
        print "selected on metadata: ",str(metadata1),"with lower and upper bounds being",str(interval11),"and",str(interval21),":"
    print len(specificPatternsList1)
    print "--- Specific patterns of length > 1 in",sampleNameList1,"---"
    if specificPatternsList1:
        for x in specificPatternsList1:
            if len(x[0]) > 1:
                print x[0]
    else:
        print "No pattern of length > 1."
    print "\n--- Total number of specific patterns in",sampleNameList2
    if metadata2:
        print "selected on metadata: ",str(metadata2),"with lower and upper bounds being",str(interval12),"and",str(interval22),":"
    print len(specificPatternsList2)
    print "--- Specific patterns of length > 1 in",sampleNameList2,"---"
    if specificPatternsList2:
        for x in specificPatternsList2:
            if len(x[0]) > 1:
                print x[0]
    else:
        print "No pattern of length > 1."
    print "\nPattern Ratio is: ",pRatio,"\n"
    print "[ If pattern ratio is superior to one, it means the two groups of samples are quite alike. Please read README ]"
    print "[ If you obtained +inf, if there are common patterns (of length 1 or superior to 1), it could mean both groups of samples contain exactly the same set of nodes. If there is no common pattern, it could mean there is no sample in both groups ]\n"
    answer = raw_input("Save the results? Y/N\n")
    if (answer == "Y"):
        data = "Pattern Ratio Results ****\nfor lists of samples " + str(sampleNameList1) + "\n"
        if metadata1:
            data += "selected on metadata: " + str(metadata1) + " with lower and upper bounds being " + str(interval11) + " and " + str(interval21) + "\n"
        data += "\nand " + str(sampleNameList2) + "\n"
        if metadata2:
            data += "selected on metadata: " + str(metadata2) + " with lower and upper bounds being " + str(interval12) + " and " + str(interval22) + "\n"
        data += "\n-> Pattern Ratio is: " + str(pRatio) + "\n\nPrinting patterns: first is the list of nodes in the pattern, then the total number of assignations in this pattern and eventually the total number of nodes in the pattern\n\n-> Common Patterns:\n"
        for x in commonPatternsList:
            data += str(x) + "\n"
        data += "\n-> Specific patterns to " + str(sampleNameList1) + ":\n"
        for x in specificPatternsList1:
            data += str(x) + "\n"
        data += "\n-> Specific patterns to " + str(sampleNameList2) + ":\n"
        for x in specificPatternsList2:
            data += str(x) + "\n"
        data += "\nEND OF FILE ****"
        writeFile(data,"","text")

#____________________________________________________________________________
        
def percentageAct(dataArray):
    uTree = raw_input("Do you to get percentage of assignments to subtrees or to bacterias themselves? subtree/bacteria \n")
    usingTree = (uTree == "subtree")
    if not (uTree == "subtree" or uTree == "bacteria"):
        print "\n/!\ ERROR: You need to answer 'bacteria' or 'subtree'."
        raise ValueError
    nodesGroup = parseListNode(raw_input("Input the list of nodes/roots of subtrees you want to consider. [ Please look at the taxonomic tree file to help you: e.g. " + sanitizeNode(dataArray[6][-3]) + ";" + sanitizeNode(dataArray[6][1]) + ";" + sanitizeNode(dataArray[6][-1]) + ". ]\n"))
    isInDatabase(nodesGroup,dataArray[6])
    sampleNameList,metadataList,interval1List,interval2List = createSampleNameList(dataArray,True)
    result = percentageAssign(dataArray[0],dataArray[1],sampleNameList,dataArray[7],nodesGroup,dataArray[2],dataArray[3],usingTree)
    print "\n[Preview.]"
    print result
    l = len(result)
    data = np.zeros(l)
    for i in range(l):
        data[i] = result[i]
    print ""
    answer = raw_input("Save the results? Y/N\n")
    if (answer == "Y"):  
        writeFile(data,"Percentage of assignments ****\nin the group of nodes: " + listNodes(nodesGroup) + listSampleInvolved(metadataList,interval1List,interval2List,sampleNameList),"array")

#_____________________________________________________________________________


#@dataArray = [samplesInfoList,infoList,samplesOccList,speciesList,paths,n,nodesList,taxoTree,sampleIDList,#similarityMatrix]
#Returns two arrays xArray and yArray, where yArray gives the value of a certain quantity depending on the values of xArray
def creatingArray(dataArray,pearson=False):
    #Available cases in Pearson function
    if pearson:
        typeInput = raw_input("Do you want to compute bacteria/bacteria or bacteria/metadatum? BB/BM [ Please read README for details. ]\n")
        if (typeInput == "BB"):
            valueInput1 = parseListNode(raw_input("Choose the first group of bacterias [ Read the taxonomic tree to help you: e.g. " + sanitizeNode(dataArray[6][-3]) + ";" + sanitizeNode(dataArray[6][1]) + ";" + sanitizeNode(dataArray[6][-1]) + " ]\n"))
            isInDatabase(valueInput1,dataArray[6])
            valueInput2 = parseListNode(raw_input("Choose the second group of bacterias [ Read the taxonomic tree to help you: e.g. " + sanitizeNode(dataArray[6][-3]) + ";" + sanitizeNode(dataArray[6][1]) + ";" + sanitizeNode(dataArray[6][-1]) + " ]\n"))
            isInDatabase(valueInput2,dataArray[6])
            return getValueBacteriaBacteria(dataArray[2],dataArray[3],dataArray[8],valueInput1,valueInput2)
        elif (typeInput == "BM"):
            valueInput1 = parseListNode(raw_input("Choose the group of bacterias [ Read the taxonomic tree to help you: e.g. " + sanitizeNode(dataArray[6][-3]) + ";" + sanitizeNode(dataArray[6][1]) + ";" + sanitizeNode(dataArray[6][-1]) + " ]\n"))
            isInDatabase(valueInput1,dataArray[6])
            print dataArray[1]
            valueInput2 = parseList(raw_input("Choose the metadatum among those printed above [ e.g. " + dataArray[1][0] + ";" + dataArray[1][-1] + " ]\n"))
            isInDatabase(valueInput2,dataArray[1])
            return getValueBacteriaMetadata(dataArray[0],dataArray[1],valueInput1,valueInput2)
        else:
            print "\nERROR: You need to answer 'BB' or 'BM', and not ",typeInput
            raise ValueError
    else:
        ()
    
#_____________________________________________________________________________

def pearsonAct(dataArray):
    pearsonType = raw_input("Do you want to compute sample Pearson coefficient or population Pearson coefficient? sample/population \n")
    if (pearsonType == "sample"):
        xArray,yArray = creatingArray(dataArray,True)
        pearson = samplePearson(xArray,yArray)
        print "\nPearson Sample coefficient is: " + str(pearson) + "\n"
    elif (pearsonType == "population"):
        xArray,yArray = creatingArray(dataArray,True)
        printProbabilityLawsList()
        probList = ["UniformProbability","UniformProbabilityProduct"]
        p1 = raw_input("Enter the law of probability for first values [among the ones above]\n")
        isInDatabase(p1,probList)
        p2 = raw_input("Enter the law of probability for second values [among the ones above]\n")
        isInDatabase(p2,probList)
        p3 = raw_input("Enter the law of probability for the product of first values with second values [among the ones above]\n")
        isInDatabase(p3,probList)
        pearson = populationPearson(xArray,yArray,p1,p2,p3)
        print "\nPearson Population coefficient is: " + str(pearson) + "\n"
    else:
        print "\n/!\ ERROR: You need to answer 'sample' or 'population', and not ",pearsonType
        raise ValueError
    answer = raw_input("Save the results? Y/N\n")
    if (answer == "Y"):
        data = "The " + pearsonType + " Pearson coefficient for values: ****\n\n" + str(xArray) + "\n and " + str(yArray) + "\n is : " + str(pearson) + "\n\nEND OF FILE ****"
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

#____________________________________________________________________________

def printTreeAct(dataArray):
    answer = raw_input("Do you want to print sample hit lists? Y/N\n")
    if not ((answer == "Y") or (answer == "N")):
        print "\n/!\ ERROR: You need to answer 'Y' or 'N'."
        raise ValueError
    printTree(dataArray[7],(answer == "Y"))

#____________________________________________________________________________

def plottingAct(dataArray):
    ()

#____________________________________________________________________________

def distanceAct(dataArray):
    ()
