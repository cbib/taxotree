from __future__ import division
import numpy as np
import re

from writeOnFiles import writeFile
from parsingMatrix import parseMatrix,importMatrix
from parsingInfo import parseInfo
from parsingTree import parseTree
from taxoTree import TaxoTree,printTree
from misc import getValueBacteriaBacteria,getValueBacteriaMetadata,mem,isInDatabase,getMaxMin,partitionSampleByMetadatumValue

from totalRatio import compute,countAssignmentsInCommon,countAssignments,totalRatio,totalRatioNormalized,diffRatio,diffRatioNormalized
from patternRatio import patternRatio,enumerateCommonPatterns,enumerateSpecificPatterns
from pearsonCorrelation import samplePearson
from percentage import percentageAssign,computeSamplesInGroup
from similarityCoefficient import similarity
from computeDiscriminatoryDistance import computeSimilarity,mostDifferentSamplesGroups
from plottingValues import plotPearsonGraph,plotGraph,plotHist,plotPie

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
    elif not (answer == "N"):
        print "/!\ You should answer 'Y' or 'N'!"

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
    elif not (answer == "N"):
        print "/!\ You should answer 'Y' or 'N'!"
        
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
    elif not (answer == "N"):
        print "/!\ You should answer 'Y' or 'N'!"
    return result,nodesGroup,sampleNameList,metadataList

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
            xArray,yArray = getValueBacteriaBacteria(dataArray[2],dataArray[3],dataArray[8],valueInput1,valueInput2)
            return xArray,yArray,typeInput,valueInput1,valueInput2
        elif (typeInput == "BM"):
            valueInput1 = parseListNode(raw_input("Choose the group of bacterias [ Read the taxonomic tree to help you: e.g. " + sanitizeNode(dataArray[6][-3]) + ";" + sanitizeNode(dataArray[6][1]) + ";" + sanitizeNode(dataArray[6][-1]) + " ]\n"))
            isInDatabase(valueInput1,dataArray[6])
            print dataArray[1]
            valueInput2 = parseList(raw_input("Choose the metadatum among those printed above [ e.g. " + dataArray[1][0] + ";" + dataArray[1][-1] + " ]\n"))
            isInDatabase(valueInput2,dataArray[1])
            xArray,yArray = getValueBacteriaMetadata(dataArray[0],dataArray[1],valueInput1,dataArray[8],dataArray[2],dataArray[3],valueInput2)
            return xArray,yArray,typeInput,valueInput1,valueInput2
        else:
            print "\nERROR: You need to answer 'BB' or 'BM', and not ",typeInput
            raise ValueError
    #Available cases for only plotting graphs
    else:
        graphTypeInput = raw_input("Do you want to plot a graph or a pie? graph/pie [Read README for details. Histograms will be available in later versions]\n")
        if graphTypeInput == "graph":
            typeInput = raw_input("Do you want to plot bacteria/bacteria or bacteria/metadatum? BB/BM [ Please read README for details. ]\n")
            if (typeInput == "BB"):
                valueInput1 = parseListNode(raw_input("Choose the first group of bacterias [ Read the taxonomic tree to help you: e.g. " + sanitizeNode(dataArray[6][-3]) + ";" + sanitizeNode(dataArray[6][1]) + ";" + sanitizeNode(dataArray[6][-1]) + " ]\n"))
                isInDatabase(valueInput1,dataArray[6])
                valueInput2 = parseListNode(raw_input("Choose the second group of bacterias [ Read the taxonomic tree to help you: e.g. " + sanitizeNode(dataArray[6][-3]) + ";" + sanitizeNode(dataArray[6][1]) + ";" + sanitizeNode(dataArray[6][-1]) + " ]\n"))
                isInDatabase(valueInput2,dataArray[6])
                return graphTypeInput,getValueBacteriaBacteria(dataArray[2],dataArray[3],dataArray[8],valueInput1,valueInput2),typeInput,valueInput1,valueInput2
            elif (typeInput == "BM"):
                valueInput1 = parseListNode(raw_input("Choose the group of bacterias [ Read the taxonomic tree to help you: e.g. " + sanitizeNode(dataArray[6][-3]) + ";" + sanitizeNode(dataArray[6][1]) + ";" + sanitizeNode(dataArray[6][-1]) + " ]\n"))
                isInDatabase(valueInput1,dataArray[6])
                print dataArray[1]
                valueInput2 = parseList(raw_input("Choose the metadatum among those printed above [ e.g. " + dataArray[1][0] + ";" + dataArray[1][-1] + " ]\n"))
                isInDatabase(valueInput2,dataArray[1])
                return graphTypeInput,getValueBacteriaMetadata(dataArray[0],dataArray[1],valueInput1,dataArray[2],dataArray[3],valueInput2),typeInput,valueInput1,valueInput2
            else:
                print "\nERROR: You need to answer 'BB' or 'BM', and not ",typeInput
                raise ValueError
        elif graphTypeInput == "pie":
                result,nodesGroup,sampleNameList,metadataList = percentageAct(dataArray)
                return graphTypeInput,result,nodesGroup,sampleNameList,metadataList
        else:
            print "\nERROR: You need to answer 'graph' or 'pie', and not ",graphTypeInput
            raise ValueError            
    
#_____________________________________________________________________________

def pearsonAct(dataArray):
    xNArray,yArray,typeInput,valueInput1,valueInput2 = creatingArray(dataArray,True)
    if typeInput == "BB":
        xArray = xNArray
        pearson = samplePearson(xArray,yArray)
    #typeInput = "BM"
    else:
        #xNArray is a list of list of (sampleID,number of assignments in this sample) pairs
        #We need thus to sum the numbers of assignments for a same value of the metadata
        xArray = []
        for ls in xNArray:
            s = 0
            for pair in ls:
                s += pair[1]
            xArray.append((ls[0][0],s))
        pearson = samplePearson(xArray,yArray)
    print "\nPearson coefficient is: " + str(pearson) + "\n"
    answer = raw_input("Save the results? Y/N\n")
    if (answer == "Y"):
        data = "The (" + typeInput + ") Pearson coefficient ****\nfor values: \n\n" + str([ x[1] for x in xArray]) + "\ncorresponding to " + str(valueInput1) + "\n\n and\n\n" + str([ y[1] for y in yArray ]) + "\ncorresponding to " + str(valueInput2) + "\n\n is : " + str(pearson) + "\n\nEND OF FILE ****"
        writeFile(data,"","text")
    elif not (answer == "N"):
        print "/!\ You should answer 'Y' or 'N'!"
    answer = raw_input("Plot the corresponding graph? Y/N\n")
    if (answer == "Y"):
        cleanedxArray = [ x[1] for x in xArray ]
        cleanedyArray = [ y[1] for y in yArray ]
        maxix,minix = getMaxMin(cleanedxArray)
        maxiy,miniy = getMaxMin(cleanedyArray)
        #It is more interesting to generate a histogram in these cases
        if len(cleanedxArray) < 4:
            plotHist(cleanedyArray,str(valueInput2[:3]) + "...",str(valueInput1[:3]) + "...",maxiy,miniy,maxix-1,minix+1,"Plotting (" + typeInput + ") Pearson coefficient and the graph of both sets of values")
        else:
            plotPearsonGraph(cleanedxArray,cleanedyArray,pearson,str(valueInput1[:3]) + "...",str(valueInput2[:3]) + "...",maxix,minix,maxiy,miniy,"Plotting (" + typeInput + ") Pearson coefficient and the graph of both sets of values")
    elif not (answer == "N"):
        print "/!\ You should answer 'Y' or 'N'!"

#_____________________________________________________________________________

def similarityAct(dataArray,iMatrix):
    print "/!\ Computing similarity matrix..."
    m = similarity(dataArray[0],dataArray[1])
    print "[Preview.]"
    print m
    answer = raw_input("Save the results? Y/N\n")
    if (answer == "Y"):  
        writeFile(m,"Similarity coefficients between patients for file meta/" + iMatrix + ".csv:\n" + listNodes(dataArray[8]),"array")
    elif not (answer == "N"):
        print "/!\ You should answer 'Y' or 'N'!"
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
    creatingArrayOutput = creatingArray(dataArray)
    if creatingArrayOutput[0] == "graph":
            graphTypeInput,xArray,yArray,typeInput,valueInput1,valueInput2 = creatingArrayOutput
            maxx,minix = getMaxMin(xArray)
            maxy,miniy = getMaxMin(yArray)
            if typeInput == "BB":
                xLabel = "Group of bacteria 1"
                yLabel = "Group of bacteria 2"
            #typeInput == "BM"
            else:
                xLabel = "Group of bacteria"
                yLabel = "Metadatum"
            plotGraph(xArray,yArray,xLabel=xLabel,yLabel=yLabel,maxx=maxx,minx=minix,maxy=maxy,miny=miniy,title="Plotting of type " + typeInput + " with values " + str(valuesInput1[:3]) + "...  and " + str(valuesInput2[:3]) + "...")
    elif creatingArrayOutput[0] == "pie":
        graphTypeInput,result,nodesGroup,sampleNameList,metadataList = creatingArrayOutput
        plotPie(sampleNameList,result,"Assignments to the group of bacterias: " + str(nodesGroup) + " depending on samples")
    else:
        print "\n/!\ ERROR: [BUG] [actions/plottingAct] Unknown type of graph."
        raise ValueError

#____________________________________________________________________________

def distanceAct(dataArray):
    answer = raw_input("Import matrix? Y/N\n")
    if answer == "Y":
        filename = raw_input("Write down the file name where the matrix is stored [ without the extension .taxotree ].\n")
        matrix = importMatrix(filename)
    else:
        if not (answer == "N"):
            print "/!\ You should answer 'Y' or 'N'!"
        print "/!\ Computing similarity matrix..."
        print "[ You may have to wait for a few minutes... ]"
        m = computeSimilarity(dataArray)
        print "[Preview.]"
        print m
        answer = raw_input("Save the results? Y/N\n")
        if (answer == "Y"):  
            writeFile(m,"Similarity coefficients between patients using previous calculi on total ratio, pattern ratio and similarity matrix\n\nNota Bene: 1e+14 stands for +inf\n","array")
        elif not (answer == "N"):
            print "/!\ You should answer 'Y' or 'N'!"
    answer = raw_input("Compute the most different groups of samples?\n")
    if (answer == "Y"):
        print dataArray[1]
        metadatum = parseList(raw_input("Choose the metadatum among those printed above [ e.g. " + dataArray[1][0] + ";" + dataArray[1][-1] + " ]\n"))
        isInDatabase(metadatum,dataArray[1])
        _,valueSampleMetadatum = partitionSampleByMetadatumValue(metadatum,dataArray[1],dataArray[0])
        valueSampleMetadatumNameOnly = []
        for sampleGroup in valueSampleMetadatum:
            sampleGroupNameOnly = []
            for sample in sampleGroup:
                sampleGroupNameOnly.append(sample[0])
            valueSampleMetadatumNameOnly.append(sampleGroupNameOnly)
        pairsList = mostDifferentSamplesGroups(matrix,dataArray[8],valueSampleMetadatumNameOnly)
        print "[ Preview. ]"
        print "List of the pairs of most different sample groups according to the similarity coefficients computed:"
        for pair in pairsList:
            print pair
        answer = raw_input("\nSave the results? Y/N\n")
        if (answer == "Y"):
            stringSamples = ""
            for group in valueSampleMetadatumNameOnly:
                stringSamples += "*" + str(group) + "\n"
            stringPairs = ""
            for pair in pairsList:
                stringPairs += "*" + str(pair) + "\n"
            data = "Most different groups of samples ****\nsorted by values of metadatum: " + metadatum[0] + "\nGroups were:\n\n" + stringSamples + "\n\nAnd the most different ones are:\n\n" + stringPairs + "\n\nEND OF FILE ****"
            writeFile(data,"","text")
        elif not (answer == "N"):
            print "/!\ You should answer 'Y' or 'N'!"
    elif not (answer == "N"):
        print "/!\ You should answer 'Y' or 'N'!"
