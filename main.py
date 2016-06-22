from parsingTree import parseTree
from parsingMatrix import parseMatrix
from parsingInfo import parseInfo
from taxoTree import TaxoTree
from actions import totalDiffRatioAct,patternRatioAct,percentageAct,pearsonAct,similarityAct,printTreeAct,plottingAct,distanceAct
from misc import getSampleIDList

def main():
    tTree = raw_input("Write down the file name of the taxonomic tree in the folder \"meta\" [ without the extension .tree ]\n")
    if (tTree == ""):
        tTree = "GGdb2015"
    oMatrix = raw_input("Write down the CSV file name of the occurrence matrix in the folder \"meta\" [ without the extension .csv ]\n")
    if (oMatrix == ""):
        oMatrix = "MGAcount_complete"
    iMatrix = raw_input("Write down the CSV file name of the data matrix in the folder \"meta\" [ without the extension .csv ]\n")
    if (iMatrix == ""):
        iMatrix = "Info"
    print "/!\ Data getting parsed..."
    try:
        samplesInfoList,infoList = parseInfo(iMatrix)
        sampleIDList = getSampleIDList(samplesInfoList)
    except IOError:
        print "\nERROR: Maybe the filename you gave does not exist in \"meta\" folder\n"
    print "..."
    try:
        samplesOccList,speciesList = parseMatrix(oMatrix)
    except IOError:
        print "\nERROR: Maybe the filename you gave does not exist in \"meta\" folder\n"
    print "..."
    try:
        paths,n,nodesList = parseTree(tTree)
    except IOError:
        print "\nERROR: Maybe the filename you gave does not exist in \"meta\" folder\n"
    print "-- End of parsing\n"
    print "/!\ Constructing the whole annotated taxonomic tree"
    print "[ You may have to wait for a few seconds... ]"
    taxoTree = TaxoTree("Root").addNode(paths,nodesList,samplesOccList)
    print "-- End of construction\n"
    dataArray = [samplesInfoList,infoList,samplesOccList,speciesList,paths,n,nodesList,taxoTree,sampleIDList]
    answer = ""
    while not ((answer == "exit") or (answer == "exit()") or (answer == "quit")):
        try:
            print "What do you want to do?"
            print "[Write down the number matching with the action required. Details are in README file]"
            print "   1: Total ratio and Diff ratio"
            print "   2: Pattern ratio"
            print "   3: Percentage of assignments in a certain group of bacterias depending on metadata"
            print "   4: Pearson correlation coefficient"
            print "   5: Similarity coefficients between patients"
            print "   6: Print the taxonomic tree"
            print "   7: Plot a graph, or a pie"
            print "   8: Compute total distance between two samples"
            print "[To quit, write down exit]"
            answer = raw_input("Your answer?\n")
            if (answer =="1"):
                totalDiffRatioAct(dataArray)
                print "-- End \n"
            elif (answer == "2"):
                patternRatioAct(dataArray)
                print "-- End \n"
            elif (answer == "3"):
                percentageAct(dataArray)
                print "-- End \n"
            elif (answer == "4"):
                pearsonAct(dataArray)
                print "-- End \n"
            elif (answer == "5"):
                matrixSim = similarityAct(dataArray,iMatrix)
                dataArray.append(matrixSim)
                print "-- End \n"
            elif (answer == "6"):
                printTreeAct(dataArray)
                print "-- End \n"
            elif (answer == "7"):
                plottingAct(dataArray)
                print "-- End \n"
            elif (answer == "8"):
                distanceAct(dataArray)
                print "-- End \n"
            elif not ((answer == "exit") or (answer == "exit()") or (answer == "quit")):
                print "/!\ ERROR: Please enter a number between 1 and 8 included, or 'exit' if you want to quit."
                raise ValueError
        except ValueError:
            print "/!\ ERROR: Please look at the line above."
            print "/!\ ERROR: If the line above is blank, it may be an uncatched ValueError.\n"
    #return dataArray
    
