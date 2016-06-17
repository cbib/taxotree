import re
from misc import sanitize

integer = re.compile("[0-9]+")

#Parse a CSV file
def parseMatrix(filename):
    speciesList = []
    samplesList = []
    file_matrix = open("meta/" + filename + ".csv","r")
    lines = file_matrix.readlines()
    file_matrix.close()
    boolean = True
    for line in lines:
        ls = line.split(",")
        #First line gives the name and rank of species in the samples
        if boolean:
            #ls is then a list of strings of type "rank:name"
            #Turns "rank:name" into (name,rank)
            for string in ls:
                ls1 = string.split(":")
                rank = sanitize(ls1[0])
                #Deletes the white space after name
                #Otherwise equality on strings does not work
                name = sanitize(ls1[-1])
                speciesList.append((name,rank))
            boolean = False
            n = len(speciesList)
        else:
            thisSampleList = []
            for number in ls:
                number = sanitize(number)
                if integer.match(number):
                    thisSampleList.append(int(number))
                else:
                    thisSampleList.append(number)
            if not (len(thisSampleList) == n):
                raise ValueError
            samplesList.append(thisSampleList)
    return samplesList,speciesList

#Provided a specie, gives the sampleHitList associated (see TaxoTree: nodesList is the list of nodes present in the tree)
def associatedData(nodesList,samplesList,name,rank):
    n = len(nodesList)
    sampleHitList = []
    #First element is "SampleID"
    i = 1
    while (i < n) and not ((nodesList[i][0] == name) and (nodesList[i][1] == rank)):
        i += 1
    if (i == n):
        return []
    else:
        for sample in samplesList:
            n = len(sample)
            if (i < n) and not (sample[i] == 0):
                #Memorizes (name of the sample, number associated to specie in this sample)
                sampleHitList.append((sample[0],sample[i]))
        return sampleHitList
