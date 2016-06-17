from misc import sanitize

def parseInfo(filename):
    samplesList = []
    file_matrix = open("meta/" + filename + ".csv","r")
    lines = file_matrix.readlines()
    file_matrix.close()
    #Data need to be sanitized
    infoListDirty = lines[0].split(",")
    infoList = []
    for info in infoListDirty:
        infoList.append(sanitize(info.split("\n")[0]))
    for line in lines[1:]:
        #Construction of the list associated to one sample
        thisSampleList = []
        lsDirty = line.split(",")
        #Checks if lsDirty is not empty
        assert len(lsDirty)
	ls = []
        for data in lsDirty:
            if not (data == ""):
                ls.append(data)
            else:
                #unknown values are remplaced by "N"
                ls.append("N")
        for data in ls:
            thisSampleList.append(sanitize(data).split("\n")[0])
        #samplesList is the list of every sample's list
	assert (len(thisSampleList) == len(infoList))
        samplesList.append(thisSampleList)
    return samplesList,infoList

