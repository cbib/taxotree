from __future__ import division
import numpy as np
import re

from misc import mem
from parsingInfo import parseInfo

integer = re.compile("[0-9]+")

#@matrix contains similarity coefficients comprised between 0 and 1.
#The more the coefficient is close to 1, the more the samples are alike.
#Only on the values of the metadata in @metadataList
def similarity(samplesList,infoList,metadataList):
    #@n is the number of samples/patients
    #First information is sample ID
    n = len(samplesList)
    m = len(infoList)
    maximum = 0
    matrix = np.zeros((n,n))
    for i in range(n):
        for j in range(i+1,n):
            s = 0
            for k in range(2,m):
                #If this metadatum belongs to the list of interest
                if mem(infoList[k],metadataList):
                    #for all i, len(samplesList[i]) == len(infoList) as asserted in parsingInfo
                    if (integer.match(samplesList[i][k]) and integer.match(samplesList[j][k])):
                        s += abs(int(samplesList[i][k]) - int(samplesList[j][k]))
            matrix[i][j] = s
            matrix[j][i] = s
            maximum = max(maximum,s)
    #The matrix is symmetric
    for i in range(n):
        matrix[i][i] = 1
    for i in range(n):
        for j in range(i+1,n):
            s = matrix[i][j]
            matrix[i][j] = (maximum-s)/maximum
            matrix[j][i] = (maximum-s)/maximum
    return matrix
