import numpy as np

from totalratio import compute,countAssignmentsInCommon,countAssignments,totalRatio,totalRatioNormalized
from patternRatio import patternRatio,enumerateCommonPatterns,enumerateSpecificPatterns
from similarityCoefficient import similarity
from normalization import normalizeSymmetricMatrix

#computes weighted sum of some of the previous calculi: total ratio, pattern ratio and similarity coefficient, between two samples
#Returns a matrix M such as M[i][j] = M[j][i] is the similarity coefficient (S-E(S))/sqrt(V(S)) where S = total ratio + pattern ratio + metadata similarity coefficient [to compute E(S) and V(S), we do not take into account infinite values]
#NB. Distance can be obtained with the formula: distance = 1/similarity
def computeSimilarity(dataArray):
    sampleIDList = dataArray[8]
    n = len(sampleIDList)
    matrix = np.zeros((n,n))
    m = similarity(dataArray[0],dataArray[1])
    for i in range(n):
        for j in range(n):
            #Total ratio computation
            common,in1,in2,_,_,_,_,_ = compute(dataArray[7],[sampleIDList[i]],[sampleIDList[j]])
            commonA = countAssignmentsInCommon(common,[sampleIDList[i]],[sampleIDList[j]])
            numberA1 = countAssignments(in1,[sampleIDList[i]])
            numberA2 = countAssignments(in2,[sampleIDList[j]])
            tratio = totalRatio(commonA,numberA1,numberA2)
            #Pattern ratio computation
            commonPatternsList = enumerateCommonPatterns(dataArray[7],[sampleIDList[i]],[sampleIDList[j]])
            specificPatternsList1 = enumerateSpecificPatterns(dataArray[7],[sampleIDList[i]],[sampleIDList[j]])
            specificPatternsList2 = enumerateSpecificPatterns(dataArray[7],[sampleIDList[i]],[sampleIDList[j]])
            pRatio = patternRatio(commonPatternsList,specificPatternsList1,specificPatternsList2)
            #Get similarity coefficient
            similarityCoefficient = m[i][j] #= m[j][i] (see similarityCoefficient.py)
            matrix[i][j] = pRatio + similarityCoefficient + tratio
            matrix[j][i] = pRatio + similarityCoefficient + tratio
    return normalizeSymmetricMatrix(matrix)

#@sampleNameList is a list of disjoint groups of samples (groups induced by metadata for example, see @computeSampleInGroup in misc).
#Gives the pairs of most different groups of samples (that is, those whose similarity coefficient is inferior to the one of the first quartile)
#The similarity coefficient of two groups of samples G1 and G2 is the sum of the similarity coefficients of any pair of samples (s1,s2) such as s1 belongs to G1, and s2 belongs to G2
def mostDifferentSamplesGroups(matrix,sampleIDList,sampleNameList):
    #Computes the similarity matrix for the groups of @sampleNameList
    m = len(sampleIDList)
    n = len(sampleNameList)
    matrixGroup = np.zeros((n,n))
    #Gets all the values (with duplicates) of the matrix (matrix is symmetric)
    valuesSet = []
    for i in range(n):
        for j in range(i,n):
            s = 0
            for samplei in sampleNameList[i]:
                for samplej in sampleNameList[j]:
                    k = 0
                    while k < m and not (samplei == sampleIDList[k]) and not (samplej == sampleIDList[k]):
                        k += 1
                    l = k
                    while l < m and not (samplei == sampleIDList[l]):
                        l += 1
                    while l < m and not (samplej == sampleIDList[l]):
                        l += 1
                    if l == m:
                        print "\n/!\ ERROR: Samples",samplei,"and",samplej,"not found."
                        raise ValueError
                    #matrix is symmetric
                    s += matrix[l][k]
            matrixGroup[i][j] = s
            matrixGroup[j][i] = s
            valuesSet.append(s)
    #Final list of pairs
    result = []
    #Number of values (matrix is symmetric)
    number = (n*n)/2
    #position of first quartile
    pos = number/4
    valuesSet = sorted(valuesSet,key=lambda x:x)
    #Value of first quartile
    quartile = valuesSet[pos]
    for i in range(n):
        for j in range(i,n):
            if (matrixGroup[i][j] <= quartile):
                #See computeSimilarity below
                result.append((sampleNameList[i],sampleNameList[j]))
    return result
