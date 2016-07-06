import numpy as np
from time import time

from totalRatio import compute,countAssignmentsInCommon,countAssignments,totalRatio,totalRatioNormalized
from patternRatio import patternRatio,enumerateCommonPatterns,enumerateSpecificPatterns
from diversityCoefficient import computeDiversityCoefficient
from misc import inf

#computes weighted sum of some of the previous calculi: total ratio, pattern ratio and similarity coefficient, between two samples
#Returns a matrix M such as M[i][j] = M[j][i] is the similarity coefficient: total ratio + pattern ratio + metadata similarity coefficient 
#NB. Distance can be obtained with the formula: distance = 1/similarity

#Adds finite non-negative values and infinite values
#returns -1 when the result is infinite
def sumOpInf(a,b):
    if (a == "+inf") or (b == "+inf"):
        return inf
    else:
        return a+b

#@dataArray = [samplesInfoList,infoList,samplesOccList,speciesList,paths,n,nodesList,taxoTree,sampleIDList,#similarityMatrix]
    
def computeSimilarity(dataArray):
    start = time()
    sampleIDList = dataArray[8]
    n = len(sampleIDList)
    matrix = np.zeros((n,n))
    for i in range(n):
        for j in range(i,n):
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
            dRatio1,_ = computeDiversityCoefficient(dataArray[5],[sampleIDList[i]],dataArray)
            dRatio2,_ = computeDiversityCoefficient(dataArray[5],[sampleIDList[j]],dataArray)
            subdRatio = abs(dRatio1 - dRatio2)
            if subdRatio:
                s = sumOpInf(pRatio,tratio) - subdRatio
            else:
                s = sumOpInf(pRatio,tratio)
            matrix[i][j] = s
            matrix[j][i] = s
    end = time()
    print "TIME:",(end-start)
    return matrix

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
                    #Positions k, l are the positions of samples in sampleIDList
                    k,l = 0,0
                    while k < m and not (samplei == sampleIDList[k]):
                        k += 1
                    if (k == m):
                        print "\n/!\ ERROR: Sample",samplei,"not found."
                        raise ValueError
                    while l < m and not (samplej == sampleIDList[l]):
                        l += 1
                    if l == m:
                        print "\n/!\ ERROR: Sample",samplej,"not found."
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
            #The "=" is important, as we may hit the infinite part of the range ("-1")
            if (matrixGroup[i][j] <= quartile):
                #See @computeSimilarity below
                #sampleNameList[i] contains the whole line associated to the ith sample
                result.append((sampleNameList[i],sampleNameList[j]))
    return result
