import numpy as np

from totalratio import compute,countAssignmentsInCommon,countAssignments,totalRatio,totalRatioNormalized
from patternRatio import patternRatio,enumerateCommonPatterns,enumerateSpecificPatterns
from similarityCoefficient import similarity
from normalization import normalizeSymmetricMatrix

#computes weighted sum of some of the previous calculi: total ratio, pattern ratio and similarity coefficient, between two samples
#Returns a matrix M such as M[i][j] = M[j][i] is the similarity coefficient (S-E(S))/sqrt(V(S)) where S = total ratio + pattern ratio + metadata similarity coefficient [to compute E(S) and V(S), we do not take into account infinite values]
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

#To find: discriminatory criterium
