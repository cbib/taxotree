#Centralization-reduction for a list of values
import numpy as np

from misc import inf

#Hypothesis of uniform probability for the occurrence of any bacteria whatever the clinic data may be (which is a strong hypothesis...)
def expectList(vArray):
    n = len(vArray)
    if not n:
        print "\n/!\ ERROR: Empty list."
        raise ValueError
    exp = 0
    for i in range(n):
        exp += vArray[i]/n
    return exp

def standardDeviationList(vArray):
    vProductArray = [x*x for x in vArray]
    expProd = expectList(vProductArray)
    exp = expectList(vArray)
    expS = exp*exp
    return np.sqrt(expProd-expS)

def normalizeList(valueList):
    exp = expectList(valueList)
    stDeviation = standardDeviationList(valueList)
    if not stDeviation:
        print "\n/!\ ERROR: Math problem (Division by zero)."
        raise ValueError
    normList = []
    for value in valueList:
        normList.append((value-exp)/stDeviation)
    return normList

#__________________________________________________________________________________________

#Idem for 2D matrix with finite and infinite values
def expectMatrix(matrix,n,m):
    exp = 0
    if not n or not m:
        print "\n/!\ ERROR: Math problem (Division by zero)."
        raise ValueError
    for i in range(n):
        for j in range(m):
            try:
            #All values are theoretically non-negative
            #Since we have to deal only with integers (numpy module in Python...)
            #-1 signifies infinite value (see computeDiscriminatoryDistance.py)
            #We ignore these infinite values in the calculus of the expectation
                if not (matrix[i][j] == inf):
                    exp += matrix[i][j]/(n*m)
            except IndexError:
                print "BUG"
                print matrix[i],matrix[j]
                from misc import getSampleIDList
                from parsingInfo import parseInfo
                sampleInfoList,_ = parseInfo("Info")
                sampleIDList = getSampleIDList(sampleInfoList)
                print sampleIDList[i],sampleIDList[j]
                print "END BUG"
    return exp

def standardDeviationMatrix(matrix,n,m):
    mProduct = [matrix[i][j]*matrix[i][j] for i in range(n) for j in range(m)]
    expProd = expectMatrix(mProduct,n,m)
    exp = expectMatrix(matrix,n,m)
    expS = exp*exp
    return np.sqrt(expProd-expS)

def normalizeSymmetricMatrix(valueMatrix):
    n,m = np.shape(valueMatrix)
    exp = expectMatrix(valueMatrix,n,m)
    stDeviation = standardDeviationMatrix(valueMatrix,n,m)
    if not stDeviation:
        print "\n/!\ ERROR: Math problem (Division by zero)."
        raise ValueError
    normMatrix = np.zeros((n,m))
    for i in range(n):
        for j in range(i,m):
            #Same remark as before
            if valueMatrix[i][j] == inf:
                normMatrix[i][j] = inf
                normMatrix[i][j] = inf
            else:
                normMatrix[i][j] = (valueMatrix[i][j]-exp)/stDeviation
                normMatrix[j][i] = (valueMatrix[i][j]-exp)/stDeviation
    return normMatrix
            
