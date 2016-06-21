#Centralization-reduction for a list of values
import numpy as np

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

#Idem for 2D matrix with finite and infinite values
def expectMatrix(matrix,n,m):
    exp = 0
    if not n or not m:
        print "\n/!\ ERROR: Math problem (Division by zero)."
        raise ValueError
    for i in range(n):
        for j in range(m):
            if not (matrix[i][j] == "+inf"):
                exp += matrix[i][j]/(n*m)
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
            normMatrix[i][j] = (valueMatrix[i][j]-exp)/stDeviation
            normMatrix[j][i] = (valueMatrix[i][j]-exp)/stDeviation
    return normMatrix
            
