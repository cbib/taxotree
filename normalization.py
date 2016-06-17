#Centralization-reduction for a list of values
import numpy as np

#Hypothesis of uniform probability for the occurrence of any bacteria whatever the clinic data may be (which is a strong hypothesis...)
def expect(vArray):
    n = len(vArray)
    exp = 0
    for i in range(n):
        exp += vArray[i]/n
    return exp

def standardDeviation(vArray):
    vProductArray = [x*x for x in vArray]
    expProd = expect(vProductArray)
    exp = expect(vArray)
    expS = exp*exp
    return np.sqrt(expProd-expS)

def normalize(valueList):
    exp = expect(valueList)
    stDeviation = standardDeviation(valueList)
    normList = []
    for value in valueList:
        normList.append((value-exp)/stDeviation)
    return normList
            
            
