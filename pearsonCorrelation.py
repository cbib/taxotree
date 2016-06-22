from __future__ import division

#Allows the measure of linear correlation between two metadata x and y
#Provided two arrays containing:
#- x values (x array being for instance the number of assigned species in a certain group of bacterias in a single sample)
#- y values associated (values for another sample)
#Returns the expectation of all the Pearson product-moment correlation coefficients (between +1 and -1 inclusive) obtained, under hypothesis of independance of samples and uniform probability (to be improved, since the sample depends a priori on all the metadata provided)

import numpy as np

from misc import truncate

def printProbabilityLawsList():
    ls = ["uniform","uniformProduct"]
    print ls
    return ls

def uniformProbability(xArray):
    n = len(xArray)
    if not n:
        print "\n/!\ ERROR: Empty list."
        raise ValueError
    probArray = []
    for i in range(n):
        probArray.append(1/n)
    return probArray

def uniformProbabilityProduct(xArray):
    n = len(xArray)
    if not n:
        print "\n/!\ ERROR: Empty list."
        raise ValueError
    probArray = []
    for i in range(n):
        probArray.append(1/(n*n))
    return probArray

def applyProb(xArray,p):
    if (p == "uniform"):
        return uniformProbability(xArray)
    elif (p == "uniformProduct"):
        return uniformProbabilityProduct(xArray)
    else:
        print "\n/!\ ERROR: Undefined probability law."
        raise ValueError

#Calculus of expectation provided the array of every probability associated to eeach value of the xArray
#xArray is an list containing (sample/metadatum,number) pairs
def expectation(xArray,probArray):
    s = 0
    n = len(xArray)
    if not (n == len(probArray)):
            print "\n/!\ ERROR: [BUG] [pearsonCorrelation/expectation] Different lengths xArray",n,"yArray",len(yArray)
    for i in range(n):
        s += xArray[i][1]*probArray[i]
    return s

#p (p1, p2 for multiple variables) is the probability function which assigns to each x value its probability
def variation(xArray,p):
    n = len(xArray)
    xSqArray = []
    for i in range(n):
        xSqArray.append((xArray[i][0],xArray[i][1]*xArray[i][1]))
    expSq = expectation(xSqArray,applyProb(xSqArray,p))
    exp = expectation(xArray,applyProb(xArray,p))
    return (expSq-exp*exp)
    
def stDeviation(xArray,p):
    return np.sqrt(variation(xArray,p))

def covariance(xArray,yArray,p1,p2,p3):
    expX = expectation(xArray,applyProb(xArray,p1))
    expY = expectation(yArray,applyProb(yArray,p2))
    n = len(xArray)
    if not (n == len(yArray)):
        print "\n/!\ ERROR: [BUG] [pearsonCorrelation/covariance] Different lengths: xArray",n,"yArray",len(yArray)
        raise ValueError
    xyArray = []
    for i in range(n):
        xyArray.append((xArray[i][0],xArray[i][1]*yArray[i][1]))
    expXY = expectation(xyArray,applyProb(xyArray,p3))
    return (expXY - expX*expY)

def populationPearson(xArray,yArray,pp1,pp2,pp3):
    cov = covariance(xArray,yArray,pp1,pp2,pp3)
    stdX = stDeviation(xArray,pp1)
    stdY = stDeviation(yArray,pp2)
    #If stdX = 0 or stdY = 0
    if not stdX or not stdY:
        #meaning there is no linear correlation exists between the two variables
        result = 0
    else:
        result = cov/(stdX*stdY)
    if result > 1 or result < -1:
        print "\n/!\ ERROR: Inconsistent value of Pearson correlation coefficient:",result,"(should be between -1 and 1)"
        raise ValueError
    return result

def mean(xArray):
    n = len(xArray)
    if not n:
        print "\n/!\ ERROR: Empty list."
        raise ValueError
    s = 0
    for i in range(n):
        s += xArray[i][1]
    return ((1/n)*s)

def samplePearson(xArray,yArray):
    #Must be uncommented if the original data has not been normalized
    #mX = mean(xArray)
    #mY = mean(yArray)
    mX,mY = 0,0
    n = len(yArray)
    s1,s2,s3 = 0,0,0
    for i in range(n):
        s1 += (xArray[i][1] - mX)*(yArray[i][1] - mY)
        s2 += (xArray[i][1] - mX)*(xArray[i][1] - mX)
        s3 += (yArray[i][1] - mY)*(yArray[i][1] - mY)
    if (s2 == 0) or (s3 == 0):
        print "\n/!\ ERROR: No assignment in these samples."
        raise ValueError
    #Truncature needed for extreme cases (1 and -1)
    result = truncate(s1/(np.sqrt(s2)*np.sqrt(s3)),1)
    if (result > 1) or (result < -1):
        print "%.25f"%result,(result > 1),(result < -1)
        print "\n/!\ ERROR: Inconsistent value of Pearson correlation coefficient:",result,"(should be between -1 and 1)"
        raise ValueError
    return result
