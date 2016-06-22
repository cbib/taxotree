from __future__ import division

#Allows the measure of linear correlation between two metadata x and y
#Provided two arrays containing:
#- x values (x array being for instance the number of assigned species in a certain group of bacterias in a single sample)
#- y values associated (values for another sample)
#Returns the expectation of all the Pearson product-moment correlation coefficients (between +1 and -1 inclusive) obtained, under hypothesis of independance of samples and uniform probability (to be improved, since the sample depends a priori on all the metadata provided)

import numpy as np

from misc import truncate

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
