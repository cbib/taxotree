import matplotlib.pyplot as plt
import numpy as np

#@pearson is the Pearson coefficient
def plotPearsonGraph(xArray,yArray,pearson,xLabel="X",yLabel="f(X)",maxx=10,minx=0,maxy=10,miny=0,title="Plotting of unknown function f"):
    mini = min(minx,miny)
    n = len(xArray)
    if not (n == len(yArray)):
        print "\n/!\ ERROR: Different lengths: xArray",n,"and yArray",len(yArray),"."
        raise ValueError
    plt.grid(True)
    plt.title(title)
    plt.xlabel(xLabel)
    plt.ylabel(yLabel)
    plt.xlim(mini,maxx)
    plt.ylim(mini,maxy)
    t = np.zeros(n)
    pearsont = np.zeros(n)
    step = int(max(maxx,maxy)/n) + 1
    currVal = 0
    for i in range(n):
        t[i] = currVal
        pearsont[i] = currVal*pearson
        currVal += step
    #Lines will be in red for the main function, blue for pearson line
    plt.plot(xArray,yArray,"ro",t,pearsont,"b--")
    plt.show()

#Draws points
#len(xArray) == len(yArray)
#xArray and yArray are the array of values for the two variables
#xLabel and yLabel are the corresponding labels
def plotGraph(xArray,yArray,xLabel="X",yLabel="f(X)",maxx=10,minx=0,maxy=10,miny=0,title="Plotting of unknown function f"):
    n = len(xArray)
    if not (n == len(yArray)):
        print "\n/!\ ERROR: Different lengths: xArray",n,"and yArray",len(yArray),"."
        raise ValueError
    plt.grid(True)
    plt.title(title)
    plt.xlabel(xLabel)
    plt.ylabel(yLabel)
    plt.xlim(minx,maxx)
    plt.ylim(miny,maxy)
    #Lines will be in red
    plt.plot(xArray,yArray,"ro")
    plt.show()

#Draws histograms
def plotHist(xArray,xLabel="X",yLabel="f(X)",maxx=10,minx=0,maxy=10,miny=0,title="Histogram of unknown function f"):
    #Green color
    plt.hist(xArray,bins=50,normed=1,facecolor="g",alpha=0.5,label=xLabel)
    plt.grid(True)
    plt.title(title)
    plt.xlabel(xLabel)
    plt.ylabel(yLabel)
    plt.xlim(minx,maxx)
    plt.ylim(miny,maxy)
    plt.show()

#@labels is the array containing the labels of the pie (can go up to 14 different labels)
#@sizes is the arrays of parts of the pie owned by the different labels
def plotPie(labels,sizes,title):
    initColors = ['gold','yellowgreen','lightcoral','lightskyblue','violet','blue','pink','red','orange','green','gray','black','brown','yellow']
    n = len(labels)
    if not (n == len(sizes)):
        print "\n/!\ ERROR: Different lengths ",len(labels),"and",len(sizes)
        raise ValueError
    if n > 14:
        print "\n/!\ ERROR: Not enough colors! Please modify plottingValues.py and restart Python"
        raise ValueError
    #explode maximum percentage
    iMax = 0
    maximum = 0
    for i in range(n):
        if maximum < sizes[i]:
            iMax = i
            maximum = sizes[i]
    explode = [0] * n
    explode[iMax] = 0.1
    labels = labels
    sizes = sizes
    colors = initColors[:n]
    plt.pie(sizes,explode=explode,labels=labels,colors=colors,autopct='%1.1f%%',shadow=True,startangle=140)
    plt.axis('equal')
    plt.title(title)
    plt.show()
