import matplotlib.pyplot as plt
import numpy as np

#Draws points
#xArray and yArray are the array of values for the two variables
#xLabel and yLabel are the corresponding labels
def plotGraph(xArray,yArray,xLabel="X",yLabel="f(X)",maxx=10,minx=0,maxy=10,miny=0,title="Plotting of unknown function f"):
    fig = plt.figure()
    plt.grid(True)
    plt.title(title)
    plt.xlabel(xLabel)
    plt.ylabel(yLabel)
    plt.xlim(minx,maxx)
    plt.ylim(miny,maxy)
    #Lines will be in red
    plt.plot(xArray,yArray,"ro")
    answer = raw_input("Show this figure? Y/N\n")
    if (answer == "Y"):
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
    answer = raw_input("Show this figure? Y/N\n")
    if (answer == "Y"):
        plt.show()
