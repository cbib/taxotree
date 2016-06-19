from misc import sanitize,containsSpecie

#Module for regular expressions
import re

#Useful macros
def getNameRankList(string):
    return [(sanitize(string.split(":")[-1]).split("\n")[0],sanitize(string.split(":")[0]).split("\n")[0])]

def boolConditions(regexpLs,numberLine,lines,n):
    #if the line matches the empty line
    b = not ((numberLine < n) and re.compile("^$").match(lines[numberLine]))
    for reg in regexpLs:
        b = b or not ((numberLine < n) and reg.match(lines[numberLine]))
    return b

#Produces a list containing the path to every leaf of the taxonomic tree
#@paths is a list of (name,rank) lists
#@allNodes to True gives every path from the root to a(n internal) node
def parseTree(filename,allNodes=False):
    paths = []
    file_tree = open("meta/" + filename + ".tree")
    lines = file_tree.readlines()
    file_tree.close()
    n = len(lines)
    #Count number of read lines
    i1 = 0
    #Regex associated with respectively reign, domain, phylum, class, order, family, genus, species
    r = re.compile("^R:")
    k = re.compile("^ K:")
    p = re.compile("^  P:")
    c = re.compile("^   C:")
    o = re.compile("^    O:")
    f = re.compile("^     F:")
    g = re.compile("^      G:")
    s = re.compile("^       S:")
    #That code might be factorizable using arrays to store regexp and lineages
    #(but would not be as readable)
    lineage0,lineage11 = [],[]
    lineage,lineage2,lineage6 = [],[],[]
    lineage3,lineage4,lineage5 = [],[],[]
    #@nodesList is the list of all nodes present in the tree
    nodesList = []
    ranks = [r,k,p,c,o,f,g,s]
    while i1 < n:
        #If the current read line matches r regexp
        if r.match(lines[i1]):
            #Adds to the current path name and rank of the bacteria
            lineage11 = getNameRankList(lines[i1])
            nodesList.append(lineage11[-1])
            i1 += 1
            #If we need every path from the root to a node
            #Or if this current bacteria is a leaf in the taxonomic tree
            if allNodes or boolConditions(ranks[1:],i1,lines,n):
                paths += [ lineage11 ]
        #If the current read line matches k regexp
        elif k.match(lines[i1]):
            lineage0 = lineage11 + getNameRankList(lines[i1])
            nodesList.append(lineage0[-1])
            i1 += 1
            if allNodes or boolConditions(ranks[2:],i1,lines,n):
                paths += [ lineage0 ]
        elif p.match(lines[i1]):
            lineage = lineage0 + getNameRankList(lines[i1])
            nodesList.append(lineage[-1])
            i1 += 1
            if allNodes or boolConditions(ranks[3:],i1,lines,n):
                paths += [ lineage ]
        elif c.match(lines[i1]):
            lineage2 = lineage + getNameRankList(lines[i1])
            nodesList.append(lineage2[-1])
            i1 += 1
            if allNodes or boolConditions(ranks[4:],i1,lines,n):
                paths += [ lineage2 ]
        elif o.match(lines[i1]):
            lineage3 = lineage2 + getNameRankList(lines[i1])
            nodesList.append(lineage3[-1])
            i1 += 1
            if allNodes or boolConditions(ranks[5:],i1,lines,n):
                paths += [ lineage3 ]
        elif f.match(lines[i1]):
            lineage4 = lineage3 + getNameRankList(lines[i1])
            nodesList.append(lineage4[-1])
            i1 += 1
            if allNodes or boolConditions(ranks[6:],i1,lines,n):
                paths += [ lineage4 ]
        elif g.match(lines[i1]):
            lineage5 = lineage4 + getNameRankList(lines[i1])
            nodesList.append(lineage5[-1])
            i1 += 1
            if allNodes or boolConditions(ranks[7:],i1,lines,n):
                paths += [ lineage5 ]
        elif s.match(lines[i1]):
            lineage6 = lineage5 + getNameRankList(lines[i1])
            nodesList.append(lineage6[-1])
            i1 += 1
            paths += [ lineage6 ]
        else:
            print "\n /!\ ERROR: [BUG] [parsingTree/parseTree] Parsing error."
            raise ValueError
    #Returns list of paths, total number of nodes, list of all nodes
    return paths,n,nodesList
