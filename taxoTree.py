#The entire taxonomy tree

#A node (assignment of a read by Tango) in such a tree contains:
#@name the name of the species
#@rank its rank (k:Kingdom, p: Phylum, c: Class, o: Order, f: Family, g: Genus, s: Species) (Domain can be added easily if needed)
#@ident an identifier for the node, unique in the tree
#@sampleHitList a list that contains (name of sample,number of assignments to this node in this sample) pairs
#@lineage the list of ancestors ((name,rank) list) of the node (from the father to the root)
#@children the list of children (TaxoTree) of the node ordered by growing @ident
#@paths keeps memory of paths to every node/leaf of the tree (to compute more easily distances between nodes of a TaxoTree)
#Other identifier for a node (ident apart) is (@name,@rank)

#Nodes and leaves are only distinguinshed by children=[]

from parsingTree import parseTree
from parsingMatrix import parseMatrix,associatedData
from misc import containsSpecie,mem,selectPath,compare

class TaxoTree(object):
    #
    #
    def __init__(self,name="Root",rank="R",ident=0,sampleHitList=None,lineage=None,children=None,paths=None):
        self.name = name
        self.rank = rank
        self.ident = ident
        self.sampleHitList = sampleHitList or []
        self.lineage = lineage or []
        self.children = children or []
        self.paths = paths or []
        #
        #
    def __eq__(self,name,rank):
        return (self.name == name) and (self.rank == rank)
    #
    #
    #Returns the tree rooted to the node labelled by name, rank
    def search(self,name,rank,DFS=False):
        if self.__eq__(name,rank):
            return self
        else:
            nodeList = []
            for ch in self.children:
                nodeList.append(ch)
            while nodeList:
                node = nodeList.pop()
                if node.__eq__(name,rank):
                    return node
                if DFS:
                    nodeList = node.children + nodeList
                else:
                    nodeList += node.children
            print "\n/!\ ERROR: Subtree rooted at: (",name,rank,") not found. It probably means that this node is not assigned in the samples you have selected."
            raise ValueError
        #
        #
    def __isMember__(self,name,rank,DFS=False):
        if self.__eq__(name,rank):
            return True
        else:
            nodeList = self.children
            while nodeList:
                node = nodeList.pop()
                if node.__eq__(name,rank):
                    return True
                if DFS:
                    nodeList = node.children + nodeList
                else:
                    nodeList += node.children
            return False
        #
        #
    #Debugging function #1
    def printLineage(self,lineage):
        print "-- LINEAGE of %s %s --"%(self.name,self.rank)
        for x in lineage:
            print "(%s,%s)" %(x[0],x[1])
        print "---"
        #
        #
    #Debugging function #2
    def printChildren(self,children):
        print "-- CHILDREN of %s %s --"%(self.name,self.rank)
        for ch in children:
            print "(%s,%s,%d)" %(ch.name,ch.rank,ch.ident)
        print "---"
    #
    #
    def addNodePreProcess(self,paths,nodesList,samplesList,ranks):
        #Firstly, sorting nodes by decreasing rank: R < K < P < C < O < F < G < S
        sortedNodesList = sorted(nodesList,cmp=lambda x,y: compare(x[1],y[1]))
        nodesNumber = len(sortedNodesList)
        #Secondly, getting paths/lineage, and sampleHitList for every node
        pathsCopy = []
        n = 0
        for path in paths:
            n += 1
            pathsCopy.append(path)
        #@selectionPath is an optimized version of @selectPath for the construction of taxonomic trees
        pathsNodes = [ selectPath(pathsCopy,node[0],node[1],n) for node in sortedNodesList ]
        pathsNodesLength = len(pathsNodes)
        sampleHitListNodes = [ associatedData(nodesList,samplesList,node[0],node[1]) for node in sortedNodesList ]
        #Thirdly, clustering brothers (of same rank) into lists for each rank in @allBrotherList
        #Index of current node in sortedNodesList
        currNodeIdent = 0
        #@allBrotherList will contain lists of lists of nodes having the same father (and having same rank)
        allBrotherList = []
        #@hashBrotherList contains the (n,m) pairs such as if hashBrotherList[i] = (n,m), it means the (i+1)th node in sortedNodesList belongs to the (m+1)th subsublist of the (n+1)th sublist of allBrotherList list 
        hashBrotherList = [ None ] * nodesNumber
        #@hashFatherList contains the (n,m) pairs such as if hashBrotherList[i] = (n,m), it means the (i+1)th node is the father of the nodes in the (m+1)th subsublist of the (n+1)th sublist of allBrotherList list 
        hashFatherList = [ None ] * nodesNumber
        currN = -1
        for rank in ranks:
            currN += 1
            #@brotherNodesList will contain lists of nodes having the same father
            brotherNodesList = []
            currM = -1
            while (currNodeIdent < nodesNumber) and pathsNodes[currNodeIdent] and (sortedNodesList[currNodeIdent][1] == rank):
                #Father of the current node will be the last node in the lineage
                currFather = pathsNodes[currNodeIdent][-1] 
                currFatherIdent = None
                for i in range(nodesNumber):
                    if (sortedNodesList[i] == currFather):
                        currFatherIdent = i
                currM += 1
                #Contains siblings nodes (their ident) and in first position the father
                thisRankNodesList = [ currFatherIdent, currNodeIdent ]
                hashBrotherList[currNodeIdent] = (currN,currM)
                currNodeIdent += 1
                hashFatherList[currFatherIdent] = (currN,currM)
                #Because root have no path
                if (currNodeIdent < pathsNodesLength) and pathsNodes[currNodeIdent]:
                    father = pathsNodes[currNodeIdent][-1]
                #If the rank of node indexed by currNodeIdent+1 is different from that of node indexed by currNodeIdent, the structure of taxonomic tree will ensure the two won't have the same father
                while (father == currFather) and (currNodeIdent < nodesNumber):
                    thisRankNodesList.append(currNodeIdent)
                    if (currNodeIdent < nodesNumber):
                        hashBrotherList[currNodeIdent] = (currN,currM)
                    currNodeIdent += 1
                    #Because root have no path (and no father)
                    if (currNodeIdent < pathsNodesLength) and pathsNodes[currNodeIdent]:
                        father = pathsNodes[currNodeIdent][-1]
                brotherNodesList.append(thisRankNodesList)
            if brotherNodesList:
                allBrotherList.append(brotherNodesList)
        return sortedNodesList,pathsNodes,sampleHitListNodes,allBrotherList,hashBrotherList,hashFatherList,pathsNodesLength,nodesNumber
    #
    #
    #Starting from the leaves:
    #(loop1) For every rank r (S to R)
    #(1) Get all the nodes identifiers of same rank r
    #(loop2) For every node identifier n of this list:
    #(2) if n has no child (we can check this through @hashFatherList: if hashFatherList[n] is to None, then the corresponding node sortedNodesList[n] does not own any child): then construct its tree and store it in @constructedTree (@constructedTree[i] matches node @sortedNodesList[i])
    #(3) else: get n's children identifiers through @hashFatherList[n], get the children's trees, construct the tree associated to n, and store it in @constructedTree[n]
    #(4) Repeat loop2 until rank R and return tree associated to Root
    def addNodeAux(self,paths,sortedNodesList,pathsNodes,samplesHitListNodes,allBrotherList,hashBrotherList,hashFatherList,pathsNodesLength,nodesNumber,ranks):
        #Initializing @constructedTree
        constructedTree = [ None ] * nodesNumber
        sLs = sortedNodesList[:: -1]
        ident = 0
        #loop1
        for r in ranks:
            #Step 1:
            sameRankedNodes = []
            name,rank = sLs.pop()
            while sLs and (rank == r):
                sameRankedNodes.append((name,rank,ident))
                ident += 1
                if sLs:
                    name,rank = sLs.pop()
            #Last node to be popped  has not the same rank as other
            #(if rank != r, then it is obvious, and if sLs is empty, it means root has been popped)
            sLs.append((name,rank))
            if (r == "R"):
                #Root is the only node of rank R
                sameRankedNodes.append((name,rank,ident))
            #loop2:
            while sameRankedNodes:
                nm,rk,idt = sameRankedNodes.pop()
                #Step 3:
                if hashFatherList[idt]:
                    n,m = hashFatherList[idt]
                    #First element is the father
                    #Get children identifiers
                    children = allBrotherList[n][m][1:]
                    childrenTrees = [ constructedTree[x] for x in children ]
                    constructedTree[idt] = TaxoTree(nm,rk,idt,samplesHitListNodes[idt],pathsNodes[idt],childrenTrees,paths)
                #Step 2:
                else:
                    constructedTree[idt] = TaxoTree(nm,rk,idt,samplesHitListNodes[idt],pathsNodes[idt],[],paths)
        constructedTree[-1].children = [ child for child in constructedTree[-1].children if child ]
        return constructedTree[-1]
    #
    #
    def addNode(self,paths,nodesList,samplesList,ranks=["S","G","F","O","C","P","K","R"]):
        from time import time
        start = time()
        sortedNodesList,pathsNodes,sampleHitListNodes,allBrotherList,hashBrotherList,hashFatherList,pathsNodesLength,nodesNumber = self.addNodePreProcess(paths,nodesList,samplesList,ranks)
        finalTree = self.addNodeAux(paths,sortedNodesList,pathsNodes,sampleHitListNodes,allBrotherList,hashBrotherList,hashFatherList,pathsNodesLength,nodesNumber,ranks)
        end = time()
        print "TIME:",(end-start),"sec"
        return finalTree
    #
    #
    def iterate(self,f,DFS=False):
        f(self)
        nodeList = []
        for ch in self.children:
            nodeList.append(ch)
        while nodeList:
            child = nodeList.pop()
            f(child)
            if DFS:
                nodeList = child.children + nodeList
            else:
                nodeList = nodeList + child.children

def printTree(tree,sampleOn=False):
    print "ROOT"
    if not tree:
        #Should not happen...
        print "\n/!\ ERROR: Tree turned out to be None."
        raise ValueError
    print "(%s,%s,%d)"%(tree.name,tree.rank,tree.ident)
    if sampleOn:
        print "SAMPLE HIT LIST"
        print tree.sampleHitList
    nodeList = []
    for ch in tree.children:
        nodeList.append(ch)
    n = len(nodeList)
    print "CHILDREN"
    name = tree.name
    rank = tree.rank
    for i in range(n):
        print "**** child #%d of %s %s"%(i,name,rank)
        tree = nodeList[i]
        printTree(tree,sampleOn)
        if tree:
            print "END #%d %s %s of %s %s ****"%(i,tree.name,tree.rank,name,rank)
    print "---"
    return 0
    
