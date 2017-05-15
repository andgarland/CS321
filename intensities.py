#Modified version of my SGRA algorithm to score using intensities

#Monoisotopic Mass Table
MMT = {
    'A': 71.03711,
    'C': 103.00919,
    'D': 115.02694,
    'E': 129.04259,
    'F': 147.06841,
    'G': 57.02146,
    'H': 137.05891,
    'I': 113.08406,
    'K': 128.09496,
    'L': 113.08406,
    'M': 131.04049,
    'N': 114.04293,
    'P': 97.05276,
    'Q': 128.05858,
    'R': 156.10111,
    'S': 87.03203,
    'T': 101.04768,
    'V': 99.06841,
    'W': 186.07931,
    'Y': 163.06333,
}

def parseInput():
    #read in the dna strings from a file
    #f = open("rosalind_sgra.txt", "r")
    f = open("tempFull2.txt", "r")
    #f = open("fullTest.txt", "r")
    #f = open("W0200003.txt", "r")
    #f = open("Z119_A3_b_0001.lst", "r")
    #f = open("temp2.txt", "r")
    #f = open("waw.txt", "r")

    weights = []
    intensities = {}

    for line in f:
        temp = line.strip().split(" ")
        weights.append(float(temp[0]))
        intensities[float(temp[0])] = float(temp[1])

    f.close()
    return (weights, intensities)

#return the amino acid who's mass is closest to the difference between two given weights
def identifyProtein(pre, suf):
    diff = suf - pre
    protein = min(MMT, key=lambda y: abs(MMT[y] - diff))

    return protein


#convert a set of weights to a protein string
def identifyProteinString(weights):
    result = ""

    #for each weight
    for i in range(len(weights)-1):
        #compare it to the next weight and find the corresponding amino acid. Add it to the current result string
        result += str(identifyProtein(weights[i], weights[i+1]))

    #return the found protein string
    return result

#if the difference between two weight corresponds to an amino acid in the MMT, return it, else return false
def findProteinMatch(weight, weightTwo):
    for aa in MMT.items():
        if abs(abs(weight - weightTwo) - aa[1]) < tolerance:
            return (aa[0], weightTwo)

    #if no match found
    return False

#create a spectrum graph for a given set of weights
def spectrumGraph(weights):
    #each weight will be a key in the dict
    graph = dict((weight, []) for weight in weights)

    #for each combination of weights
    for weight in weights:
        for target in weights:
            #if the first weight < the target weight
            if target > weight:
                #see if there could be an amino acid between them
                aa = findProteinMatch(weight, target)
                #if so, create an edge from weight to target
                if aa:
                    graph[weight].append(aa)

    return graph

#determine a topological order for the graph
def topOrder(graph, edges):
    nodes = []
    temp = []

    #get all the nodes in the graph
    for node in graph.keys():
        nodes.append((node, 0))
        temp.append((node, 0))

    order = []

    #recursive function visit
    def visit(node):
        #if temporarily marked, return
        if nodes[node][1] == 1:
            return

        #if unmarked, temporarily mark the node
        if nodes[node][1] == 0:
            nodes[node] = (nodes[node][0], 1)

            #visit each of the nodes neighbors
            targets = [edge[1] for edge in graph[nodes[node][0]]]
            for target in filter(lambda x: x[0] in targets, nodes):
                visit(nodes.index(target))

            #set the node to be permanently marked
            nodes[node] = (nodes[node][0], 2)
            #add it to the order
            order.insert(0, nodes[node][0])
            #remove it from the list of unmarked nodes
            temp[node] = None

    #while there are still unmarked nodes
    while not all(node is None for node in temp):
        #visit an unmarked node
        visit(temp.index(next(node for node in temp if node is not None)))

    #return the order
    return order

#find the longest/highest scoring path
def longestPath(graph, intensities):
    dist = dict((key, 0) for key in graph.keys())
    prev = dict((key, 0) for key in graph.keys())

    #get the graph's edges
    edges = []
    for item in graph.items():
        if item[1] != []:
            for edge in item[1]:
                edges.append((item[0], edge[1]))

    #for each node in the topolocgical order
    for node in topOrder(graph, edges):

        #for each outgoing edge in node
        for edge in filter(lambda x: x[0] == node, edges):

            #if the neighbor node's score can be increased, increase it and store that it got its value/came from this node
            if dist[edge[1]] <= dist[node] + intensities[node]:
                dist[edge[1]] = dist[node] + intensities[node]
                prev[edge[1]] = node

    #get the best/highest scoring node
    longest = max(dist, key=dist.get)
    path = []

    #backtrack to find the path
    while prev[longest] != 0:
        path.insert(0, longest)
        longest = prev[longest]

    path.insert(0, longest)

    #return the path
    return path

inputs = parseInput()
weights = inputs[0]
intensities = inputs[1]
tolerance = .8
#for i in range(5):
graph = spectrumGraph(weights)
#paths = findPaths()
path = longestPath(graph, intensities)
protein = identifyProteinString(path)
print(str(tolerance)+ ": " + protein)
#tolerance += .2