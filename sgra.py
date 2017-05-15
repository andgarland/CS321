#Using the Spectrum Graph to Infer Peptides
#Rosalind Problem: 'sgra'
#Given: A list L (of length at most 100) containing positive real numbers.
#Return: The longest protein string that matches the spectrum graph of L (if multiple solutions exist, you may output any one of them). Consult the monoisotopic mass table.

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
    f = open("tempFull.txt", "r")
    #f = open("fullTest.txt", "r")
    #f = open("W0200003.txt", "r")
    #f = open("Z119_A3_b_0001.lst", "r")
    #f = open("temp2.txt", "r")

    weights = []

    for line in f:
        weights.append(float(line.strip()))

    #for specific file
    weights.pop(0)

    f.close()
    return weights

def identifyProtein(pre, suf):
    diff = suf - pre
    protein = min(MMT, key=lambda y: abs(MMT[y] - diff))

    return protein

def identifyProteinString(weights):
    result = ""

    for i in range(len(weights)-1):
        result += str(identifyProtein(weights[i], weights[i+1]))

    return result

def findProteinMatch(weight, weightTwo):
    for aa in MMT.items():
        if abs(abs(weight - weightTwo) - aa[1]) < tolerance:
            return (aa[0], weightTwo)

    #if no match found
    return False

def spectrumGraph(weights):
    graph = dict((weight, []) for weight in weights)

    for weight in weights:
        for target in weights:
            if target > weight:
                aa = findProteinMatch(weight, target)
                if aa:
                    graph[weight].append(aa)

    return graph

def topOrder(graph, edges):
    nodes = []
    temp = []
    for node in graph.keys():
        nodes.append((node, 0))
        temp.append((node, 0))

    order = []

    def visit(node):
        if nodes[node][1] == 1:
            return

        if nodes[node][1] == 0:
            nodes[node] = (nodes[node][0], 1)

            targets = [edge[1] for edge in graph[nodes[node][0]]]
            for target in filter(lambda x: x[0] in targets, nodes):
                visit(nodes.index(target))

            nodes[node] = (nodes[node][0], 2)
            order.insert(0, nodes[node][0])
            #del temp[node]
            temp[node] = None

    while not all(node is None for node in temp):
        visit(temp.index(next(node for node in temp if node is not None)))

    return order

def longestPath(graph):
    dist = dict((key, 0) for key in graph.keys())
    prev = dict((key, 0) for key in graph.keys())

    edges = []
    for item in graph.items():
        if item[1] != []:
            for edge in item[1]:
                edges.append((item[0], edge[1]))

    for node in topOrder(graph, edges):
        for edge in filter(lambda x: x[0] == node, edges):
            if dist[edge[1]] <= dist[node] + 1:
                dist[edge[1]] = dist[node] + 1
                prev[edge[1]] = node

    longest = max(dist, key=dist.get)
    path = []

    while prev[longest] != 0:
        path.insert(0, longest)
        longest = prev[longest]

    path.insert(0, longest)

    return path

weights = parseInput()
tolerance = .1
for i in range(5):
    graph = spectrumGraph(weights)
    path = longestPath(graph)
    protein = identifyProteinString(path)
    print(str(tolerance)+ ": " + protein)
    tolerance += .2