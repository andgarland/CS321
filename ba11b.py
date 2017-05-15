#Implement DecodingIdealSpectrum
#ba11b

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

#Parse Rosalind's input
def parseInput():
    f = open("rosalind_ba11b.txt", "r")
    #f = open("temp.txt", "r")

    weights = f.readline().strip().split(" ")
    #turn the input into an array of ints
    weights = [int(weight) for weight in weights]

    f.close()
    return weights

#if the difference between two weight corresponds to an amino acid in the MMT, return it, else return false
def findProteinMatch(weight, weightTwo):
    for aa in MMT.items():
        if abs(abs(weight - weightTwo) - int(aa[1])) < 0.001:
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

#print the spectrum graph in a format accepted by Rosalind
def printSpectrumGraph(graph):
    keys = sorted(graph.keys())

    for key in keys:
        if graph[key] != []:
            for prot in graph[key]:
                print(str(key) + "->" + str(prot[1]) + ":" + str(prot[0]))

#recursively find all paths in a graph from a source node to the destination
def findPaths(graph, src, dest, path = []):
    #create a new path with the source, or append the current source node to the existing current path
    path = path + [src]

    #if we're at the dest, return the path
    if src[1] == dest:
        return [path]

    #if we're at a node not in the graph, return an empty array
    if src[1] not in graph.keys():
        return []

    paths = []

    #for each of the current node's outgoing edges
    for node in graph[src[1]]:
        #if the node isn't already in our path
        if node not in path:
            #recursively collect the new group of possible paths the follow this node
            newPaths = findPaths(graph, node, dest, path)

            #add all of the newfound paths to the current collection of paths
            for newPath in newPaths:
                paths.append(newPath)

    #return all the paths found
    return paths

#get the weight of a protein string (sum of the values in the MMT corresponding to each aa in the string)
def getStringWeights(protein):
    total = 0
    for i in range(len(protein)):
        total += int(MMT[protein[i]])

    return total

#get all prefixes and suffixes of a protein
def getAllFixes(protein):
    fixes = [protein]
    for i in range(1,len(protein)):
        fixes.append(protein[i:])
        fixes.append(protein[:i])

    return fixes

#collect the entire set of protein weights
def proteinToWeightSet(protein):
    wSet = []

    #get all prefixes and suffixes
    fixes = getAllFixes(protein)

    #for each pre/suffix
    for fix in fixes:

        #find the weight of the string and append it to the set
        wSet.append(getStringWeights(fix))

    return wSet

#find the ideal peptide based on a given spectrum
def decodingIdealSpectrum(spectrum):

    #create the spectrum graph
    graph = spectrumGraph(spectrum)

    #find all paths from 0 to the last element in the spectrum
    paths = findPaths(graph, ("start", 0), spectrum[-1], [])
    spectrum.pop(0)

    #for each found path
    for path in paths:

        #assemble the peptide that corresponds to the path
        peptide = ""
        for node in path[1:]:
            peptide += node[0]

        #generate the spectrum of the found protein
        spec = sorted(proteinToWeightSet(peptide))

        #if it matches the given spectrum, retrun the peptide
        if spec == spectrum:
            return peptide

weights = parseInput()
weights.insert(0, 0)
print(decodingIdealSpectrum(weights))