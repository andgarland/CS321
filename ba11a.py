#Construct the Graph of a Spectrum
#ba11a

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

#parse Rosalind's given file
def parseInput():
    f = open("rosalind_ba11a.txt", "r")
    #f = open("temp.txt", "r")

    weights = f.readline().strip().split(" ")
    #convert to an array of ints
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

#code to be run
weights = parseInput()
weights.insert(0, 0)
graph = spectrumGraph(weights)
printSpectrumGraph(graph)