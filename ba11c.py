#Convert a Peptide into a Peptide Vector
#ba11c

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

    'X': 4,
    'Z': 5,
}

#get all the prefixes of a protein string
def getPrefixes(protein):
    fixes = [protein]
    for i in range(1,len(protein)):
        fixes.append(protein[i:])

    return fixes

#get the weight of a protein string (sum of the values in the MMT corresponding to each aa in the string)
def getStringWeights(protein):
    total = 0
    for i in range(len(protein)):
        total += int(MMT[protein[i]])

    return total

#convert a peptide to a "peptide vector"
def peptideToVector(peptide):
    weights = []

    #find all of the peptide's prefixes
    prefixes = getPrefixes(peptide)

    #collect the weights of all the prefixes
    for pf in prefixes:
        weights.append(getStringWeights(pf))

    #create an empty vector who's length is equal to the weight of the entire peptide
    vector = [0 for i in range(weights[0])]

    #for each prefix weight, change the index corresponding to the value of the weight to 1
    for weight in weights:
        vector[weight-1] = 1

    #return the binary peptide vector
    return vector

vector = peptideToVector("HVIAGEDAWDDPPYIILQSTWEPYWMC")
print(" ".join([str(i) for i in vector]))
