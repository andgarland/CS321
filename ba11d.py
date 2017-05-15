#Convert Peptide Vector into Peptide
#ba11d

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

#convert a peptide vector to a peptide
def vectorToPeptide(vector):
    #turn the string into an array
    vector = vector.split(" ")

    #convert all non-zero values in the vector to the correct weight
    weights = [int(i)+1 for i, x in enumerate(vector) if x != "0"]
    weights.insert(0, 0)

    #return the protein string corresponding to the set of weights
    return identifyProteinString(sorted(weights))

given = "0 0 0 4 -2 -3 -1 -7 6 5 3 2 1 9 3 -8 0 3 1 2 1 0"
print(vectorToPeptide(given))