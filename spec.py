#Inferring Protein from Spectrum
#Rosalind Problem: 'spec'
#Given: A list L of n (n≤100) positive real numbers.
#Return: A protein string of length n−1 whose prefix spectrum is equal to L (if multiple solutions exist, you may output any one of them). Consult the monoisotopic mass table.

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

#parse the initial file
def parseInput():
    #read in the data from a file
    f = open("rosalind_spec.txt", "r")
    #f = open("temp.txt", "r")

    weights = []

    for line in f:
        weights.append(float(line.strip()))

    f.close()

    return weights

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

weights = parseInput()
print(identifyProteinString(weights))
