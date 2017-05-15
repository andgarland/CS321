#Find a Highest-Scoring Peptide in a Proteome against a Spectrum:
#ba11f

#parse Rosalind's input
def parseInput():
    f = open("rosalind_ba11f.txt", "r")

    spec = f.readline().strip()
    prote = f.readline().strip()

    f.close

    #return the given spectrum and proteome
    return [spec, prote]

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

#get the weight of a protein string (sum of the values in the MMT corresponding to each aa in the string)
def weight(peptide):
    weight = 0
    for i in range(len(peptide)):
        weight += MMT[peptide[i]]

    return weight

#score a peptide using a spectrum as weight
def scorePeptide(peptide, spectrum):
    #score to be returned
    score = 0
    #intialize an index variable
    index = -1

    #for each amino acid in the peptide
    for i in range(len(peptide)):
        #increment the index by the mass of that peptide
        index += int(MMT[peptide[i]])
        #increment the score by the weight of the spectrum at that index
        score += spectrum[index]

    return score

#find the highest scoring peptide in a proteome given a spectrum
def maxProteome(spectrum, proteome):
    #variables to be overwritten
    bestScore = -1000000000
    bestPeptide = ""

    #iterate through all possible substrings in the proteome
    for i in range(len(proteome)):
        for j in range(1, len(proteome)):

            #get the weight of the current substring
            wgt = int(weight(proteome[i:j+1]))

            #spectrum length should = weight of peptide
            if wgt == len(spectrum):

                #if this substring/peptide has a better score than the current best peptide/score, replace it
                if scorePeptide(proteome[i:j+1], spectrum) > bestScore:
                    bestScore = scorePeptide(proteome[i:j+1], spectrum)
                    bestPeptide = proteome[i:j+1]

            #as incrementing j will always increase the weight, if we've already surpassed the target weight, break
            if wgt > len(spectrum):
                break

    return bestPeptide

inputs = parseInput()
#spec = "0 0 0 4 -2 -3 -1 -7 6 5 3 2 1 9 3 -8 0 3 1 2 1 8"
spec = inputs[0]
spec = spec.split(" ")
spec = [int(i) for i in spec]
#prote = "XZZXZXXXZXZZXZXXZ"
prote = inputs[1]


print(maxProteome(spec, prote))