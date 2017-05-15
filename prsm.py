#Matching a Spectrum to a Protein
#Rosalind Problem: 'prsm'
#Given: A positive integer n followed by a collection of n protein strings s1, s2, ..., sn and a multiset R of positive numbers (corresponding to the complete spectrum of some unknown protein string).
#Return: The maximum multiplicity of R⊖S[sk] taken over all strings sk, followed by the string sk for which this maximum multiplicity occurs (you may output any such value if multiple solutions exist).

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

#Parse Rosalind's Input
def parseInput():
    #read in the dna strings from a file
    f = open("rosalind_prsm.txt", "r")
    #f = open("temp.txt", "r")

    n = int(f.readline().strip())
    weights = []
    strings = []

    #collect n protein strings
    for i in range(n):
        strings.append(f.readline().strip())

    #collect multiset of weights
    for line in f:
        weights.append(float(line.strip()))

    f.close()
    return [strings, weights]

#return the "Minkowski Difference", or a set that holds the result of set 1 - set 2 (each item in s1 -  each item in s2. Multiplicities included)
def minkowskiDifference(set1, set2):
    newSet = {}

    #for each combination of numbers between the two sets
    for num1 in set1:
        for num2 in set2:

            #subtract the second number from the first
            diff = round(num1 - num2, 6)

            #if the result hasn't been seen before, add it to the new set. Else increment the count for that result
            if diff in newSet:
                newSet[diff] += 1
            else:
                newSet[diff] = 1

    return newSet

#get the weight of a protein string (sum of the values in the MMT corresponding to each aa in the string)
def getStringWeights(protein):
    total = 0
    for i in range(len(protein)):
        total += MMT[protein[i]]

    return total

#get all the prefixes and subfixes of a protein string
def getAllFixes(protein):
    fixes = []
    for i in range(1,len(protein)):
        fixes.append(protein[i:])
        fixes.append(protein[:i])

    return fixes

#return the set of all the prefix and suffix string weights of a protein string
def proteinToWeightSet(protein):
    wSet = []

    #get all pre/suffixes
    fixes = getAllFixes(protein)

    #add all the prefix/suffix weights to the set of weights
    for fix in fixes:
        wSet.append(getStringWeights(fix))

    return wSet

#find the protein with the largest number of shared masses to an unknown protein (weights). Return the found protein and the maximum shared multiplicity
def match(proteins, weights):
    weightSets = []

    #for each protein, collect it's weight set
    for protein in proteins:
        weightSets.append(proteinToWeightSet(protein))

    maxMult = -1
    maxProt = ""

    #for each weight set
    for i in range(len(weightSets)):

        #get the minkowski difference of it and the reference set
        diff = minkowskiDifference(weights, weightSets[i])

        #find the largest multiplicity in the set
        mult = max(diff.values())

        #if it is higher than the current highest multiplicity, replace it and save the protein that corresponded to the weight set
        if mult > maxMult:
            maxMult = mult
            maxProt = proteins[i]

    return (maxMult, maxProt)

inputs = parseInput()
results = match(inputs[0], inputs[1])
print(results[0])
print(results[1])