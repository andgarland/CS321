#Inferring Peptide from Full Spectrum
#Rosalind Problem: 'full'
#Given: A list L containing 2n+3 positive real numbers (n≤100). The first number in L is the parent mass of a peptide P, and all other numbers represent the masses of some b-ions and y-ions of P (in no particular order). You may assume that if the mass of a b-ion is present, then so is that of its complementary y-ion, and vice-versa.
#Return: A protein string t of length n for which there exist two positive real numbers w1 and w2 such that for every prefix p and suffix s of t, each of w(p)+w1 and w(s)+w2 is equal to an element of L. (In other words, there exists a protein string whose t-prefix and t-suffix weights correspond to the non-parent mass values of L.) If multiple solutions exist, you may output any one.

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
    #f = open("rosalind_full.txt", "r")
    f = open("tempFull.txt", "r")

    weights = []

    for line in f:
        weights.append(float(line.strip()))

    f.close()

    return weights

#given a weight and a set of weights, see if there is a weight in that set who's difference from the first weight corresponds to the mass of an amino acid. Return false if no weight was found
def findProteinMatch(weight, weights):
    for other in weights:
        diff = abs(weight - other)

        #extra check -- no amino acid has a weight greater than 195
        if diff < 195:
            for aa in MMT.items():
                if abs(diff - aa[1]) < 0.8:
                    return (aa[0], other)

    #if no match found
    return False

#assemble a protein string based off of a set of weights
def assembleProtein(weights):
    #it is given that the string will be n length
    n = (len(weights) - 3) / 2
    protein = ""

    #we start with the first weight (weight[0] is the parent weight)
    total = weights[1]
    weights = weights[2:]

    #while our assembled protein is not n length
    while len(protein) < n:
        #find the next possible amino acid
        aa = findProteinMatch(total, weights)

        #if we find a match
        if aa:

            #add it to the protein string
            protein += aa[0]

            #increase the current total mass of our string by its mass
            total += MMT[aa[0]]

            #eliminate any weights smaller than the found weight
            weights = weights[weights.index(aa[1]):]

    return protein

weights = parseInput()
print(weights)
print(assembleProtein(weights))