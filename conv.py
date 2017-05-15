#Comparing Spectra with Spectral Convolution
#Rosalind Problem: 'conv'
#Given: Two multisets of positive real numbers S1 and S2. The size of each multiset is at most 200.
#Return: The largest multiplicity of S1⊖S2, as well as the absolute value of the number x maximizing (S1⊖S2)(x) (you may return any such value if multiple solutions exist).

def parseInput():
    f = open("rosalind_conv.txt", "r")
    #f = open("temp.txt", "r")

    one = f.readline().strip().split(" ")
    two = f.readline().strip().split(" ")

    setOne = []
    setTwo = []

    #cast the first set to floats
    for num in one:
        setOne.append(float(num))

    #cast the second set to floats
    for num in two:
        setTwo.append(float(num))

    f.close()
    return [setOne, setTwo]

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

sets = parseInput()
newSet = minkowskiDifference(sets[0], sets[1])
lm = max(newSet, key=newSet.get)
#print the highest multiplicity and it's value
print(newSet[lm])
print(abs(lm))