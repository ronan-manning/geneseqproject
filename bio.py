dros = open("Drosophila_proteome.fasta (2)", 'r')

singList = ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"]
doubList = []
tripList = []

"""Rewrite with no line breaks so extra duos/triples are counted"""

for first in singList:
    for second in singList:
        doubList.append(first + second)

for duo in doubList:
    for char in singList:
        tripList.append(duo + char)

freqDict = {} # freqDict contains the number of occurrences of each single character, double, and pair

for char in singList:
    freqDict[char] = 0
for char in doubList:
    freqDict[char] = 0
for char in tripList:
    freqDict[char] = 0


for line in dros:
    if line[0] != ">":
        for char in line:
            if char in singList:
                single = char
                freqDict[single] = freqDict[single] + 1
            if (line.index(char) + 2 <= len(line)) and char in singList and line[line.index(char) + 1] in singList:   # Tests if there is room for a pair starting with this character
                pair = line[line.index(char):line.index(char) + 2]
                freqDict[pair] = freqDict[pair] + 1
            if (line.index(char) + 3 <= len(line)) and char in singList and line[line.index(char) + 1] in singList and line[line.index(char) + 2] in singList:   # Tests if there is room for a trio starting with this character
                trio = line[line.index(char):line.index(char) + 3]
                freqDict[trio] = freqDict[trio] + 1

# Creating xFreqDoubs
singleRateDict = {}
totalSingles = 540911+142129+382375+475895+265254+446730+195129+365958+415090+678742+173789+347244+407241+382103+408981+618306+423587+437051+74309+221876

for single in singList:
    singleRateDict[single] = freqDict[single] / totalSingles  # singleRateDict contains each single character with the proportion of the total characters that are it

xFreqDoubs = {}

for pair in doubList:
    one = singleRateDict[pair[0]]
    two = singleRateDict[pair[1]]
    xFreqDoubs[pair] = one * two

doubleRateDict = {} # doubleRateDict contains each pair with the proportion of the total pairs that are it

for pair in doubList:
    doubleRateDict[pair] = freqDict[pair] / (totalSingles - 1)  # The -1 should be minus the number of rows, but this will be irrelevant when pairs and trios across lines are fixed

doubleRateOverExpected = {} # Maybe change this from raw occurrence rate minus expected occurrence rate to percent difference?

for pair in doubList:
    doubleRateOverExpected[pair] = doubleRateDict[pair] - xFreqDoubs[pair]

#Creating xFreqTrios
xFreqTrios = {}

for trio in tripList:
    one = singleRateDict[trio[0]]
    two = singleRateDict[trio[1]]
    three = singleRateDict[trio[2]]
    xFreqTrios[trio] = one * two * three

tripleRateDict = {}

for trio in tripList:
    tripleRateDict[trio] = freqDict[trio] / (totalSingles - 1)

tripleRateOverExpected = {}

for trio in tripList:
    tripleRateOverExpected[trio] = tripleRateDict[trio] - xFreqTrios[trio]

print(tripleRateOverExpected)
