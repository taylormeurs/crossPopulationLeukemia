import sys

openFile = sys.argv[1]
finalFile = sys.argv[2]

tissues = []
with open("tissue_list.txt", "r") as tissues_file:
    for line in tissues_file:
        tissues.append(line.strip())

genes = set()
badGenes = []
lineDict = {}
popGenes = {}

for tissue in tissues:
    lineDict[tissue] = {}
    popGenes[tissue] = []
with open(openFile) as genes_fa:
    header = genes_fa.readline().rstrip()
    line = genes_fa.readline()
    while line != "":
        row = line.rstrip().split(',')
        genes.add(row[1])
        if row[1] not in lineDict[row[2]].keys():
            lineDict[row[2]][row[1]] = line.rstrip()
        # handles genes with two isoforms
        # same expression = keep the both
        # different expression = throw out both isoforms
        else:
            # print(row[2] + " " + row[1])
            badGeneFirstCopy = lineDict[row[2]][row[1]].split(',')
            if badGeneFirstCopy[7] != row[7]:
                popGenes[row[2]].append(row[1])
                badGenes.append(badGeneFirstCopy[1])
                # badGenes.append((badGeneFirstCopy[0], row[0]))
            else:
                print("Isoform " + row[0] + " " + row[1] + " is a double with equal expression")
        line = genes_fa.readline()
print(set(badGenes))

with open(openFile) as firstFile, open(finalFile,"w") as finalfile:
    for line in openFile:
        row = line.rstrip().split(",")
        if row[1] not in popGenes[row[2]]:
            finalfile.write(line)
