import csv
totalPops = 0
snpDict = {}
with open("subPropTestAF.csv") as snp_file:
    snp_file.readline()
    for line in snp_file:
        row = line.rstrip().split(",")
        totalPops += 2
        if row[2] not in snpDict.keys():
            snpDict[row[2]] = 1
        else:
            snpDict[row[2]] += 1
        if row[3] not in snpDict.keys():
            snpDict[row[3]] = 1
        else:
            snpDict[row[3]] += 1

with open("subPopCounts.csv", "w") as write_file:
    write_file.write(",".join(["Pop", "Count\n"]))
    for pop, count in snpDict.items():
        write_file.write(",".join([pop, str(count) + "\n"]))
