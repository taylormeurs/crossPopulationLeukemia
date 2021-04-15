import csv

snpDict = {}
with open("superPropTest.csv") as snp_file:
    snp_file.readline()
    for line in snp_file:
        row = line.rstrip().split(",")
        if row[0] not in snpDict.keys():
            snpDict[row[0]] = {}
            snpDict[row[0]][row[1]] = 1
        elif row[1] not in snpDict[row[0]].keys():
            snpDict[row[0]][row[1]] = 1
        else:
            snpDict[row[0]][row[1]] += 1

with open("superSigSNPCounts.csv", "w") as write_file:
    write_file.write(",".join(["SNP", "NUC", "Count\n"]))
    for snp, nucDic in snpDict.items():
        for nuc, num in nucDic.items():
            write_file.write(",".join([snp, nuc, str(num) + "\n"]))
