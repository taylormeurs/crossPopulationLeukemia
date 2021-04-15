totalPops = 0
snpDict = {}
with open(r"C:\Users\taylo\Documents\Winter 2021\Genetics of Disease\subTTestCorrected.csv") as snp_file:
    snp_file.readline()
    for line in snp_file:
        row = line.rstrip().split(",")
        totalPops += 2
        if row[2] + " and " + row[3] in snpDict.keys():
            snpDict[row[2] + " and " + row[3]] += 1
        elif row[3] + " and " + row[2] in snpDict.keys():
            snpDict[row[3] + " and " + row[2]] += 1
        else:
            snpDict[row[2] + " and " + row[3]] = 1

with open(r"C:\Users\taylo\Documents\Winter 2021\Genetics of Disease\pairSubPopCounts.csv", "w") as write_file:
    write_file.write(",".join(["Pop", "Count\n"]))
    for pop, count in snpDict.items():
        write_file.write(",".join([pop, str(count) + "\n"]))
