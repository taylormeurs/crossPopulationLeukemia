snpDict = {}
with open("gwas-association-downloaded_2021-03-03-EFO_0000094-withChildTraits.tsv") as snp_file:
    snp_file.readline()
    for line in snp_file:
        row = line.rstrip().split(",")
        if row[11] not in snpDict.keys():
            snpDict[row[11]] = [row[21]]
        else:
            snpDict[row[11]].append(row[21])

with open("selectedSNPs", "w") as write_file:
    for chrmNum in snpDict.keys():
        print("looking at chromosome " + str(chrmNum))
        chrmSNPList = snpDict[chrmNum]
        with open("ALL.chr" + str(chrmNum) + ".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz") as chrmFile:
            line = chrmFile.readline()
            while line[0] != "#":
                line = chrmFile.readline()
            while line != "":
                row = line.rstrip().split(",")
                if row[2] in chrmSNPList:
                    write_file.write(line)
                line = chrmFile.readline()
