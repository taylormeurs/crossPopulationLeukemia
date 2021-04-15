sigSNPS = set()
sigStudies = set()
with open("superPropTestAF.csv") as snp_file:
    snp_file.readline()
    for line in snp_file:
        row = line.rstrip().split(",")
        sigSNPS.add(row[0])

with open("subPropTestAF.csv") as snp_file:
    snp_file.readline()
    for line in snp_file:
        row = line.rstrip().split(",")
        sigSNPS.add(row[0])

studydict = {}
with open("gwas-association-downloaded_2021-03-03-EFO_0000094-withChildTraits.tsv", "r") as studyFile:
    studyFile.readline()
    for line in studyFile:
        row = line.rstrip().split(",")
        if row[21] not in studydict.keys():
            studydict[row[21]] = [row[6]]
        else:
            studydict[row[21]].append(row[6])
for SNP in sigSNPS:
    for study in studydict[SNP]:
        sigStudies.add(study)

print(len(sigSNPS))
print(len(sigStudies))
