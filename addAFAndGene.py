import sys
import numpy as np
from statsmodels.stats.proportion import proportions_ztest

alphaSub = 0.05/(365*325)
alphaSup = 0.05/3650
subpopulations = ['GBR', 'FIN', 'CHS', 'PUR', 'CDX', 'CLM', 'IBS', 'PEL', 'PJL', 'KHV', 'ACB',
                  'GWD', 'ESN', 'BEB', 'MSL', 'STU', 'ITU', 'CEU', 'YRI', 'CHB', 'JPT', 'LWK',
                  'ASW', 'MXL', 'TSI', 'GIH']
superpopulations = ['America', 'East Asia', 'South Asia', 'Africa', 'Europe']


def TukeyToCSV(d, csvName):
    header = ['SNP', 'risk allele', 'population1', 'population2', 'gene', 'allele frequency 1', 'total samples 1',
              'allele frequency 2', 'total samples 2', 'p_value\n']
    with open(csvName, 'w') as tukey_csv_file:
        tukey_csv_file.write(",".join(header))
        for SNP, nucDict in d.items():
            for nuc, p_val_dict in nucDict.items():
                if len(p_val_dict.keys()) > 0:
                    for p_val, popsList in p_val_dict.items():
                        for pops in popsList:
                            if nuc in closestGenes[SNP].keys():
                                this_gene = closestGenes[SNP][nuc]
                            else:
                                this_gene = "?"
                            af0 = sum(AFdict[SNP][nuc][pops[0]])
                            tot0 = len(AFdict[SNP][nuc][pops[0]])
                            af1 = sum(AFdict[SNP][nuc][pops[1]])
                            tot1 = len(AFdict[SNP][nuc][pops[1]])
                            this_gene = ";".join(this_gene)
                            strToWrite = ",".join([SNP, nuc, pops[0], pops[1], this_gene, str(af0), str(tot0), str(af1), str(tot1), str(p_val) + "\n"])
                            tukey_csv_file.write(strToWrite)


def t_test(curDict, sub):
    td = {}
    if len(curDict['Africa']) != 0:
        if sub:
            subpop = [pop for pop in subpopulations]
            for pop1 in subpopulations:
                subpop.pop(0)
                for pop2 in subpop:
                    count = np.asarray([sum(curDict[pop1]), sum(curDict[pop2])])
                    nobs = np.asarray([len(curDict[pop1]), len(curDict[pop2])])
                    stat, P = proportions_ztest(count, nobs)
                    if P < alphaSub:
                        if P not in td.keys():
                            td[P] = [[pop1, pop2]]
                        else:
                            td[P].append([pop1, pop2])
        else:
            suppop = [pop for pop in superpopulations]
            for pop1 in superpopulations:
                suppop.pop(0)
                for pop2 in suppop:
                    count = np.asarray([sum(curDict[pop1]), sum(curDict[pop2])])
                    nobs = np.asarray([len(curDict[pop1]), len(curDict[pop2])])
                    stat, P = proportions_ztest(count, nobs)
                    if P < alphaSup:
                        if P not in td.keys():
                            td[P] = [[pop1, pop2]]
                        else:
                            td[P].append([pop1, pop2])
    return td


# make population map
subpopMap = {}
superpopMap = {}
with open('sampleKeyPop.csv') as popKey_file:
    for line in popKey_file:
        row = line.rstrip().split(",")
        subpopMap[row[0]] = row[1]
        superpopMap[row[0]] = row[2]

# get list of leukemia SNPs
allSNPs = {}
with open("leukSNPs") as SNP_file:
    SNP_file.readline()
    for line in SNP_file:
        row = line.rstrip().split(",")
        if row[0] not in allSNPs.keys():
            allSNPs[row[0]] = set(row[1])
        else:
            allSNPs[row[0]].add(row[1])

closestGenes = {}
with open("closestGene") as gene_file:
    gene_file.readline()
    for line in gene_file:
        row = line.rstrip().split(",")
        rs = row[1].split("-")
        allele = rs[1]
        rs = rs[0]
        if rs not in closestGenes.keys():
            closestGenes[rs] = {}
            closestGenes[rs][allele] = [row[0]]
        elif allele not in closestGenes[rs].keys():
            closestGenes[rs][allele] = [row[0]]
        else:
            closestGenes[rs][allele].append(row[0])

num_tests = 0

subTukeyResults = {}
superTukeyResults = {}
runForSub = True

AFdict = {}

# fill dictionaries
for SNP in allSNPs.keys():
    print(SNP)
    SNP_tab = "\t" + SNP + "\t"
    found = False
    # reinitialize nucDicts
    nucleotideDicts = {}
    if "?" in allSNPs[SNP]:
        nucs = allSNPs[SNP].union({'A', 'T', 'C', 'G'})
        nucs.remove("?")
    else:
        nucs = allSNPs[SNP]
    for nucleotide in nucs:
        nucleotideDicts[nucleotide] = {}
        for superpopulation in superpopulations:
            nucleotideDicts[nucleotide][superpopulation] = []
        for subpopulation in subpopulations:
            nucleotideDicts[nucleotide][subpopulation] = []
    with open(sys.argv[1]) as csv_file:  # arg[1] is the selected data from 1000 genomes
        head = csv_file.readline().rstrip().split("\t")
        for line in csv_file:
            if line.find(SNP_tab) != -1:
                row = line.rstrip().split("\t")
                for nucleotide in nucs:
                    # define allele type
                    if row[3] == nucleotide:
                        # alleleType = "ref"
                        alleleMatch = 0
                    elif nucleotide in row[4].split(","):
                        # alleleType = "alt"
                        alleleMatch = row[4].split(",").index(nucleotide) + 1
                    else:
                        # alleleType = "NA"
                        alleleMatch = -1

                    # add all sample data to population lists if we have allele
                    i = 9
                    if alleleMatch != -1:
                        found = True
                        num_tests += 1
                        for sample in head[9:]:
                            allele1 = int(int(row[i][0]) == alleleMatch)
                            allele2 = int(int(row[i][-1]) == alleleMatch)
                            nucleotideDicts[nucleotide][subpopMap[sample]] += [allele1, allele2]
                            nucleotideDicts[nucleotide][superpopMap[sample]] += [allele1, allele2]
                            i += 1
    if found:
        subTukeyResults[SNP] = {}
        superTukeyResults[SNP] = {}
        for nuc in nucs:
            if len(nucleotideDicts[nuc]["Europe"]) > 0:
                if SNP not in AFdict.keys():
                    AFdict[SNP] = {}
                    AFdict[SNP][nuc] = nucleotideDicts[nuc]
                else:
                    AFdict[SNP][nuc] = nucleotideDicts[nuc]
                subTukeyResults[SNP][nuc] = t_test(nucleotideDicts[nuc], runForSub)
                superTukeyResults[SNP][nuc] = t_test(nucleotideDicts[nuc], not runForSub)

print("number of tests: " + str(num_tests))
#print(AFdict)
TukeyToCSV(superTukeyResults, 'superPropTestAF.csv')
TukeyToCSV(subTukeyResults, 'subPropTestAF.csv')
