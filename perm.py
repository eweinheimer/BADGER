import sys
import re
import math
import itertools
import statistics
import numpy as np
import scipy.stats as stats
import argparse 

parser = argparse.ArgumentParser(description='Bivariate Analysis of Differential Gene Expression Reactions')
parser.add_argument('-c', '--control', required=True, help='control condition as indicated in experimental design file')
parser.add_argument('-f', '--factors', required=True, help='CSV file of experimental design factors')
parser.add_argument('-i', '--input', required=True, help='CSV file of normalized read counts')
args = parser.parse_args()
control = args.control
factors = args.factors
input = args.input

def processExpDesign(file):
    with open(file) as design_file:
        header = False
        coldict = {}
        times = []
        treatments = []
        groups = []
        for line in design_file:
            if header:
                newline = line.split("\n")[0].split(",")
                times.append(newline[coldict["Time"]])
                treatments.append(newline[coldict["Treatment"]])
                groups.append([newline[coldict["Time"]],newline[coldict["Treatment"]],newline[coldict["Individual"]],newline[coldict["SampleID"]]])
            else:
                header = line.split(",")
                for i in range(len(header)):
                    if re.search("time", header[i], re.IGNORECASE):
                        coldict["Time"] = i
                    elif re.search("treat", header[i], re.IGNORECASE):
                        coldict["Treatment"] = i
                    elif re.search("sample", header[i], re.IGNORECASE):
                        coldict["SampleID"] = i
                    elif re.search("indiv", header[i], re.IGNORECASE):
                        coldict["Individual"] = i
    return groups,np.unique(times),np.unique(treatments)

def permutationTest(file):
    timecomps = [[i,j] for (i,j) in itertools.combinations(times,2) if i<j]
    treatcomps = [[i,j] for (i,j) in itertools.combinations(treats,2) if i!=j]
    for i in timecomps:
        perms = False
        for j in treatcomps:
            indivs = []
            if j[0]==control:
                controllabel=j[0]
                treatmentlabel=j[1]
            else:
                controllabel=j[1]
                treatmentlabel=j[0]
            for k in range(len(groups)):
                if (groups[k][1]==controllabel or groups[k][1]==treatmentlabel) and groups[k][2] not in indivs:
                    indivs.append(groups[k][2])    
            perms = list(itertools.combinations(list(range(0,len(indivs))),math.ceil(len(indivs)/2)))
            controlids_t0 = {}
            controlids_t1 = {}
            treatids_t0 = {}
            treatids_t1 = {}
            
            for k in range(len(groups)):
                if groups[k][0]==i[0] and groups[k][1]==controllabel:
                    controlids_t0[groups[k][2]] = groups[k][3]
                elif groups[k][0]==i[1] and groups[k][1]==controllabel:
                    controlids_t1[groups[k][2]] = groups[k][3]
                elif groups[k][0]==i[0] and groups[k][1]==treatmentlabel:
                    treatids_t0[groups[k][2]] = groups[k][3]
                elif groups[k][0]==i[1] and groups[k][1]==treatmentlabel:
                    treatids_t1[groups[k][2]] = groups[k][3]
            
            with open(file) as infile, open("permutation_output_comp{}.csv".format(i[0] + "v" + i[1] + "." + j[0] + "v" + j[1]),"w") as f:
                f.write("Gene" + "," + "PropSigPermutations" + "," + "CritPVal" + "," + "Mean_{}".format(treatmentlabel) + "-Mean_{}".format(controllabel) + "," + "MaxIndivExpression" + "," + "SDD" + "\n")
                samples = False
                for line in infile:
                    if samples:
                        sigreps = 0
                        controldiffs = []
                        treatdiffs = []
                        countssub = []
                        for k in controlids_t0:
                            controldiffs.append(float(line.split(",")[samples.index(controlids_t1[k])]) - float(line.split(",")[samples.index(controlids_t0[k])]))
                            countssub.append(float(line.split(",")[samples.index(controlids_t1[k])]))
                            countssub.append(float(line.split(",")[samples.index(controlids_t0[k])]))
                        for k in treatids_t0:
                            treatdiffs.append(float(line.split(",")[samples.index(treatids_t1[k])]) - float(line.split(",")[samples.index(treatids_t0[k])]))
                            countssub.append(float(line.split(",")[samples.index(treatids_t1[k])]))
                            countssub.append(float(line.split(",")[samples.index(treatids_t0[k])]))
                        if len(indivs) % 2 == 0:
                            reps = perms[0:(int(len(perms)/2))]
                        else:
                            reps = perms
                        merge = treatdiffs + controldiffs
                        critTstat, critPval = stats.ttest_ind(controldiffs,treatdiffs)
                        for k in range(len(reps)):
                            exDiff1 = list()
                            exDiff2 = list()
                            order = [int(x) for x in reps[k]]
                            recip = [x for x in list(range(0,len(indivs))) if x not in order]
                            exDiff1 = [[merge[m]] for m in order]
                            exDiff2 = [[merge[m]] for m in recip]
                            permTstat, permPval = stats.ttest_ind(exDiff1,exDiff2)
                            if permPval < critPval:
                                sigreps+=1
                        if max(countssub) == 0:
                            entry = [line.split(",")[0],sigreps/len(reps),critPval,0,0,0]
                        else:
                            entry = [line.split(",")[0],sigreps/len(reps),critPval,(statistics.mean(treatdiffs)-statistics.mean(controldiffs)),max(countssub),(statistics.mean(treatdiffs)-statistics.mean(controldiffs))/max(countssub)]
                        for m in range(len(entry)):
                            f.write(str(entry[m]) + ",") if m < len(entry)-1 else f.write(str(entry[m]) + "\n")
                    else:
                        samples = line.split("\n")[0].split(",")
    return ''

groups,times,treats = processExpDesign(factors)
permutationTest(input)


