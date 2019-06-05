#!/usr/bin/env python3

import sys
import os
import getopt
import numpy as np
import numpy.polynomial.polynomial as poly


'''
sort the file according to SpecEvalue first
sort -s -g -t$'\t' -k13,13 file

the script used the method described in the following paper to estimate subgroup specific FDR. 
Transferred subgroup false discovery rate for rare post-translational modifications detected by mass spectrometry.
Mol Cell Proteomics. 2014 May

'''

psm_qval = 0.1
if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
    print("Warning! wrong command!")
    print("Example: python subgroupFDR.py --input PSM_filename --output output_filename --decoy_prefix XXX_ --group_target SAAV --group_decoy XXX_SAAV")
else:
    options, remainder = getopt.getopt(sys.argv[1:],'', ['input=','output=','decoy_prefix=','group_target=','group_decoy=','psm_qval='])
    for opt, arg in options:
        if opt == '--input': input_file=arg
        elif opt == '--output': output_file=arg
        elif opt == '--decoy_prefix': decoy_prefix=arg
        elif opt == '--group_target': group_target=arg
        elif opt == '--group_decoy': group_decoy=arg
        elif opt == '--psm_qval': psm_qval = float(arg)
        else:
            print("Warning! Command-line argument: %s not recognized. Exiting..." % opt); sys.exit()


input=open(input_file,'r')# the input tsv file need to be sorted by SpectEval in assending order first
output=open(output_file,'w')

header = input.readline().strip().split("\t")
header += ["group-PSMFDR","group-PepFDR"]

output.write("\t".join(header)+"\n")

score = []
score_dic = {}

group_targetcount=0
group_decoycount=0

decoycount=0
targetcount=0

pep_dic={}
grouppep_dic={}

for line in input:
    row=line.strip().split('\t')
    pro=row[9]
    specEval = -np.log10(float(row[12]))
    
    if decoy_prefix in pro:
        decoycount+=1
        if group_decoy in pro:
            group_decoycount+=1
    else:
        targetcount+=1
        if group_target in pro:
            group_targetcount+=1
 
    score_dic[specEval] = [targetcount,decoycount,group_targetcount,group_decoycount]

input.close()

x=[]
y=[]

for score in score_dic:
     counts=score_dic[score]
     if counts[3]>0: # if group_decoycount>0
         y.append(counts[3]/counts[1])
         x.append(score)

coefs = poly.polyfit(np.array(x),np.array(y), 1)

intercept = coefs[0]
slope = coefs[1]

input2=open(input_file,'r')

pep_col = header.index("Peptide")
prot_col = header.index("Protein")

for line in input2:
    row=line.strip().split('\t')
    pep=row[pep_col]
    pro=row[prot_col]
    if group_target not in pro:
        continue;
    
    specEval = -np.log10(float(row[12]))
    counts = score_dic[specEval]   

    targetcount=float(counts[0])
    decoycount=float(counts[1])
    FDR=decoycount/targetcount

    group_targetcount=float(counts[2])
    gamma = slope*specEval+intercept
    groupFDR = FDR*gamma*(targetcount/group_targetcount)

    if pep not in grouppep_dic:
        grouppep_dic[pep]= groupFDR
   
    row.append(str(groupFDR))
    row.append(str(grouppep_dic[pep]))
   
    if groupFDR < psm_qval:
        if decoy_prefix not in pro: #write only target PSMs
            output.write("\t".join(row)+"\n")

print ("Hits from the subgroup: targe,decoy",group_targetcount,group_decoycount)
print ("intercept and slope is:",intercept,slope)

input2.close()
output.close()
