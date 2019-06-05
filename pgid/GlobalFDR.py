#!/usr/bin/env python
import sys
import os
import getopt

'''
sort the input file according to SpectEvalue in assending order first
sort -s -g -t$'\t' -k13,13 file
the script uses  T-TDC method in following paper to calculate FDR
Improved False Discovery Rate Estimation Procedure for Shotgun Proteomics.
Uri Keich, Attila Kertesz-Farkas,and William Stafford Noble. 2015 JPR

'''

decoy_prefix = "XXX_"
psm_qval = 0.01

if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
    print("Warning! wrong command!")
    print("Example: python GlobalFDR.py --input PSM_filename --output output_filename --decoy_prefix XXX_  --psm_qval 0.01")
else:
    options, remainder = getopt.getopt(sys.argv[1:],'', ['input=','output=','decoy_prefix=','psm_qval='])
    for opt, arg in options:
        if opt == '--input': input_file=arg
        elif opt == '--output': output_file=arg
        elif opt == '--decoy_prefix': decoy_prefix=arg
        elif opt == '--psm_qval': psm_qval = float(arg)
        else:
            print("Warning! Command-line argument: %s not recognized. Exiting..." % opt); sys.exit()



input=open(input_file,'r')
output=open(output_file,'w')

header=input.readline().strip().split("\t")
header=header+["PSM-FDR","Pep-FDR"]
output.write("\t".join(header)+"\n")

targetcount=0
decoycount=0

pep_dic={}

for line in input:
    row=line.strip().split('\t')
    pep=row[8]
    acc=row[9]

    if decoy_prefix in acc:
        decoycount+=1
    else:
        targetcount+=1

    FDR=float(decoycount)/float(targetcount)
 
    if pep not in pep_dic:
       pep_dic[pep] = FDR
    
    row.append(str(FDR)) # PSM-FDR
    row.append(str(pep_dic[pep])) #Pep-FDR

    if FDR<psm_qval:
        if decoy_prefix not in acc: #write only target PSMs
            output.write("\t".join(row)+"\n")
    else:
        break;

print "target hits",targetcount
print "decoy hits", decoycount

input.close()
output.close()
