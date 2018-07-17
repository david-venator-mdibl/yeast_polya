#!/usr/bin/env python3

import numpy as np
import pandas as pd
import scipy
from scipy import stats
import pickle

INPuT = "steinmetz_paProb.txt"

suffix = "paProb"

detail = ""

control = "wt"

dg = pd.read_csv(INPuT, delimiter = "\t") #creates a pandas table from the text file

gene_array = dg["#gene"].tolist()

inpUt = "Steinmetz_design.txt"

#creates a pandas file from design file
df = pd.read_csv(inpUt, delimiter = "\t")

genetic_array = df["genetics"].tolist()

sample_array = df["sample"].tolist()

genetics = ["cft2-1", "pcf11-2", "ipa1-1", "wt"]

#remove dashes, python doesnt like the - sign
Genetics = ["cft2_1", "pcf11_2", "ipa1_1", "wt"]

sample_names = []

#this loop finds which samples belong to each genetic type, and puts the list of names into a .data file containing a list of samples for the type. They are stored using pickle
for s in range(0, len(genetics)):
  for i in range(0, len(genetic_array)):
    if genetic_array[i] == genetics[s]:
      sample_names.append(sample_array[i])
  print(sample_names)
  with open(Genetics[s]+'.data', 'wb') as fp:
    pickle.dump(sample_names, fp)
  del sample_names[:]

#this loop is supposed to combine use the sample titles found 
#in the first loop to create one bug array for each genetic 
#type - which worked, based on the print(bigarray) and np.mean(bigarray) working). 
#The issue was that the big array will not save into a .data file, even though
#im using essentially the same code from the previous loop.

for i in range(0, len(Genetics)) :
  with open(Genetics[i]+'.data', 'rb') as fp:
    Genetics[i] = []
    Genetics[i].append(pickle.load(fp))
    listitem = Genetics[i]
    print(Genetics[i])
    bigarray = []
    for s in range(0, len(listitem[0])):
      array = dg[listitem[0][s]+detail].tolist()
      with open(genetics[i]+listitem[0][s]+'_single_array.data', 'wb') as fp:
        pickle.dump(array, fp)
      print(np.mean(array))
      bigarray  = np.hstack((bigarray,array))
      del array[:]
    print(genetics[i], "mean - ", np.mean(bigarray))
    with open(genetics[i]+'_combined.data', 'wb') as fp:
      pickle.dump(bigarray, fp)  				

for n in range(0, len(gene_array)) :
  for i in range(0, len(Genetics)) :
    with open(Genetics[i]+'.data', 'rb') as fp:   #this line is broken
      Genetics[i] = []
      Genetics[i].append(pickle.load(fp))
      listitem = Genetics[i]
      gene_genetic_array = []
      for s in range(0, len(listitem[0])):
        with open(genetics[i]+listitem[0][s]+'_single_array.data', 'rb') as fp:
          array = pickle.load(fp)
          gene_genetic_array.append(array[n])
          with open(gene_array[n]+'_'+genetics[i], 'wb') as fp:
            pickle.dump(gene_genetic_array, fp)
      #del Genetics[i][:]

#with open(ipa1-1lane3TS1248II_30nt_Ttrim_single_array.data, 'rb') as fp:
  #blah = pickle.load(fp)

#meant to assign a variable to the control array, once that works
with open(control+'_combined.data', 'rb') as fp:
  control_array = pickle.load(fp)

#this is the t-test function, meant to compare each of the mutations to the wildtype. doesnt work yet because of the issues storing the bigarray in the previous loop.
for a in range(0, len(Genetics)) :
   if Genetics[a] == control:
     stop
   else:
     with open(genetics[a]+'_combined.data', 'rb') as fp:
       Genetics[a] = []
       Genetics[a].append(pickle.load(fp))
       answer = scipy.stats.ttest_ind(control_array, Genetics[a][0], equal_var = False)
       with open('T-Test_'+suffix+detail+'.txt', 'a+') as f:
           print(genetics[a], "vs WT unequal variance t-test result        ", answer, file=f)


