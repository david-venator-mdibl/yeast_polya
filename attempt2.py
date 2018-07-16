#!/usr/bin/env python3

import numpy as np
import pandas as pd
import scipy
from scipy import stats
import pickle

INPuT = "steinmetz_paProb.txt"

suffix = "paProb"

control = "wt"

dg = pd.read_csv(INPuT, delimiter = "\t") #creates a pandas table from the text file

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
for listitem in Genetics :
  with open(listitem+'.data', 'rb') as fp:
    listitem = []
    listitem.append(pickle.load(fp))
    for s in range(0, len(listitem)):
      array = dg[listitem[0][s]].tolist()
      bigarray = np.hstack((array))
      del array[:]
    print(bigarray)
    print(np.mean(bigarray))
    with open(listitem+'_combined.data', 'wb') as fp:
      pickle.dump(bigarray, fp)  

#meant to assign a variable to the control array, once that works
with open(control+'combined.data', 'rb') as fp:
  control_array = pickle.load(fp)

print(control_array)


#print(np.mean(control_array))

#with open('cft2_1_combined.data', 'rb') as fp:
  #cft2_1_array = pickle.load(fp)

#print(cft2_1_array)

#print(np.mean(cft2_1_array))

#thing = scipy.stats.ttest_ind(control_array, cft2_1_array, equal_var = False)

#print(thing)


#the above hashtagged sutff was to check iff the bigarrays were saving as a .data file


#this is the t-test function, meant to compare each of the mutations to the wildtype. doesnt work yet because of the issues storing the bigarray in the previous loop.
for listitem in Genetics :
   if listitem != control:
     with open(listitem+'_combined.data', 'rb') as fp:
       listitem = []
       listitem.append(pickle.load(fp))
       answer = scipy.stats.ttest_ind(control_array, listitem[0], equal_var = False)
       with open('T-Test_'+suffix+'.txt', 'a+') as f:
           print(str(listitem), "vs WT unequal variance t-test result        ", answer, file=f)


