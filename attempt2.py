#!/usr/bin/env python3

import numpy as np
import pandas as pd
import scipy
from scipy import stats
import pickle

#inputs

table = "steinmetz_paProb.txt" 
design = "Steinmetz_design.txt"
suffix = "paProb" #names files how you want
detail = ""  #if name in design file differs from column titles in table, write the addendem here
genetics = ["cft2-1", "pcf11-2", "ipa1-1", "wt"] #the different groups
Genetics = ["cft2_1", "pcf11_2", "ipa1_1", "wt"] #the different groups, replacing - with _ for file naming purposes
control = "wt" #name the control, t-tests are run against it

dg = pd.read_csv(table, delimiter = "\t") #creates a pandas table from the text file for table
df = pd.read_csv(design, delimiter = "\t") #creates a pandas table from the text tile for design

gene_array = dg["#gene"].tolist()      	
extra_array = dg["position"].tolist()    #print any added information, such as position on gene
genetic_array = df["genetics"].tolist()  #the differentiator, such as genetics, as it appears as a column title in the design file
sample_array = df["sample"].tolist()

sample_names = []
t_test_array = []

#this loop finds which samples belong to each genetic type, and puts the list of names into a .data file containing a list of samples for the type. They are stored using pickle
for s in range(0, len(genetics)):
  for i in range(0, len(genetic_array)):
    if genetic_array[i] == genetics[s]:
      sample_names.append(sample_array[i])
  with open(Genetics[s]+'.data', 'wb') as fp:
    pickle.dump(sample_names, fp)
  del sample_names[:]


#reads the pickle file, and combines the arrays for each sample of a given genetic type, 


for i in range(0, len(Genetics)) :
  with open(Genetics[i]+'.data', 'rb') as fp4:
    x  = []
    x.append(pickle.load(fp4))
    listitem = x
    bigarray = []
    for s in range(0, len(listitem[0])):
      array = dg[listitem[0][s]+detail].tolist()
      with open(genetics[i]+listitem[0][s]+'_single_array.data', 'wb') as fp1:
        pickle.dump(array, fp1)
      with open(suffix+'_means.txt', 'a+') as f:
        print(listitem[0][s], "	", np.mean(array), file=f)
      bigarray  = np.hstack((bigarray,array))
      del array[:]
    with open(suffix+'_means.txt', 'a+') as f:
      print("	", file=f)
    with open(suffix+'_means.txt', 'a+') as f:
      print(genetics[i], "mean - ", np.mean(bigarray), file=f)
    with open(suffix+'_means.txt', 'a+') as f:
      print("	", file=f)
    with open(genetics[i]+'_combined.data', 'wb') as fp2:
      pickle.dump(bigarray, fp2)  				

fp.close()
fp1.close()
fp2.close()

with open(control+'_combined.data', 'rb') as fp:
  Control_array = pickle.load(fp)


for a in range(0, len(Genetics)) :
   if Genetics[a] == control:
     print("nope")
   else:
     with open(genetics[a]+'_combined.data', 'rb') as fp3:
       f = []
       f.append(pickle.load(fp3))
       answer = scipy.stats.ttest_ind(Control_array, f[0], equal_var = False)
       with open('T-Test_'+suffix+detail+'.txt', 'a+') as f:
           print(genetics[a], "vs WT unequal variance t-test result        ", answer, file=f)

for n in range(0, len(gene_array)) : 
  t_test_array = []
  control_array = []
  for i in range(0, len(Genetics)) :
    with open(Genetics[i]+'.data', 'rb') as fp5:   #this line is broken
      h = []
      h.append(pickle.load(fp5))
      listitem = h
      if genetics[i] == control:
        control_genetic_array = []
        for s in range(0, len(listitem[0])):
          with open(genetics[i]+listitem[0][s]+'_single_array.data', 'rb') as fp1:
            Control = pickle.load(fp1)
            control_genetic_array.append(Control[n])
        control_array.append(control_genetic_array)
      else:
        gene_genetic_array = []
        for s in range(0, len(listitem[0])):
          with open(genetics[i]+listitem[0][s]+'_single_array.data', 'rb') as fp1:
            array = pickle.load(fp1)
            gene_genetic_array.append(array[n])
        t_test_array.append(gene_genetic_array)
  for a in range(0, len(t_test_array)) :
     genetics.remove(control)
     answer = scipy.stats.ttest_ind(control_array[0], t_test_array[a], equal_var = False)
     with open('Gene_based_T-Test_'+suffix+detail+'.txt', 'a+') as f:
        print(gene_array[n], extra_array[n], genetics[a], "vs wt unequal variance t-test result        ", answer, file=f)  
     genetics.append(control)

