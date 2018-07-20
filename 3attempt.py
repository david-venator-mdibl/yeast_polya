#!/usr/bin/env python3

import numpy as np
import pandas as pd
import scipy
from scipy import stats
import pickle

#inputs

table = "steinmetz_paGene.txt" 
design = "Steinmetz_design.txt"
suffix = "paGene" #names files how you want
detail = "_AvgL3pUTR"  #if name in design file differs from column titles in table, write the addendem here
genetics = ["cft2-1", "pcf11-2", "ipa1-1", "wt"] #the different groups
Genetics = ["cft2_1", "pcf11_2", "ipa1_1", "wt"] #the different groups, replacing - with _ for file naming purposes
#_genetics_ = ["Cft2_1", "Pcf11_2", "Ipa1_1". "Wt"]
control = "wt" #name the control, t-tests are run against it

dg = pd.read_csv(table, delimiter = "\t") #creates a pandas table from the text file for table
df = pd.read_csv(design, delimiter = "\t") #creates a pandas table from the text tile for design

gene_array = dg["#gene"].tolist()      	
extra_array = dg["SampleTotal"].tolist()    #print any added information, such as position on gene
genetic_array = df["genetics"].tolist()  #the differentiator, such as genetics, as it appears as a column title in the design file
sample_array = df["sample"].tolist()

sample_names = []
t_test_array = []

Associative_Array = {}

#this loop finds which samples belong to each genetic type, and puts the list of names into a .data file containing a list of samples for the type. They are stored using pickle
for s in range(0, len(genetics)):
  sample_names = []
  for i in range(0, len(genetic_array)):
    if genetic_array[i] == genetics[s]:
      sample_names.append(sample_array[i])
  Associative_Array[Genetics[s]] = sample_names
  with open(Genetics[s]+'.data', 'wb') as fp:
    pickle.dump(sample_names, fp)
  #del sample_names[:]

print(Associative_Array)

print(Associative_Array[Genetics[0]])

#reads the pickle file, and combines the arrays for each sample of a given genetic type, 

AArray2 = {}

AArray3 = {}

for i in range(0, len(genetics)) :
   mutation = Associative_Array[Genetics[i]]
   bigarray = []
   array = []
   print(mutation)
   for s in range(0, len(mutation)) :
     array = dg[mutation[s]+detail].tolist()
     AArray2[mutation[s]] = array
     with open(suffix+'_meansy.txt', 'a+') as f:
       print(mutation[s], "	", np.mean(array), file=f)
     bigarray  = np.hstack((bigarray,array))
   with open(suffix+'_meansy.txt', 'a+') as f:
     print("	", file=f)
   with open(suffix+'_meansy.txt', 'a+') as f:
     print(genetics[i], "mean - ", np.mean(bigarray), file=f)
   with open(suffix+'_meansy.txt', 'a+') as f:
     print("	", file=f)
   AArray3[genetics[i]] = bigarray
   
print(AArray2)
print(AArray3)
	
print(Genetics[0])

Control_array = AArray3[control]

for a in range(0, len(Genetics)) :
   if Genetics[a] == control:
     print("nope")
   else:
     mutation = AArray3[genetics[a]] 
     answer = scipy.stats.ttest_ind(Control_array, mutation, equal_var = False)
     with open('new_T-Test_'+suffix+detail+'.txt', 'a+') as f:
         print(genetics[a], "vs WT unequal variance t-test result        ", answer, file=f)

with open('GGene_based_T-Test_'+suffix+detail+'.txt', 'a+') as f:
  print("Gene", "additional", "name", "p_value", file = f) 

for n in range(0, len(gene_array)) : 
  t_test_array = []
  control_array = []
  for i in range(0, len(Genetics)) :
    type = Associative_Array[Genetics[i]]   
    if genetics[i] == control:
      control_genetic_array = []
      for s in range(0, len(type)):
        Control = AArray2[type[s]]
        control_genetic_array.append(Control[n])
      control_array.append(control_genetic_array)
    else:
      gene_genetic_array = []
      for s in range(0, len(type)):
        array = AArray2[type[s]]
        gene_genetic_array.append(array[n])
      t_test_array.append(gene_genetic_array)
  for a in range(0, len(t_test_array)) :
     genetics.remove(control)
     answer = scipy.stats.ttest_ind(control_array[0], t_test_array[a], equal_var = False)
     with open('GGene_based_T-Test_'+suffix+detail+'.txt', 'a+') as f:
        print(gene_array[n], extra_array[n], genetics[a], answer[1], file=f)  
     genetics.append(control)

dv = pd.read_csv('GGene_based_T-Test_'+suffix+detail+'.txt', delimiter = " ") #creates a pandas table from the text file for table

signif = dv.query('p_value < 0.05')

with open('PPandas_pvalues.txt', 'w') as f:
  print(signif, file=f)


