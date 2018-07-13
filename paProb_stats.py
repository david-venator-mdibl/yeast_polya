#!/usr/bin/env python

import numpy as np
import pandas as pd
import scipy
from scipy import stats

INPuT = "steinmetz_paProb.txt"

suffix = "paProb"

df = pd.read_csv(INPuT, delimiter = "\t") #creates a pandas table from the text file

with open('pandas_table_'+suffix+'.txt', 'w') as f: #prints the df file, not necessary but helpful for visualation
   print(df, file=f)

column_titles = list(df.columns.values) #creates a list of column titles based off if the pandas file

array = np.genfromtxt(INPuT, skip_header=1, unpack=True) #creates an array for each column of the text file and and defines them by outputdata

for column in range(1, len(column_titles)):              #loops through, averaging each array and finding each column titles, and listening them together
   Mean = np.mean(array[column])
   title = column_titles[column] 
   with open('averages_'+suffix+'.txt', 'a+') as f:
      print(title, "	",  Mean, file=f)

with open('averages_'+suffix+'.txt', 'a+') as f:   #provides spacer in output document
   print("	", file=f)
   print("Combined_Averages_UTR", file=f)
   print("	", file=f)

WT_a  = np.hstack((array[10],array[11],array[12],array[20]))  #combines sample arrays to one array for each genetic type

ipa1_1_a = np.hstack((array[13],array[14],array[21]))

pcf11_2_a = np.hstack((array[15],array[16],array[22],array[23]))

cft2_1_a = np.hstack((array[11],array[12],array[13],array[18],array[19]))
    
for genetic in range(0, 4, 1):                                     #generates and lables averages for each genetic makeup
   mutation = (WT_a, ipa1_1_a, pcf11_2_a, cft2_1_a)
   mutation_name = ("WT_a", "ipa1_1_a", "pcf11_2_a", "cft2_1_a")
   Strain_Mean = np.mean(mutation[genetic])
   with open('averages_'+suffix+'.txt', 'a+') as f:
     print(mutation_name[genetic], "	", Strain_Mean, file=f)
   
for genetic in range(0, 3, 1):                                      #does t-test for each mutation vs WT
   mutation = (ipa1_1_a, pcf11_2_a, cft2_1_a)
   mutation_name = ("ipa1_1_a", "pcf11_2_a", "cft2_1_a")
   answer = scipy.stats.ttest_ind(WT_a, mutation[genetic], equal_var = False)
   with open('T-Test_'+suffix+'.txt', 'a+') as f:
     print(mutation_name[genetic], "vs WT unequal variance t-test result	", answer, file=f)
