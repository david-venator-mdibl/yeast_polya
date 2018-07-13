#!/usr/bin/env python

import numpy as np
import pandas as pd
import scipy
from scipy import stats

inpUt = "Steinmetz_design.txt"

df = pd.read_csv(inpUt, delimiter = "\t")

genetic_array = df["genetics"].tolist()

sample_array = df["sample"].tolist()

#for genetic in (pcf11-2, cft2-1, ipa1-1)

genetics = ["cft2-1", "pcf11-2", "ipa1-1", "wt"]

sample_names = []

sample_names_all = []

for s in range(0, len(genetics)):
  for i in range(0, len(genetic_array)):
    if genetic_array[i] == genetics[s]:
      sample_names.append(sample_array[i])
  print(sample_names)
  with open('array_of_arrays.txt', 'a+') as f:
    print(sample_names, file=f)
  del sample_names[:]


#for number in range(0, 4):
  #for sample in range(0, len(array_of_arrays.txt[number])):
    #array1 = df["sample"].tolist()
  #print(array1)


#print(sample_names)

#print(genetic_array)

#print(sample_array)


