#!/usr/bin/env python3

import numpy as np
import pandas as pd
import scipy
from scipy import stats
import statsmodels.api as sm

#inputs

table = "pcf11_paGene.txt"
design = "New_Pcf11_DRS_design.txt"
suffix = "pcf11_paGene" #names files how you want
detail = "_AvgL3pUTR"  #if name in design file differs from column titles in table, write the addendem here
columns_want = ["rgenetics", "media"]
tidbit = "rgenetics_"

db = pd.read_csv(table, delimiter = "\t") #creates a pandas table from the text file for table
df = pd.read_csv(design, delimiter = "\t") #creates a pandas table from the text tile for design

sample_array = df["sample"].tolist()
genetic_array = df["rgenetics"].tolist()
batch_array = df["media"].tolist()

Genetics = []
Batch = []

for n in range(len(sample_array)) :
  if genetic_array[n] not in Genetics :
    Genetics.append(genetic_array[n])
  if batch_array[n] not in Batch :
    Batch.append(batch_array[n])

print(Genetics)
print(Batch)

column_titles = list(df.columns.values)

dg = db.query('UTR3count > 1')

gene_array = dg["#gene"].tolist()
UTR_array = dg["UTR3count"].tolist()

for n in range(0, len(column_titles)) :
  if column_titles[n] not in columns_want :
    df.drop([column_titles[n]], axis=1, inplace=True)
  else :
    continue

du = pd.get_dummies(df)

column_titles = list(du.columns.values)

du["samples"] = sample_array

print(du)

print(column_titles)
  
with open(tidbit+"OutPutTable.txt", "a+") as f:
  print("#gene", "3UTRcount", "rsquared", "rsquared_adj", "fvalue", "f_pvalue", file=f)

column_titles.insert(0, "const")
with open(tidbit+"OutPutTable2.txt", "a+") as f:
  for k in range(len(column_titles)):	
    print(column_titles[k], end =" ", file=f)
    print(column_titles[k]+"_stderror", end =" ", file=f)
    print(column_titles[k]+"_p_value", end =" ", file=f) 
  print("", file=f)
column_titles.remove("const")

for n in range(0, len(gene_array)) :
  arrayy = []
  justincase = {}
  work = {}
  Values = []
  for s in range(len(sample_array)) :
     array = dg[sample_array[s]+detail].tolist()
     value = array[n]
     Values.append(value)
  work['UTR'] = Values
  dy = pd.DataFrame(data=work)
  X = du[column_titles]
  X = sm.add_constant(X)
  y = dy['UTR']
  model = sm.OLS(y.astype(float), X).fit() 
  predictions = model.predict(X)
  with open(tidbit+"OutPutReg.txt", "a+") as f:
    print(gene_array[n], file=f)
    print(model.summary(), file=f)
  with open(tidbit+"OutPutTable.txt", "a+") as f:
    print(gene_array[n], UTR_array[n], model.rsquared, model.rsquared_adj, model.fvalue, model.f_pvalue, file=f)
  column_titles.insert(0, "const")
  with open(tidbit+"OutPutTable2.txt", "a+") as f:
    for i in range(len(column_titles)) :
      print(model.params[i], end =" ", file=f)
      print(model.bse[i], end =" ", file=f)
      print(model.pvalues[i], end =" ", file=f)
    print("", file=f)
  column_titles.remove("const")

dn = pd.read_csv(tidbit+"OutPutTable2.txt", sep = " ")
  
dq = pd.read_csv(tidbit+"OutPutTable.txt", sep = " ")

result = pd.concat([dq, dn], axis=1)

with open(tidbit+"pandas_version.txt", "w") as f:
  print(result, file=f)

result.to_csv(tidbit+"tab_delimited_result2.txt", header = True, index = True, sep = "\t")
