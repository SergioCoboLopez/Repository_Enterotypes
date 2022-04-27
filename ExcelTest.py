import pandas as pd
import numpy as np
import os

#0.FUNCTIONS
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++

#0.1. Read Data
#===============================================
def read_file(inFileName):
    lines = open(inFileName).readlines()
    d = [line.strip().split() for line in lines]
    return d
#===============================================

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++

path='../Input_Data/Leave_One_Out_'
Datasets=['S-8', 'V-10', 'V-22', 'V-23_24', 'V-25']

tuples=[]
for Dataset in Datasets:
    folder=path+Dataset+'/'
    print folder
    Patients=read_file(folder+'Patient_Record_'+Dataset+'.txt')
    #----------------------------------------------------
    for line in Patients:
        Patient=line[0]
    #----------------------------------------------------
    	for Fold in range(5):

            for measure in ["Likelihood", "Accuracy"]:
                tuples.append((str(Dataset), str(Patient),str(Fold), measure))

print tuples


index = pd.MultiIndex.from_tuples(tuples, names=\
                ['Dataset', 'Patient','Fold','Magnitude'])

df = pd.DataFrame(index=index, columns=list('13579'))
print df.loc['S-8', '14','1', 'Accuracy']['1']
df.to_csv('../Reports/Results_Sheet.csv')

