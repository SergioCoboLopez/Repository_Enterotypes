#13/076/2018. Este codigo hace un heatmap de accuracies para diferentes combinaciones de parametros
#K y L. Toma como entrada el fichero Accuracies.txt

import numpy as np
import pandas as pd
import math
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib


#0. FUNCTIONS
#++++++++++++++++++++++++++++++++++++++++++++++++++

#0.1. Reading data
#================================================
def read_file(inFileName):
    lines = open(inFileName).readlines()
    d = [line.strip().split() for line in lines]
    return d
#================================================

#++++++++++++++++++++++++++++++++++++++++++++++++++


#1. READ DATA
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
path='/export/home/shared/Projects/Microbiome/Output_Data/'
File=path+'Accuracies.txt'
Data=read_file(File)
Data=Data[1:]
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#2. CONVERT DATA TO LISTS
#+++++++++++++++++++++++++++++++++++++++++++
K=[];L=[];BestAccuracy=[]
for line in Data:
    K.append( int(line[0]) )
    L.append( int(line[1]) )
    BestAccuracy.append( float(line[2]) )
#+++++++++++++++++++++++++++++++++++++++++++

#3. PLOTS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#Define areas of the dots proportional to Accuracy
#===========================================================================
area = [ ( 130*number/sum(BestAccuracy) )**3 for number in BestAccuracy]
#===========================================================================

plt.scatter(K, L, c=BestAccuracy , s=area )
plt.xlabel("K (#Groups of Patients)")
plt.ylabel("L (#Groups of Microbes)  ")

#Plot colorbar
#==================================================
cbar= plt.colorbar()
cbar.set_label("Predictive Accuracy", labelpad=+10)
#==================================================
plt.show()


