#4/6/2019. Calculamos la shannon alpha diversity de los pacientes enfermos y pacientes sanos y sacamos una
#distribucion

import numpy as np
import random as rnd
import copy as cp
import sys
import matplotlib.gridspec as gridspec
import matplotlib
from pylab import *
from matplotlib import colors
import matplotlib.patches as patches

import seaborn as sns
import matplotlib.pyplot as plt
import operator
import time
import collections


#0. FUNCTIONS
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#0.1. Read Data
#===============================================
def read_file(inFileName):
    lines = open(inFileName).readlines()
    d = [line.strip().split() for line in lines]
    return d
#===============================================

#0.2. Shannon Entropy
#===========================================
def Shannon(v):
    Entropy=0
    for i in v:
        if i==0:
            Entropy+=0
        else:
            Entropy+=-i*math.log(i,math.e)

    return Entropy
#===========================================

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#1. PARAMETERS FOR INPUT/OUTPUT DATA
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#1.1. EXTERNAL PARAMETERS
#=====================
Dataset=sys.argv[1]#==
number_of_links=3  #==
rnd.seed(1)        #==
#=====================

#1.2. INTERNAL PARAMETERS
#===================================================================

#::::::::::::
K=10;L=20 #::
#::::::::::::


Patients_per_Dataset={'S-8_Tot':107, 'V-10_Tot':92,\
'V-22_Tot':467, 'V-23_24_Tot':222,'V-25_Tot':883}

Microbes_per_Dataset={'S-8_Tot':128, 'V-10_Tot':137,\
'V-22_Tot':134, 'V-23_24_Tot':118,'V-25_Tot':144}
#===================================================================

#1.3. PATHS FOR FILES
#===============================================================
path_plot='../../Plots_Paper/'
Path_in= '../../Input_Data/' 
#===============================================================


#1.4. ENTROPIES FOR PATIENTS
#=========================================================
Patients=Patients_per_Dataset[Dataset]
Microbes=Microbes_per_Dataset[Dataset]

filename=Dataset + '.txt'
Data=read_file(Path_in + filename)

Entropies=[]
for line in np.transpose(Data):
    line=[float(element) for element in line]
    Entropy_Patient=Shannon(line)
    Entropies.append(Entropy_Patient)

#=========================================================

#1.5. PLOT DISTRIBUTIONS
#===============================================================
sns.distplot(Entropies,bins=10,kde=True)

plt.xlabel('Shannon ' + r'$\alpha$ diversity') 
plt.ylabel('Density')

plt.savefig(path_plot +'Shannon_Diversity_' +Dataset +'_Shannon_Diversity.pdf',dpi=300)
plt.savefig(path_plot +'Shannon_Diversity_' +Dataset +'_Shannon_Diversity.png',dpi=300)

plt.show()
#===============================================================



    
