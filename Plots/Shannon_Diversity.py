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
Diseases=0         #==
#=====================

#1.2. INTERNAL PARAMETERS
#===================================================================

#::::::::::::::::::::
K=10;L=20

if "L13" in Dataset:
    L=13
#::::::::::::::::::::


Patients_per_Dataset=\
{'S-8_Tot':107, 'V-10_Tot':92, 'V-22_Tot':467, 'V-23_24_Tot':222,'V-25_Tot':883,\
 'He_CD':47, 'He_Healthy':39, 'Nielsen_IBD':78, 'Nielsen_Healthy':58, \
 'Vogtmann_CRC': 47 , 'Vogtmann_Healthy':30 , 'Zeller_CRC': 135 , 'Zeller_Healthy':42 ,\
 'Zhang_RA': 92,'Zhang_Healthy':55,\
 'He_CD_L13':47, 'He_Healthy_L13':39, 'Nielsen_IBD_L13':78, 'Nielsen_Healthy_L13':58, \
 'Vogtmann_CRC_L13': 47 , 'Vogtmann_Healthy_L13':30 , 'Zeller_CRC_L13': 135 , 'Zeller_Healthy_L13':42 ,\
 'Zhang_RA_L13': 92,'Zhang_Healthy_L13':55, 'Test':5}


Microbes_per_Dataset=\
{'S-8_Tot':128, 'V-10_Tot':137, 'V-22_Tot':134, 'V-23_24_Tot':118,'V-25_Tot':144,\
 'He_CD':313, 'He_Healthy':313, 'Nielsen_IBD':313 , 'Nielsen_Healthy':313,\
 'Vogtmann_CRC':313 , 'Vogtmann_Healthy':313 , 'Zeller_CRC':313 , 'Zeller_Healthy':313 ,\
 'Zhang_RA':313,'Zhang_Healthy':313,\
 'He_CD_L13':313, 'He_Healthy_L13':313, 'Nielsen_IBD_L13':313 , 'Nielsen_Healthy_L13':313,\
 'Vogtmann_CRC_L13':313 , 'Vogtmann_Healthy_L13':313 , 'Zeller_CRC_L13':313 , 'Zeller_Healthy_L13':313 ,\
 'Zhang_RA_L13':313,'Zhang_Healthy_L13':313, 'Test':5}
#===================================================================

                       
if Diseases==1:

    #1.3. PATHS FOR FILES
    #===============================================================
    path_plot='../../Plots_Paper/'
    #===============================================================

    #1.4. HEALTHY
    #===============================================================
    Dataset_Name=Dataset.split('_')[0]
    Patients_Healthy=Patients_per_Dataset[Dataset_Name + '_Healthy']
    Microbes_Healthy=Microbes_per_Dataset[Dataset_Name + '_Healthy']
    Path_in= '../../Input_Data/' + 'Datasets_Diseases/' + 'raw/'
    
    filename_healthy=Dataset_Name + '_Healthy.txt'
    Data_Healthy=read_file(Path_in + filename_healthy)

    Entropies_H=[]
    for line in np.transpose(Data_Healthy):
        line=[float(element) for element in line]
        Entropy_Patient_H=Shannon(line)
        Entropies_H.append(Entropy_Patient_H)
    #===============================================================

    #1.5. DISEASED
    #===============================================================
    Patients_Diseased=Patients_per_Dataset[Dataset]
    Microbes_Diseased=Microbes_per_Dataset[Dataset]
    
    filename_diseased=Dataset.split('_')[0] +\
                       '_' + Dataset.split('_')[1] + '.txt'
    Data_Diseased=read_file(Path_in + filename_diseased)

    Entropies_D=[]
    for line in np.transpose(Data_Diseased):
        line=[float(element) for element in line]
        Entropy_Patient_D=Shannon(line)
        Entropies_D.append(Entropy_Patient_D)
    #===============================================================

    #1.6. PLOT DISTRIBUTIONS
    #===============================================================
    sns.distplot(Entropies_H,bins=10, color='r',kde=True)
    sns.distplot(Entropies_D,bins=10,kde=True)
    
    plt.xlabel('Shannon ' + r'$\alpha$ diversity') 
    plt.ylabel('Density')
    plt.title( Dataset )
    plt.legend(['Healthy','Diseased'])

    plt.savefig(path_plot +'Shannon_Diversity_'  + Dataset +'_Shannon_Diversity.png',dpi=300)
    
    plt.show()
    #===============================================================

else:

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



    
