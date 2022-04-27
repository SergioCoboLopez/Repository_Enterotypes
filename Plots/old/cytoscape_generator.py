#23/01/2019. Este codigo genera los ficheros de datos para pintar las meta-redes de grupos de microbios y de
#grupos de personas en cytoscape.

import sys
import numpy as np
import csv

#0. FUNCTIONS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#0.1. Read Data
#===============================================
def read_file(inFileName):
    lines = open(inFileName).readlines()
    d = [line.strip().split() for line in lines]
    return d
#===============================================


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



#1. PARAMETERS FOR INPUT/OUTPUT DATA
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#1.1. EXTERNAL PARAMETERS
#=====================
Dataset=sys.argv[1]#==
number_of_links=3  #==
#=====================

#1.2. INTERNAL PARAMETERS
#==================================================================
K=10;L=20

Patients_per_Dataset=\
{'S-8_Tot':107, 'V-10_Tot':92, 'V-22_Tot':467, 'V-23_24_Tot':222,\
 'V-25_Tot':883}

Microbes_per_Dataset=\
{'S-8_Tot':128, 'V-10_Tot':137, 'V-22_Tot':134, 'V-23_24_Tot':118,\
'V-25_Tot':144}

Patients=Patients_per_Dataset[Dataset]
Microbes=Microbes_per_Dataset[Dataset]
#==================================================================

#1.3. PATHS FOR FILES
#==============================================================
Path_out='../../Output_Data/'+ 'Leave_One_Out_' + Dataset + '/'
Path_in= '../../Input_Data/' + 'Leave_One_Out_' + Dataset + '/'
#==============================================================

#1.4. CHOOSE BEST LIKELIHOOD
#======================================================================
Likelihood_Seed={}
for Seed in range(1,10,2):
    likelihood_file=\
    'Dataset_' + Dataset + '_Seed_' + str(Seed) + '_LogLikelihood.txt'
    likelihoodfinal=read_file(Path_out + likelihood_file)[-1]
    Likelihood_Seed[float(likelihoodfinal[0])]=Seed

    print Likelihood_Seed
    print max(Likelihood_Seed)
    print Likelihood_Seed[max(Likelihood_Seed)]
    Seed=Likelihood_Seed[max(Likelihood_Seed)]

file='Dataset_' + Dataset + '_Seed_' + str(Seed) + '_Parameters.txt'
Data=read_file(Path_out+file)
#======================================================================

#1.6. COMPUTE THETAS, ETAS AND AGGREGATED
#=================================================
theta=[];contador_theta=0
eta=[];contador_eta=0

for line in Data[:Patients]:
    line.insert(0, contador_theta )
    contador_theta+=1

theta=Data[:Patients]

for line in Data[Patients+1:Patients+1+Microbes]:
    line.insert(0, contador_eta )
    contador_eta+=1

eta=Data[Patients+1:Patients+1+Microbes]

theta_aggregate=[0]*K
for lineT in theta:
    for i in range(K):
        theta_aggregate[i]+=float(lineT[i+1])

eta_aggregate=[0]*L
for lineE in eta:
    for j in range(L):
        eta_aggregate[j]+=float(lineE[j+1])
#================================================= 



#1.5. READ PMATRICES AND FORMAT INTO MATRICES
#======================================================================
counter=0
with open('cytoscape_' + str(Dataset) + '.csv','w') as network: 
    for line in Data:
        if 'slice' and '0]' in line:
            print line
            for row in range(K):
                for col in range(L):
                    Real_Link=1 - float(Data[counter+row+1][col])
                    # if Real_Link==0:
                    #     continue
                    # else:
                    network.write('Patients_'+str(row) + ' ' + \
                                      'Microbes_'+str(col) + ' ' + \
                                      str(Real_Link)+' ' + str(theta_aggregate[row]) + ' '+\
                                      str(2.5*eta_aggregate[col]) +'\n' )
    
        counter+=1
#======================================================================

print theta_aggregate, '\n'
print eta_aggregate
