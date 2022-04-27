#17/04/2022. Adaptation and corrections of Save_Intra-Inter_Distances_Hosts.py for metadata with hosts

import sys
import os
import numpy as np
import random as rnd
import pandas as pd
import copy as cp
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
import scipy.cluster
from scipy import stats
from scipy.spatial import distance
import math
import operator
from collections import Counter
import time
from itertools import product
import pickle
import seaborn as sns
start_time = time.time()

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

#0.2. Measure Intra Euclidean distances
#=================================================================================
def Measure_Intra_Distances(Dict1):
    Dict2={}
    #Run over phylum
    for phylum in Dict1:
        Microbes_per_phylum=len(Dict1[phylum])
        #Run over pairs of microbes
        for microbe_a in range(Microbes_per_phylum):
            for microbe_b in range(microbe_a+1,Microbes_per_phylum):
               a=np.array(Dict1[phylum][microbe_a][1])
               b=np.array(Dict1[phylum][microbe_b][1])
               dist = np.linalg.norm(a - b)


               Ids_distance_tuple=(\
               Dict1[phylum][microbe_a][0], Dict1[phylum][microbe_b][0], dist)

               try:
                   Dict2[phylum].append(Ids_distance_tuple)

               except KeyError:
                   Dict2[phylum]=[Ids_distance_tuple]
               
    return Dict2
#=================================================================================

#0.3. Measure Inter Euclidian distances
#=================================================================================
def Measure_Inter_Distances(Dict1):

    Dict_Inter={}    
    Orders=Dict1.keys()

    for order_0 in Dict1:

        Orders.remove(order_0)

    	for order_1 in Orders:
        
            for microbe_0 in Dict1[order_0]:
                for microbe_1 in Dict1[order_1]:

                   Dist_0=np.array(microbe_0[1])
                   Dist_1=np.array(microbe_1[1])
                   Dist = np.linalg.norm(Dist_0 - Dist_1)

                   Ids_distance_tuple=(microbe_0[0], microbe_1[0], Dist)

                   try:
                       Dict_Inter[order_0+'-'+order_1].append(Ids_distance_tuple)

                   except KeyError:
                       Dict_Inter[order_0+'-'+order_1]=[Ids_distance_tuple]
    
    return Dict_Inter
#=================================================================================

#0.4. Extract distances from dictionary
#======================================================================
def Extract_Distances(Dict_Distances,Hosts=0):

    Distances_fun=[]
    Means=[];SEMs=[];stds=[]
    xticks=[]; Sizes_Plot=[]

    for phylum in Dict_Distances:

        if Hosts==0:
            #Extract names of phyla and sizes of groups for plotting purposes
            #----------------------------------------------
            phylum_ticks=phylum.split('__')
            xticks.append(phylum_ticks[1])
            Sizes_Plot.append( len(Dict_Distances[phylum]))
            #----------------------------------------------

        else:
            xticks.append(phylum)
        #----------------------------------------------
        Distances_per_Phylum_fun=[]
        for distances in Dict_Distances[phylum]:
            Distances_per_Phylum_fun.append(distances[2])
        Distances_fun.append(Distances_per_Phylum_fun)
        Mean=sum(Distances_per_Phylum_fun)/float(len(Distances_per_Phylum_fun))
        std=np.std(Distances_per_Phylum_fun)
        SEM=stats.sem(Distances_per_Phylum_fun,ddof=0)
        #----------------------------------------------

        #-------------------
        SEMs.append(SEM)
        stds.append(std)
        Means.append(Mean)
        #-------------------
    
    return Distances_fun, Means, SEMs, stds, xticks, Sizes_Plot
#======================================================================

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#1. PARAMETERS FOR INPUT/OUTPUT DATA
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#1.1. DATASETS PARAMETERS
#==================================================================
Datasets=['V-10_Tot','V-23_24_Tot']
K=10;L=20

Patients_per_Dataset=\
{'S-8_Tot':107, 'V-10_Tot':92, 'V-22_Tot':467, 'V-23_24_Tot':222,\
 'V-25_Tot':883}

Microbes_per_Dataset=\
{'S-8_Tot':128, 'V-10_Tot':137, 'V-22_Tot':134, 'V-23_24_Tot':118,\
 'V-25_Tot':144}
#================================================================== 

#1.2. PATHS FOR FILES
#=============================
path_In= '../../Input_Data/'
#=============================

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#2. READ DATA AND COMPUTE PARAMETERS
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Taxonomic_Orders=['Kingdom','Phylum','Class','Order','Family','Genus','Species']
Type_of_Metadata_Hosts=['Sex']


Intra_Averages=[]; Intra_SEMs=[]
Inter_Averages=[]; Inter_SEMs=[]

Intra_Averages_NM=[]; Intra_SEMs_NM=[]
Inter_Averages_NM=[]; Inter_SEMs_NM=[]

Total_Means_Intra_NM=[]; Total_SEMs_Intra_NM=[]
Total_Means_Inter_NM=[]; Total_SEMs_Inter_NM=[]

Statistics_Distance={}
Statistics_Distance_Hosts={}

for Data_read in Datasets:

    #2.1. Choose best likelihood
    #======================================================================
    path_Out='../../Output_Data/'+ 'Leave_One_Out_' + Data_read + '/'
    
    Likelihood_Seed={}
    for Seed in range(1,10,2):
        Likelihood_File=\
        'Dataset_' + Data_read + '_Seed_' + str(Seed) + '_LogLikelihood.txt'
        LikelihoodFinal=read_file(path_Out + Likelihood_File)[-1]
        Likelihood_Seed[float(LikelihoodFinal[0])]=Seed
        Seed=Likelihood_Seed[max(Likelihood_Seed)]
    #======================================================================

    #2.2. Parameters of Data
    #======================================================================
    Patients=Patients_per_Dataset[Data_read]
    Microbes=Microbes_per_Dataset[Data_read]
    
    File_Data='Dataset_' + Data_read + '_Seed_' + str(Seed) + '_Parameters.txt'
    Data=read_file(path_Out+File_Data)
    #======================================================================

    #2.3. Read and save membership vectors
    #=========================================
    theta=[];eta=[]
    theta=Data[:Patients]
    eta=Data[Patients+1:Patients+1+Microbes]
    #=========================================

    #2.4. Format microbe membership vectors to float
    #================================================================
    for microbe in range(len(eta)):
        eta[microbe] = [ float(membership) for membership in eta[microbe] ]

    for host in range(len(theta)):
        theta[host] = [ float(membership) for membership in theta[host] ]
    #================================================================


    #2.5. Sex metadata Data
    #=====================================================================

    #2.5.1. Read main file
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    Metadata_Hosts='Metadata_Hosts_' + Data_read + '.csv'
    Metadata_Hosts_Raw=pd.read_csv(path_In + Metadata_Hosts)
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    Sex_Data=Metadata_Hosts_Raw['Sex']
    
    #1. GROUP HOSTS BY SEX
    #-----------------------------------------------------------------------------
    Sex_Hosts_Dict={}
    for host_1 in range( len(Sex_Data) ):
        Id_theta_tuple=( host_1, theta[host_1] )
        
        try:
            Sex_Hosts_Dict[ Sex_Data[host_1] ].append( Id_theta_tuple )
        except KeyError:
            Sex_Hosts_Dict[ Sex_Data[host_1] ]=[ Id_theta_tuple ]

    #-----------------------------------------------------------------------------

    #2.5.2.3. COMPUTE EUCLIDEAN INTRA-DISTANCES FOR SAME SEX
    #-------------------------------------------------------------------------
    Intra_Host_Dict={}          
    Intra_Host_Dict=Measure_Intra_Distances(Sex_Hosts_Dict)
    Distances,Mean_Distances,SEM_Distances,std_Distances,phylums_xticks,sizes=\
    Extract_Distances(Intra_Host_Dict,1)

    Intra_Flat_List = [item for sublist in Distances for item in sublist]
    Intra_Flat_List = np.asarray(Intra_Flat_List)

    #Compute Average and SEM
    #..............................................................
    Intra_Average=sum(Intra_Flat_List)/float(len(Intra_Flat_List))
    Intra_SEM=scipy.stats.sem(Intra_Flat_List,ddof=0)
    #..............................................................

    #-----------------------------------------------------------------------------

    #3.3.2.3. COMPUTE EUCLIDEAN INTER-DISTANCES FOR PAIRS OF DIFFERENT SEX
    #-------------------------------------------------------------------------
    Inter_Host_Dict={}
    Inter_Host_Dict=Measure_Inter_Distances(Sex_Hosts_Dict)

    Distances_Inter,Mean_Distances_Inter,SEM_Distances_Inter,\
    std_Distances_Inter,phylums_xticks_Inter,sizes_Inter=Extract_Distances(Inter_Host_Dict,1)

    Inter_Flat_List = [item for sublist in Distances_Inter for item in sublist]
    Inter_Flat_List = np.asarray(Inter_Flat_List)

    #Compute Average and SEM
    #..............................................................
    Inter_Average=sum(Inter_Flat_List)/float(len(Inter_Flat_List))
    Inter_SEM=scipy.stats.sem(Inter_Flat_List,ddof=0)
    #..............................................................

    #-----------------------------------------------------------------------------

    #2.5.2.4. Save intra/inter means and SEMs
    #-------------------------------------------------------------------------
    Statistics_Distance_Hosts[Data_read]={\
    'Same Sex':Intra_Flat_List,'Opposite Sex':Inter_Flat_List}
    #-------------------------------------------------------------------------


#Convert dictionary to dataframe for plots
#--------------------------------------------------------------------------
List_df=[]
for key1 in Statistics_Distance_Hosts:
    print(key1)
    for key2 in Statistics_Distance_Hosts[key1]:
        print(key2)
        if key1=='V-10_Tot':
            Name='Qin'
        else:
            Name='H/Ll-P'
        for element in Statistics_Distance_Hosts[key1][key2]:
        	List_df.append([Name,key2,element])


Dataframe_Test=pd.DataFrame(List_df,columns=["Dataset","Type of Distance","Euclidian Distance"])
#--------------------------------------------------------------------------
print(Dataframe_Test)
#Plot
#--------------------------------------------------------------------------
sns.set(font_scale=1.2)
sns.set_style("ticks")

sns.barplot(x="Dataset", y="Euclidian Distance", hue="Type of Distance", data=Dataframe_Test)
plt.legend(loc='upper left',frameon=False)
sns.despine(top=True, right=True, left=False, bottom=False)
plt.show()
#--------------------------------------------------------------------------
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#3. SAVE DATA
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
pickle.dump( Statistics_Distance_Hosts, open( "Pairs_of_Distances_Hosts.p", "wb"))

elapsed_time = time.time() - start_time
print elapsed_time

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
