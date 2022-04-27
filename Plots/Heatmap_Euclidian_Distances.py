#13/05/2019. Hacemos un heatmap que calcula las distancias euclideas entre microbios de diferentes ordenes 
#taxonomicos. Es un nuevo intento de ver si existen diferencias entre las distancias intra e interclade.

import sys
import os
import numpy as np
import random as rnd
import pandas as pd
import copy as cp
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
import seaborn as sns
import scipy.cluster
from scipy import stats
from scipy.spatial import distance
import math
import operator
from collections import Counter
import time

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

#0.2. Measure pairwise distance
#=======================================================================
def Pairwise_Distance(vector1,vector2):
    vector1=np.array(vector1)
    vector2=np.array(vector2)
    Distance=np.linalg.norm(vector1 - vector2)

    return Distance
#=======================================================================

#0.2. Measure Euclidean distances
#============================================================================================
def Measure_Distances(Dict1):
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
#============================================================================================

#0.3. Measure Inter-Order distances
#============================================================================================
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

                   Ids_distance_tuple=(\
                   microbe_0[0], microbe_1[0], Dist)

                   try:
                       Dict_Inter[order_0+'-'+order_1].append(Ids_distance_tuple)

                   except KeyError:
                       Dict_Inter[order_0+'-'+order_1]=[Ids_distance_tuple]
    
    return Dict_Inter
#============================================================================================

#0.3. Extract distances from dictionary
#======================================================================
def Extract_Distances(Dict_Distances):

    Distances_fun=[]
    Means=[];SEMs=[];stds=[]
    xticks=[]; Sizes_Plot=[]

    for phylum in Dict_Distances:
        
        #Extract names of phyla and sizes of groups for plotting purposes
        #----------------------------------------------
        phylum_ticks=phylum.split('__')
        xticks.append(phylum_ticks[1])
        Sizes_Plot.append( len(Dict_Distances[phylum]))
        #----------------------------------------------

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

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#1. PARAMETERS FOR INPUT/OUTPUT DATA
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#1.1. DATASETS PARAMETERS
#==================================================================

#1.1. EXTERNAL PARAMETERS
#=====================
Dataset=sys.argv[1]#==
Diseases=0         #==
#=====================

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

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#2. PHYLOGENETIC DATA
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Taxonomic_Orders=['Kingdom','Phylum','Class','Order','Family','Genus','Species']

Taxonomic_level=6
Taxonomic_Order=Taxonomic_Orders[Taxonomic_level]
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#3. READ DATA AND STORE DATA
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#3.1. Choose best likelihood
#======================================================================
path_Out='../../Output_Data/'+ 'Leave_One_Out_' + Dataset + '/'

Likelihood_Seed={}
for Seed in range(1,10,2):
    Likelihood_File=\
    'Dataset_' + Dataset + '_Seed_' + str(Seed) + '_LogLikelihood.txt'
    LikelihoodFinal=read_file(path_Out + Likelihood_File)[-1]
    Likelihood_Seed[float(LikelihoodFinal[0])]=Seed
    Seed=Likelihood_Seed[max(Likelihood_Seed)]
#======================================================================

#3.2. Parameters of Data
#=======================================
Patients=Patients_per_Dataset[Dataset]
Microbes=Microbes_per_Dataset[Dataset]
#=======================================


#3.3. Taxonomic Data
#=====================================================================================================
File_Taxon='Names_Microbes_' + Dataset + '.csv'
Taxon_Data=pd.read_csv(path_In+File_Taxon,header=None)

Species_Data=Taxon_Data[Taxon_Data.columns[Taxonomic_level]] #We choose the phylum classification
Species_List=Species_Data.tolist()

Taxonomy_Hierarchy={}
for index, row in Taxon_Data.iterrows():
    Kingdom=row[0];Phylum=row[1];Class=row[2];Order=row[3];
    Family=row[4];Genus=row[5];Species=row[6];

    #Create taxonomic superdictionary
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    try:
        Taxonomy_Hierarchy[Phylum][Class][Order][Family][Genus].append(Species)
    except KeyError:
        try:
            Taxonomy_Hierarchy[Phylum][Class][Order][Family][Genus]=[Species]
        except KeyError:
            try:
                Taxonomy_Hierarchy[Phylum][Class][Order][Family]={ Genus:[Species] }
            except KeyError:
                try:
                    Taxonomy_Hierarchy[Phylum][Class][Order]={ Family: {Genus:[Species]} }
                except KeyError:
                    try:
                        Taxonomy_Hierarchy[Phylum][Class]={ Order:{ Family: {Genus:[Species]} }}
                    except KeyError:
                        Taxonomy_Hierarchy[Phylum]={Class: {Order:{ Family: {Genus:[Species]} }}}
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


Ordered_Species_RAW=[]
for Phylum in Taxonomy_Hierarchy:
    for Class in Taxonomy_Hierarchy[Phylum]:
        for Order in Taxonomy_Hierarchy[Phylum][Class]:
            for Family in Taxonomy_Hierarchy[Phylum][Class][Order]:
                
                # if "noname" in Family:
                #     continue

                for Genus in Taxonomy_Hierarchy[Phylum][Class][Order][Family]:
                    
                    # if "noname" in Genus:
                    #     continue
                    
                    Ordered_Species_RAW.append( Taxonomy_Hierarchy[Phylum][Class][Order][Family][Genus] )

print Ordered_Species_RAW, '\n'
Ordered_Species = [item for sublist in Ordered_Species_RAW  for item in sublist]
print Ordered_Species
print len(Ordered_Species)
#=====================================================================================================

#3.4. Model data
#==========================================================================    
File_Data='Dataset_' + Dataset + '_Seed_' + str(Seed) + '_Parameters.txt'
Data=read_file(path_Out+File_Data)
print Dataset
print Patients
print Microbes

#3.4.1. Read and save membership vectors
#::::::::::::::::::::::::::::::::::::::::::
theta=[];eta=[]
theta=Data[:Patients]
eta=Data[Patients+1:Patients+1+Microbes]
#::::::::::::::::::::::::::::::::::::::::::

#3.4.2. Format microbe membership vectors to float
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
for microbe in range(len(eta)):
    eta[microbe] = [float(membership) for membership in eta[microbe] ]
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#==========================================================================

#3.5 Dictionary of membership vectors
#=====================================================================
Microbial_Species_Dict={}
for microbe_1 in range( len(Species_List) ):
    Microbial_Species_Dict[ Species_List[microbe_1] ]= eta[microbe_1]
#=====================================================================
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#4. MAKE HEATMAP
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Distances=np.zeros(( len(Ordered_Species),len(Ordered_Species) ))

counter_i=0;counter_j=0
for microbe_i in Ordered_Species:
    for microbe_j in Ordered_Species:

        vector_i=Microbial_Species_Dict[microbe_i]
        vector_j=Microbial_Species_Dict[microbe_j]
        Distance_ij=Pairwise_Distance(vector_i, vector_j)

        Distances[counter_i][counter_j]=Distance_ij

        counter_j+=1

    counter_j=0
    counter_i+=1

print Distances
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#7. PLOTS
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#Heatmap
#=======================================
sns.heatmap(Distances,vmin=0,vmax=math.sqrt(2))
#=======================================

#title
#=======================================
plt.title("Euclidean Distances "+Dataset )
#=======================================

#Saving
#=======================================
path_out='/export/home/shared/Projects/Microbiome/Plots/First_Path/'
plt.savefig(path_out + 'Heatmap_Distances_'+ Dataset+\
             '.pdf',dpi=300)
#=======================================

plt.show()
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
