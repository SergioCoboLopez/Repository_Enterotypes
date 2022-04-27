#25/03/2019. En este codigo sacamos, para cada dataset, la media de sus distancias intra-clade y de sus distancias
#inter-clade con sus correspondientes SEMs. Comparamos con un null model para comprobar la robustness del
#resultado.

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

#0.4. Extract distances from dictionary
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
Datasets=['S-8_Tot','V-10_Tot','V-22_Tot','V-23_24_Tot','V-25_Tot' ]
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

#2. READ DATA AND COMPUTE PARAMETERS
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Taxonomic_Orders=['Kingdom','Phylum','Class','Order','Family','Genus','Species']


Intra_Averages=[]; Intra_SEMs=[]
Inter_Averages=[]; Inter_SEMs=[]

Intra_Averages_NM=[]; Intra_SEMs_NM=[]
Inter_Averages_NM=[]; Inter_SEMs_NM=[]

Total_Means_Intra_NM=[]; Total_SEMs_Intra_NM=[]
Total_Means_Inter_NM=[]; Total_SEMs_Inter_NM=[]

Statistics_Distance={}

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
    print Data_read
    print Patients
    print Microbes
    #======================================================================

    #2.3. Read and save membership vectors
    #=========================================
    theta=[];eta=[]
    theta=Data[:Patients]
    eta=Data[Patients+1:Patients+1+Microbes]
    #=========================================

    #2.4. Format microbe membership vectors to float
    #=====================================================================
    for microbe in range(len(eta)):
        eta[microbe] = [ float(membership) for membership in eta[microbe] ]
    #=====================================================================


    #2.5. Taxonomic Data
    #=====================================================================

    #2.5.1. Read main file
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    File_Taxon='Names_Microbes_' + Data_read + '.csv'
    Taxon_Data_Raw=pd.read_csv(path_In + File_Taxon,header=None)
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #2.5.2. Run Over Taxonomic Layers
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    for Clade in range(1,len(Taxonomic_Orders)-1):

        #2.5.2.1. Read data and format to list
        #--------------------------------------------------------
        Taxon_Data=Taxon_Data_Raw[Taxon_Data_Raw.columns[Clade]] #We choose the phylum classification
        Taxon_List=Taxon_Data.tolist()
        #--------------------------------------------------------

        #2.5.2.2. GROUP MICROBES BY TAXONOMIC ORDER
        #--------------------------------------------------------------------------------
        Taxonomy_Microbes_Dict={}
        for microbe_1 in range( len(Taxon_List) ):
            Id_eta_tuple=( microbe_1, eta[microbe_1] )

            try:
                Taxonomy_Microbes_Dict[ Taxon_List[microbe_1] ].append( Id_eta_tuple )
            except KeyError:
                Taxonomy_Microbes_Dict[ Taxon_List[microbe_1] ]=[ Id_eta_tuple ]
        #--------------------------------------------------------------------------------

        #2.5.2.3. COMPUTE EUCLIDEAN INTRA-DISTANCES FOR PAIRS OF ALL CLADES
        #--------------------------------------------------------------------------------
        Distances_Intra_Dict={}
        Distances_Intra_Dict=Measure_Distances(Taxonomy_Microbes_Dict)

        Distances, Mean_Distances, SEM_Distances, std_Distances, phylums_xticks, sizes=\
        Extract_Distances(Distances_Intra_Dict)

        Intra_Flat_List = [item for sublist in Distances for item in sublist]
        Intra_Flat_List = np.asarray(Intra_Flat_List)

        #Compute Average and SEM
        #..............................................................
        Intra_Average=sum(Intra_Flat_List)/float(len(Intra_Flat_List))
        Intra_SEM=scipy.stats.sem(Intra_Flat_List,ddof=0)
        #..............................................................

        #--------------------------------------------------------------------------------

        #3.3.2.3. COMPUTE EUCLIDEAN INTER-DISTANCES FOR PAIRS OF DIFFERENT CLADES
        #--------------------------------------------------------------------------------
        Distances_Inter_Dict={}
        Distances_Inter_Dict=Measure_Inter_Distances(Taxonomy_Microbes_Dict)

        Distances_Inter,Mean_Distances_Inter,SEM_Distances_Inter,\
        std_Distances_Inter,phylums_xticks_Inter,sizes_Inter=Extract_Distances(Distances_Inter_Dict)

        Inter_Flat_List = [item for sublist in Distances_Inter for item in sublist]
        Inter_Flat_List = np.asarray(Inter_Flat_List)

        #Compute Average and SEM
        #..............................................................
        Inter_Average=sum(Inter_Flat_List)/float(len(Inter_Flat_List))
        Inter_SEM=scipy.stats.sem(Inter_Flat_List,ddof=0)
        #..............................................................

        #--------------------------------------------------------------------------------

        #2.5.2.4. Save intra/inter means and SEMs
        #---------------------------------------------------------------------------------------------------------
        try:
            Statistics_Distance[Data_read][Taxonomic_Orders[Clade]]={\
            'Intra':Intra_Flat_List,'Inter':Inter_Flat_List}

        except KeyError:
            Statistics_Distance[Data_read]={ Taxonomic_Orders[Clade]:\
            {'Intra':Intra_Flat_List,'Inter':Inter_Flat_List} }

        #---------------------------------------------------------------------------------------------------------
    #=============================================================================================================

    #2.6. NULL MODEL
    #===========================================================================================================
    iterations=10 #Number of randomizations
    Means_of_randomizations_Intra=[]; SEMs_of_randomizations_Intra=[]
    Means_of_randomizations_Inter=[]; SEMs_of_randomizations_Inter=[]
    
    for iteration in range(iterations):
        
        #2.6.1. Shuffle membership microbes vector
        #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        eta_null_model=cp.deepcopy(eta)
        rnd.shuffle(eta_null_model)

        Taxonomy_Microbes_Null_Model_Dict={}
        for microbe_2 in range( len(Taxon_List) ):
            Id_eta_tuple_Null_Model=( microbe_2, eta_null_model[microbe_2] )

            try:
                Taxonomy_Microbes_Null_Model_Dict[ Taxon_List[microbe_2] ].append( Id_eta_tuple_Null_Model )
            except KeyError:
                Taxonomy_Microbes_Null_Model_Dict[ Taxon_List[microbe_2] ]=[ Id_eta_tuple_Null_Model ]    
        #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        #2.6.1. Shuffle membership microbes vector
        #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        Distances_Intra_Null_Model_Dict={}
        Distances_Intra_Null_Model_Dict=Measure_Distances(Taxonomy_Microbes_Null_Model_Dict)

        Distances_NM, Mean_Distances_NM, SEM_Distances_NM, std_Distances_NM, phylums_xticks_NM, sizes_NM=\
        Extract_Distances(Distances_Intra_Null_Model_Dict)

        Intra_Flat_List_Null_Model = [item for sublist in Distances_NM for item in sublist]
        Intra_Flat_List_Null_Model = np.asarray(Intra_Flat_List_Null_Model)

        #2.6.3. Compute statistics
        #--------------------------------------------------------------
        Intra_Average_NM=sum(Intra_Flat_List_Null_Model)/float(len(Intra_Flat_List_Null_Model))
        Means_of_randomizations_Intra.append(Intra_Average_NM)
        #--------------------------------------------------------------

        #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        #2.6.3. Inter-distances
        #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        Distances_Inter_Null_Model_Dict={}
        Distances_Inter_Null_Model_Dict=Measure_Inter_Distances(Taxonomy_Microbes_Null_Model_Dict)

        Distances_Inter_NM,Mean_Distances_Inter_NM,SEM_Distances_Inter_NM,std_Distances_Inter_NM,\
        phylums_xticks_Inter_NM,sizes_Inter_NM= Extract_Distances(Distances_Inter_Null_Model_Dict)

        Inter_Flat_List_Null_Model = [item for sublist in Distances_Inter_NM for item in sublist]
        Inter_Flat_List_Null_Model = np.asarray(Inter_Flat_List_Null_Model)


        #2.6.3.1. Compute statistics
        #-------------------------------------------------------------
        Inter_Average_NM=sum(Inter_Flat_List_Null_Model)/float(len(Inter_Flat_List_Null_Model))
        Means_of_randomizations_Inter.append(Inter_Average_NM)
        #-------------------------------------------------------------
        
        #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        
    #2.6.4. Compute Statistics for all randomizations
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    Mean_of_randomizations_Intra=sum(Means_of_randomizations_Intra)/float(len(Means_of_randomizations_Intra))
    Mean_of_randomizations_Inter=sum(Means_of_randomizations_Inter)/float(len(Means_of_randomizations_Inter))

    SEM_of_randomization_Intra=scipy.stats.sem(Means_of_randomizations_Intra, ddof=0)
    SEM_of_randomization_Inter=scipy.stats.sem(Means_of_randomizations_Inter, ddof=0)
    
    Total_Means_Intra_NM.append(Mean_of_randomizations_Intra)
    Total_SEMs_Intra_NM.append(SEMs_of_randomizations_Intra)
        
    Total_Means_Inter_NM.append(Mean_of_randomizations_Inter)
    Total_SEMs_Inter_NM.append(SEMs_of_randomizations_Inter)
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #===========================================================================================================
    
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#3. SAVE DATA
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
pickle.dump( Statistics_Distance, open( "Pairs_of_Distances.p", "wb" ) )

elapsed_time = time.time() - start_time
print elapsed_time

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
