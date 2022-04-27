#25/03/2019. Hacemos diferentes plots para ver la diferencia de distancias euclideas intra-clade y inter-clade. Es
#decir, comparamos los parametros del modelo con la informacion taxonomica.

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

#1.1. EXTERNAL PARAMETERS
#===================================
Dataset=sys.argv[1]              #==
Type_of_Randomization=sys.argv[2]#==
#===================================

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
path_Out='../../Output_Data/'+ 'Leave_One_Out_' + Dataset + '/'
path_In= '../../Input_Data/'
#==============================================================

#1.4. Choose best likelihood
#======================================================================
Likelihood_Seed={}
for Seed in range(1,10,2):
    Likelihood_File=\
    'Dataset_' + Dataset + '_Seed_' + str(Seed) + '_LogLikelihood.txt'
    LikelihoodFinal=read_file(path_Out + Likelihood_File)[-1]
    Likelihood_Seed[float(LikelihoodFinal[0])]=Seed

print Likelihood_Seed
print max(Likelihood_Seed)
print Likelihood_Seed[max(Likelihood_Seed)]
Seed=Likelihood_Seed[max(Likelihood_Seed)]
#======================================================================

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#2. READ PHYLOGENETIC DATA
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Taxonomic_Orders=['Kingdom','Phylum','Class','Order','Family','Genus','Species']

Taxonomic_level=2
Taxonomic_Order=Taxonomic_Orders[Taxonomic_level]

File_Taxon='Names_Microbes_' + Dataset + '.csv'
Taxon_Data=pd.read_csv(path_In+File_Taxon,header=None)
Taxon_Data=Taxon_Data[Taxon_Data.columns[Taxonomic_level]] #We choose the phylum classification
Taxon_List=Taxon_Data.tolist()
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#3. READ MODEL PARAMETERS
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
File_Data='Dataset_' + Dataset + '_Seed_' + str(Seed) + '_Parameters.txt'
Data=read_file(path_Out+File_Data)

#3.1. Read and save membership vectors
#:::::::::::::::::::::::::::::::::::::::::
theta=[];eta=[]
theta=Data[:Patients]
eta=Data[Patients+1:Patients+1+Microbes]
#:::::::::::::::::::::::::::::::::::::::::

#3.2. Format microbe membership vectors to float
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
for microbe in range(len(eta)):
    eta[microbe] = [float(membership) for membership in eta[microbe] ]
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#4. GROUP MICROBES BY TAXONOMIC ORDER
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Taxonomy_Microbes_Dict={}
for microbe_1 in range( len(Taxon_List) ):
    Id_eta_tuple=( microbe_1, eta[microbe_1] )
    
    try:
        Taxonomy_Microbes_Dict[ Taxon_List[microbe_1] ].append( Id_eta_tuple )
    except KeyError:
        Taxonomy_Microbes_Dict[ Taxon_List[microbe_1] ]=[ Id_eta_tuple ]
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#5. Intra-distances
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#5.1. Compute distances
#========================================================
Distances_Dict={}
Distances_Dict=Measure_Distances(Taxonomy_Microbes_Dict)
#========================================================

#5.2 Extract distances as vectors
#===============================================================================
Distances, Mean_Distances, SEM_Distances, std_Distances, phylums_xticks, sizes=\
Extract_Distances(Distances_Dict)

Distances_Flat_List = [item for sublist in Distances for item in sublist]

for List in Distances:
    List=np.asarray(List)
#===============================================================================

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#6. COMPUTE EUCLIDEAN DISTANCES OF PAIRS INTER-GENUS
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Distances_Inter_Dict={}
Orders=Taxonomy_Microbes_Dict.keys()

#6.1. Run by pairs of different clades to compute distances
#==========================================================================================
for order_0 in Taxonomy_Microbes_Dict:
    Orders.remove(order_0)

    for order_1 in Orders:
        
        contador_parejas=0
        for microbe_0 in Taxonomy_Microbes_Dict[order_0]:
            for microbe_1 in Taxonomy_Microbes_Dict[order_1]:

               Dist_0=np.array(microbe_0[1])
               Dist_1=np.array(microbe_1[1])
               Dist = np.linalg.norm(Dist_0 - Dist_1)

               Ids_distance_tuple=(\
               microbe_0[0], microbe_1[0], Dist)

               try:
                   Distances_Inter_Dict[order_0+'-'+order_1].append(Ids_distance_tuple)

               except KeyError:
                   Distances_Inter_Dict[order_0+'-'+order_1]=[Ids_distance_tuple]

               
               contador_parejas+=1
#==========================================================================================

#6.2. Extract distances as vectors
#==============================================================================================================
Distances_Inter,Mean_Distances_Inter,SEM_Distances_Inter,std_Distances_Inter,phylums_xticks_Inter,sizes_Inter=\
Extract_Distances(Distances_Inter_Dict)
#==============================================================================================================

Inter_Flat_List = [item for sublist in Distances_Inter for item in sublist]
Inter_Flat_List = np.asarray(Inter_Flat_List)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#7. RANDOMIZE PAIRS OF MICROBES
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#7.1. Choose number of randomizations,  fix seed, and type of randomization
#==========================================================================
Randomizations=10
rnd.seed(1)
#Normal, Intra-Inter
#Type_of_Randomization="Normal"
#==========================================================================

#7.2. Define dictionaries to store the statistics
#======================================================
Distances_rand_Dict={};Mean_Distances_rand_Dict={};
SEM_Distances_rand_Dict={}; std_Distances_rand_Dict={}
#======================================================

Randomized_Dataframes= [[] for _ in range(Randomizations)]

#7.3. Intra-Inter Randomization
#=======================================================================================
if Type_of_Randomization=="Intra-Inter":

    #7.3.1. Declare dictionaries for Intra vs Inter Randomization
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    Intra_Inter_Distances_Dict={}
    Number_of_Pairs={ key:len(Distances_Dict[key]) for key in Distances_Dict }
    print Number_of_Pairs    
    Instances_Taxonomic_Order=Number_of_Pairs.keys()
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #7.3.2. Do many randomizations
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    for iteration in range(Randomizations):

        Intra_Inter_Distances_Dict={}

        #7.3.2.1. Generate pool of inter microbes and dictionary of intra-inter pairs
        #----------------------------------------------------------------------------
        for Instance in Instances_Taxonomic_Order:

            #Pool of microbes from inter orders
            #.....................................................
            Pool_of_Orders=cp.copy(Instances_Taxonomic_Order)
            Pool_of_Orders.remove(Instance)

            Pool_of_InterMicrobes=[]
            for Order in Pool_of_Orders:
                for microbe_2 in Taxonomy_Microbes_Dict[Order]:
                    Pool_of_InterMicrobes.append( microbe_2 )
            #.....................................................

            #Dictionary of Intra-Inter pairs with distances
            #.........................................................................
            for Pair in range(Number_of_Pairs[Instance]):

                #Take an intra and inter microbe at random
                #-----------------------------------------------------------
                Microbe_Intra=rnd.choice( Taxonomy_Microbes_Dict[Instance] )
                Microbe_Inter=rnd.choice( Pool_of_InterMicrobes )
                #-----------------------------------------------------------

                #Compute euclidean distance
                #----------------------------
                a=np.array(Microbe_Intra[1])
                b=np.array(Microbe_Inter[1])
                dist = np.linalg.norm(a - b)
                #----------------------------

                #Save to dictionary
                #-------------------------------------------------------------------
                Ids_distance_tuple=(Microbe_Intra[0], Microbe_Inter[0], dist)

                try:
                    Intra_Inter_Distances_Dict[Instance].append(Ids_distance_tuple)

                except KeyError:
                    Intra_Inter_Distances_Dict[Instance]=[Ids_distance_tuple]
                #-------------------------------------------------------------------

        #----------------------------------------------------------------------------

            
        #7.3.2.2. Convert Dictionary to Dataframes
        #--------------------------------------------------------------------------
        for phylum in Intra_Inter_Distances_Dict:

            phylum_full_name=phylum.split('__')
            phylum1=phylum_full_name[1]

            Distance_for_dataframe_rand=[]
            for distance in Intra_Inter_Distances_Dict[phylum]:
                Distance_for_dataframe_rand.append(distance[2])

            df_random = pd.DataFrame({'Euclidean Distance': Distance_for_dataframe_rand,
                        Taxonomic_Order: phylum1,
                        'Type of Data': "Randomized" } )
            Randomized_Dataframes[iteration].append(df_random)

        Randomized_Dataframes[iteration]=pd.concat(Randomized_Dataframes[iteration])
        #---------------------------------------------------------------------------

        #7.3.3.4. Extract distances from dictionaries
        #-------------------------------------------------------------------
        Distances_rand, Mean_Distances_rand, SEM_Distances_rand, std_Distances_rand, phylums_xticks_rand,Sizes=\
        Extract_Distances(Intra_Inter_Distances_Dict)

        for List_rand in Distances_rand:
            List_rand=np.asarray(List_rand)
        #-------------------------------------------------------------------

        #7.3.3.5. Allocate distances in Dictionaries
        #-------------------------------------------------------------------
        Distances_rand_Dict[iteration]=Distances_rand;Mean_Distances_rand_Dict[iteration]= Mean_Distances_rand;
        SEM_Distances_rand_Dict[iteration]=SEM_Distances_rand;
        std_Distances_rand_Dict[iteration]=std_Distances_rand
        #-------------------------------------------------------------------

        #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#7.4. Normal Randomization
#=======================================================================================
elif Type_of_Randomization=="Normal":

    #7.4.1. Pool of all microbes
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    Pool_of_Microbes = [item for sublist in Taxonomy_Microbes_Dict.values() for item in sublist]
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #7.4.2. Do many randomizations
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    for iteration in range(Randomizations):
        
        #7.4.2.1. Choose random samples of microbes
        #-------------------------------------------------------------------
        Taxonomy_Microbes_Dict_Random={}
        for phylum in Taxonomy_Microbes_Dict:
            Taxonomy_Microbes_Dict_Random[phylum]=rnd.sample(Pool_of_Microbes,len(Taxonomy_Microbes_Dict[phylum]))
        #-------------------------------------------------------------------

        #7.4.2.2. Measure euclidean distances of randomized pairs for each phylum
        #-------------------------------------------------------------------
        Distances_Dict_Randomized={}
        Distances_Dict_Randomized=Measure_Distances(Taxonomy_Microbes_Dict_Random)
        #-------------------------------------------------------------------

        #7.4.2.3. Save distances to dataframes for plotting
        #-------------------------------------------------------------------
        for phylum in Distances_Dict_Randomized:

            phylum_full_name=phylum.split('__')
            phylum1=phylum_full_name[1]

            Distance_for_dataframe_rand=[]
            for distance in Distances_Dict_Randomized[phylum]:
                Distance_for_dataframe_rand.append(distance[2])

            df_random = pd.DataFrame({'Euclidean Distance': Distance_for_dataframe_rand,
                        Taxonomic_Order: phylum1,
                        'Type of Data': "Randomized" } )
            Randomized_Dataframes[iteration].append(df_random)

        Randomized_Dataframes[iteration]=pd.concat(Randomized_Dataframes[iteration])
        #-------------------------------------------------------------------

        #7.4.2.3. Extract distances from dictionaries
        #-------------------------------------------------------------------
        Distances_rand, Mean_Distances_rand, SEM_Distances_rand, std_Distances_rand, phylums_xticks_rand,Sizes=\
        Extract_Distances(Distances_Dict_Randomized)

        for List_rand in Distances_rand:
            List_rand=np.asarray(List_rand)
        #-------------------------------------------------------------------

        #7.4.2.5. Allocate distances in Dictionaries
        #-------------------------------------------------------------------
        Distances_rand_Dict[iteration]=Distances_rand;Mean_Distances_rand_Dict[iteration]= Mean_Distances_rand;
        SEM_Distances_rand_Dict[iteration]=SEM_Distances_rand;
        std_Distances_rand_Dict[iteration]=std_Distances_rand
        #-------------------------------------------------------------------

#=======================================================================================
        
#7.5. Concatenate randomized frameworks
#====================================================
Randomized_Dataframe=pd.concat(Randomized_Dataframes)
#====================================================

#7.6. Do general statistics
#========================================================================

#7.6.1. Calculate the mean of the means
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
All_Means= [[] for _ in range(len(Mean_Distances))]
for iteration_1 in Mean_Distances_rand_Dict:
    for phylum in range(len(All_Means)):
        (All_Means[phylum]).append(Mean_Distances_rand_Dict[iteration_1][phylum])

Total_Mean_Rand=[];Total_std_Rand=[];Total_SEM_Rand=[]
for Total_List in All_Means:
    Total_Mean_Rand.append( sum(Total_List)/len(Total_List) )
    Total_std_Rand.append(np.std(Total_List))
    Total_SEM_Rand.append(stats.sem(Total_List,ddof=0))
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#========================================================================

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#8. Convert data to pandas dataframes
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#8.5.1. Convert OBSERVED Data to dataframe
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Observed_Dataframes=[]
for phylum in Distances_Dict:
    
    phylum_full_name=phylum.split('__')
    phylum1=phylum_full_name[1]
    
    Distance_for_dataframe=[]
    Phylum=[]
    for distance in Distances_Dict[phylum]:
        Distance_for_dataframe.append(distance[2])

    df1 = pd.DataFrame({'Euclidean Distance': Distance_for_dataframe,
                        Taxonomic_Order: phylum1,
                        'Type of Data': "Observed" } )
    Observed_Dataframes.append(df1)

Observed_Dataframe=pd.concat(Observed_Dataframes)
Total_Dataframe=pd.concat([Observed_Dataframe,Randomized_Dataframe])
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#8.5.2. Convert MEANS OF RANDOMIZED Data to dataframe
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Random_Mean_Distances=[];Random_phyla=[]
for iteration in Mean_Distances_rand_Dict:
    
    for random_phylum,distance in zip(phylums_xticks_rand,Mean_Distances_rand_Dict[iteration]):
        Random_Mean_Distances.append(distance)
        Random_phyla.append(random_phylum)

Randomized_Mean_Dataframe=pd.DataFrame({'Euclidean Distance': Random_Mean_Distances,
                  Taxonomic_Order: Random_phyla,
                  'Type of Data': "Randomized"})
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#9. COMPUTE P-VALUES
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#9.1. Values over the mean
#=================================================================
Counter_Significance=[0]*len(Mean_Distances)
for iteration in Mean_Distances_rand_Dict:

    for phylum in range(len(Mean_Distances_rand_Dict[iteration])):
        Random=Mean_Distances_rand_Dict[iteration][phylum]
        Observed=Mean_Distances[phylum]
        if Observed > Random:
            Counter_Significance[phylum]+=1
#=================================================================

#9.2. P-Values
#=======================================
PValues=[]
for value in Counter_Significance:
    PValues.append( value/float(1000))
#=======================================

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#10. PLOTS
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
xticks_pos   =[tick for tick in range(1,len(phylums_xticks) + 1 )]
xticks_pos_violin   =[tick for tick in range(0,len(phylums_xticks))]
xticks_labels=phylums_xticks
path_out='/export/home/shared/Projects/Microbiome/Plots/'
path_out_OVD='../../Plots/'

#10.1. Vector de tamanyos de puntos para el errorplot
#===================================================
Visible_Sizes=[]
for Size in Sizes:
    Normalized=Size/float(sum(Sizes))
    Visible_Sizes.append( Normalized*1000)
#===================================================
    
Information=["Double_Violin_Plot", "Total_Distribution","Mean+Randomized_Means"]
Plot=Information[1]

if Plot=="Double_Violin_Plot":
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    fig=plt.figure(figsize=(20,9))

    ax = sns.violinplot(x=Taxonomic_Order, y="Euclidean Distance",hue="Type of Data",
                    data=Total_Dataframe, palette="muted",split=True,cut=0)

    #Introduce Text in Plot
    #------------------------------------------------------------------------------------------------------------
    pos_x=0.3;pos_y=0.57
    for phylum_0 in Taxonomy_Microbes_Dict:
        if len(Taxonomy_Microbes_Dict[phylum_0])>1:
            ax.text(pos_x, pos_y, str(phylum_0[3:]) +': '+ str(len(Taxonomy_Microbes_Dict[phylum_0]))+' species')
            pos_y=pos_y - 0.07
    #------------------------------------------------------------------------------------------------------------
    
    
    plt.savefig(path_out + Dataset + '_' + Taxonomic_Order +'_'+ Type_of_Randomization +\
            '_Distance_Double_Vplot.pdf',dpi=300)
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if Plot=="Total_Distribution":
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    sns.distplot(Inter_Flat_List, color='r',kde=True)
    sns.distplot(Distances_Flat_List, kde=True)
    plt.xlabel('Euclidean Distance')
    plt.ylabel('Density')
    plt.legend(['Inter-Distances','Intra-Distances'])
    plt.title('Inter vs Intra distances for taxonomic level: ' + str(Taxonomic_Order))
    
    plt.savefig(path_out_OVD + Dataset + '_'+ Taxonomic_Order+\
                '_Intra-Inter_Total_Distribution.pdf',dpi=300)
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if Plot=="Mean+Randomized_Means":
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    fig=plt.figure(figsize=(20,9))

    #VIOLIN PLOT
    #-----------------------------------------------------------------------
    ax1=sns.violinplot(x=Taxonomic_Order, y="Euclidean Distance",#hue="Type of Data",
                    data=Randomized_Mean_Dataframe, palette="muted")
    
    #-----------------------------------------------------------------------

    #Hlines
    #-----------------------------------------------------------------------
    x0=[x-0.3 for x in xticks_pos_violin]
    x1=[x+0.3 for x in xticks_pos_violin]
    ax2=plt.hlines(Mean_Distances,x0,x1,colors='r')
    #-----------------------------------------------------------------------

    #Introduce Text in Plot
    #-------------------------------------------------------------------------------------
    plt.ylabel('Average Euclidian Distance between pairs of microbes',fontsize=15)
    plt.xlabel(Taxonomic_Order,fontsize=15)

    pos_x=-0.2;pos_y=min(Mean_Distances)-0.1
    for element in PValues:
    	ax1.text(pos_x, pos_y, 'P-Value: '+ str(element) , style='italic',
                 bbox={'facecolor':'red', 'alpha':0.5, 'pad':10})
        pos_x=pos_x + 1
    #-------------------------------------------------------------------------------------
    plt.ylim((0,1.4))
    plt.savefig(path_out + Dataset + '_'+ Taxonomic_Order+'_'+Type_of_Randomization+\
                '_Distance_Averages.pdf',dpi=300)
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

plt.show()
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
