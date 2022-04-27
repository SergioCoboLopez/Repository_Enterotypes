#23/01/2020. En este codigo miramos si los microbios generalistas lo son en todos los datasets. Para ello, miramos la correlacion entre \
#las medias de los rangos (filas que ocupan en la matriz nested) y \
#las SEMs que tienen.

import numpy as np
import random as rnd
import copy as cp
import sys
import matplotlib.gridspec as gridspec
import matplotlib
from pylab import *
from matplotlib import colors
import matplotlib.patches as patches
import pandas as pd
import scipy.cluster
from scipy import stats
from sklearn import datasets, linear_model
from sklearn.metrics import mean_squared_error, r2_score
from scipy.special import expit
from collections import OrderedDict


import seaborn as sns
import matplotlib.pyplot as plt
import operator
import time
import collections
#import jellyfish

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

#0.2. Compute Nestedness
#===================================================================
def Nestedness_fun(Matrix):

    #0.2.1. Reorder matrix by generalist-to-specialist criteria
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #Dictionary of rows: sum of columns
    #-------------------------------------------
    sum_cols={}
    for row in range(rows):
        sum_cols[row]=sum(Matrix[row])
    #-------------------------------------------

    #Order the dictionary by ascending sum of columns
    #-------------------------------------------
    sorted_sum_cols = \
    sorted(sum_cols.items(), key=operator.itemgetter(1),reverse=True)
    #-------------------------------------------

    #Reorder matrix rows by the sum of the columns
    #-------------------------------------------
    Dict_row_membership_fun={}
    Matrix_row_ordered=np.zeros(( rows,columns ))
    for new_row in range(rows):
        Matrix_row_ordered[new_row]=\
        Matrix[sorted_sum_cols[new_row][0] ]
        Dict_row_membership_fun[new_row]=sorted_sum_cols[new_row][0]
    #-------------------------------------------

    #Dictionary of cols: sum of rows
    #-------------------------------------------
    sum_rows={}
    for col in range(columns):
        sum_rows[col]=sum( np.transpose(Matrix_row_ordered)[col] )
    #-------------------------------------------

    #Order the dictionary by ascending sum of rows
    #------------------------------------------------
    sorted_sum_rows = \
    sorted(sum_rows.items(), key=operator.itemgetter(1),reverse=True)
    #------------------------------------------------ 

    #Reorder matrix rows by the sum of the columns
    #------------------------------------------------
    Dict_column_membership_fun={}
    Matrix_ordered=np.zeros(( rows,columns ))
    for new_col in range(columns):
        np.transpose(Matrix_ordered)[new_col] = \
        np.transpose(Matrix_row_ordered)[sorted_sum_rows[new_col][0]]
        Dict_column_membership_fun[new_col]=\
                sorted_sum_rows[new_col][0]
    #------------------------------------------------
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #0.2.2. Compute Nestedness from ordered matrix
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    d_ij_microbes=0;d_ij_patients=0;
    Min_Interactions_microbes=0;Min_Interactions_patients=0;

    Index=max(rows, columns)

    for row_i in range(Index-1):
        for row_j in range(row_i+1,Index):

            #Numerator
            #--------------------------------------------------------

            #Patients
            #........................................................
            if row_j<rows:
                product_patients=np.multiply(\
                Matrix_ordered[row_i], Matrix_ordered[row_j])
                shared_interactions_patients=np.sum(product_patients)
                d_ij_patients+=shared_interactions_patients
            #....................................................... 

            #Microbes
            #........................................................
            if row_j<columns:
                product_microbes=\
                np.multiply(np.transpose(Matrix_ordered)[row_i],\
                np.transpose(Matrix_ordered)[row_j])
                shared_interactions_microbes=np.sum(product_microbes)
                d_ij_microbes+=shared_interactions_microbes
            #........................................................

            #--------------------------------------------------------

            #Denominator
            #--------------------------------------------------------

            #Patients
            #........................................................
            if row_j<rows:
                Interactions_patients_i=np.sum(Matrix_ordered[row_i])
                Interactions_patients_j=np.sum(Matrix_ordered[row_j])
                Min_Interactions_patients+=\
                min(Interactions_patients_i, Interactions_patients_j)
            #........................................................

            #Microbes
            #........................................................
            if row_j<columns:
                Interactions_microbes_i=\
                np.sum(np.transpose(Matrix_ordered)[row_i])
                Interactions_microbes_j=\
                np.sum(np.transpose(Matrix_ordered)[row_j])
                Min_Interactions_microbes+=\
                min(Interactions_microbes_i, Interactions_microbes_j)
            #.......................................................

            #-------------------------------------------------------

    Nestedness=float(d_ij_microbes+d_ij_patients)/\
    float(Min_Interactions_microbes+Min_Interactions_patients)
    
    return Matrix_ordered, Nestedness, Dict_row_membership_fun, Dict_column_membership_fun
#===================================================================

#0.3. Randomized Matrix (Hard Null Model)
#====================================================================
def Null_Model_fun(Probability_Matrix):
    rows_null=np.shape(Probability_Matrix)[0]
    cols_null=np.shape(Probability_Matrix)[1]

    Null_Model=np.zeros((rows_null,cols_null))

    for row_null in range(rows_null):
        for col_null in range(cols_null):
            Candidate=rnd.random()

            if Candidate <= Probability_Matrix[row_null][col_null]:
                Null_Model[row_null][col_null]=1


    Null_Model_ordered,Nestedness_Null_Model,Dict_row_membership_fun, Dict_column_membership_fun=\
    Nestedness_fun(Null_Model)

    return Null_Model_ordered,Nestedness_Null_Model
#====================================================================

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#1. PARAMETERS FOR INPUT/OUTPUT DATA
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#1.1. EXTERNAL PARAMETERS
#=====================
number_of_links=3  #==
rnd.seed(1)        #==
#=====================

Microbe_Presence=4

#1.2. INTERNAL PARAMETERS
#===================================================================

#::::::::::::::::::::::::::::::::::::
K=10;L=20
Nested_Matrices={};All_Nestedness={}
Triangles={};Squares={}
Null_Models={}
#::::::::::::::::::::::::::::::::::::

Datasets=['S-8_Tot','V-10_Tot','V-22_Tot','V-23_24_Tot','V-25_Tot']

Patients_per_Dataset=\
{'S-8_Tot':107,'V-10_Tot':92,'V-22_Tot':467,'V-23_24_Tot':222,'V-25_Tot':883}

Microbes_per_Dataset=\
{'S-8_Tot':128,'V-10_Tot':137,'V-22_Tot':134,'V-23_24_Tot':118,'V-25_Tot':144}
#===================================================================

Microbes_sorted={};Microbes_size={}

for Dataset in Datasets:
    
    #1.1. Adjacency Matrix and taxonomies
    #===============================================================
    Patients=Patients_per_Dataset[Dataset]
    Microbes=Microbes_per_Dataset[Dataset]
    
    #1.1.1. Taxonomic Information
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    File_Taxon='Names_Microbes_' + Dataset + '.csv'
    Taxon_Data_Raw=pd.read_csv('../../Input_Data/'\
                        + File_Taxon,header=None)
    Species_Names=Taxon_Data_Raw[6]
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #1.1.2. Genetic Information
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    File_Genomics='Genomes_Full.csv'
    Genomic_Data_Raw=pd.read_csv('../../Input_Data/'+File_Genomics)
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #1.1.3. Missing microbes file
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    File_Missing_Microbes="Missing_Microbes.csv"
    Missing_Microbes_Record=\
    pd.read_csv('../../Input_Data/'+File_Missing_Microbes)
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
    
    AdjacencyMatrix_str=read_file(\
    '../../Input_Data/Binarized_Datasets/Binarized_'+Dataset+'.txt') 
    AdjacencyMatrix=[]
    for line in AdjacencyMatrix_str:
        new_line=[float(block) for block in line]
        AdjacencyMatrix.append(new_line)
        
    AdjacencyMatrix=np.asarray(AdjacencyMatrix)

    #Trasponemos por consistencia con los latent enterotypes
    #---------------------------------------------
    AdjacencyMatrix=np.transpose(AdjacencyMatrix)
    #---------------------------------------------

    global rows;   rows=AdjacencyMatrix.shape[0]
    global columns;columns=AdjacencyMatrix.shape[1]
    #===============================================================

    #1.2. Get Nestedness and ranks of microbes
    #================================================================
    NestedMatrix,Nestedness,Dict_row_membership,\
    Dict_column_membership=Nestedness_fun(AdjacencyMatrix)

    Inverted_Dict_column_membership = \
    dict((v,k) for k,v in Dict_column_membership.iteritems())
    #================================================================

    #1.1.2. Taxonomic Information
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    counter_microbes=0
    Missing_microbes_record_dataset=[]
    for microbe in Species_Names:
        
        #Cogemos informacion del dataframe
        #------------------------------------------------------------
        Matching_Microbes=Genomic_Data_Raw.loc\
        [Genomic_Data_Raw['#Organism Name']==microbe]

        Sizes_Raw=Matching_Microbes['Size(Mb)_y']

        #Si no encuentra la especie miramos a ver si se ha equivocado
        #............................................................
        if Matching_Microbes.empty==True:
            
            #We start looking at Vinod's File
            #------------------------------------------------------
            Alter_Info=Missing_Microbes_Record\
            [Missing_Microbes_Record['Total']==microbe]

            #If we find something on Vinod's file we enter there
            #-------------------------------------------------------
            if Alter_Info.empty==False:
                Alter_Info=Alter_Info.reset_index(drop=True)
                Alter_Name=Alter_Info['Alternative species names'][0]
                In_Dataset=Alter_Info['Found in dataset'][0]

                #If the microbe exists in our dataset with other
                #name
                if In_Dataset==1:

                    #Cambiar espacios por guion bajo y poner s__
                    #.............................................
                    Alter_Name=Alter_Name.replace("\xc2\xa0", " ")
                    Alter_Name=Alter_Name.replace(" ", "_")
                    Alter_Name="s__"+Alter_Name
                    #.............................................
                    
                    Matching_Microbes=Genomic_Data_Raw.loc\
                    [Genomic_Data_Raw['#Organism Name']==Alter_Name]
                    
                    Sizes_Raw=Matching_Microbes['Size(Mb)_y']

                #If the microbe does not exist in our dataset
                else:

                    Sizes_Raw=Alter_Info['Size']
            #-------------------------------------------------------

            #If we don't find anything on Vinod's file we look
            #somewhere else
            elif Alter_Info.empty==True:
                        
                #Pasar por todos los nombres del dataset de genomica
                #....................................................
                Names_Genomics=Genomic_Data_Raw['#Organism Name']
                for Name_Genomic in Names_Genomics:
                    microbe_splits=microbe.split("_")
                    Name_Genomic_splits=Name_Genomic.split("_")

                #Descompongo nombre en genus, species y compruebo que
                #son iguales. Ademas, compruebo que la specie no sea
                #"bacterium".Me quedo con el nombre resultante porque
                #no me consta que haya otro
                #....................................................
                    if microbe_splits[2]==Name_Genomic_splits[2]:
                        if ("_"+microbe_splits[3]+"_" in Name_Genomic and microbe_splits[3]!="bacterium"):
                            Matching_Microbes=Genomic_Data_Raw.loc\
                            [Genomic_Data_Raw['#Organism Name']==\
                             Name_Genomic]

                            Sizes_Raw=Matching_Microbes['Size(Mb)_y']
                            
        #------------------------------------------------------------
        
        #Compute mean of all microbes genetic info found
        #-----------------------------------------------
        Sizes=[Size for Size in Sizes_Raw]
        MeanSize=np.mean(Sizes)
        Microbes_size[microbe]=MeanSize
        #-----------------------------------------------
            
        try:
            Microbes_sorted[microbe][Dataset]=\
            Inverted_Dict_column_membership[counter_microbes]

        except KeyError:
            Microbes_sorted[microbe]={Dataset:\
            Inverted_Dict_column_membership[counter_microbes]}

        counter_microbes+=1
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
    print "Nestedness",Nestedness, "\n"

    Nested_Matrices[Dataset]=NestedMatrix
    All_Nestedness[Dataset]=Nestedness
    #================================================================
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#Compute average ranks and errors and save microbes by rank
#====================================================================
Microbes_sorted_means={};Microbes_ranks={};Rankings=[]
for microbe in Microbes_sorted:

    if len(Microbes_sorted[microbe])>4:

        Rankings.append( (np.mean(Microbes_sorted[microbe].values()),
        scipy.stats.sem(Microbes_sorted[microbe].values() ,ddof=0),\
                        np.std(Microbes_sorted[microbe].values()),
                        Microbes_size[microbe] ))
#====================================================================

#Save most generalist and most specialist to file
#====================================================================
sorted_ranks=\
sorted(Microbes_ranks.items(),key=operator.itemgetter(1))

Top_Specialists=sorted_ranks[:10];Top_Generalists=sorted_ranks[-10:]
with open('Specialists_and_Generalists.txt', 'w') as f:
    print >> f, "SPECIALISTS: \n"
    for Specialist in Top_Specialists:
        print >> f, Specialist[0][3:]

    print >> f, "\nGENERALISTS: \n"
    for Generalist in Top_Generalists:
        print >> f, Generalist[0][3:]

Names_Specialists=[]
for Specialist in Top_Specialists:
    Names_Specialists.append(Specialist[0][3:])

Names_Generalists=[]
for Generalist in Top_Generalists:
    Names_Generalists.append(Generalist[0][3:])

sorted_by_first = sorted(Rankings, key=lambda tup: tup[0])
#====================================================================

#Hacer Dataframe
#================================
Means=[];SEMs=[];stds=[];sizes=[]
for pair in sorted_by_first:
    Means.append(pair[0])
    SEMs.append(pair[1])
    stds.append(pair[2])
    sizes.append(pair[3])

d={"Mean Rank": Means, "SEM": SEMs, "Size(Mb)": sizes}
df=pd.DataFrame(d)
df.to_csv('Test.csv')
#================================

#Plot
#====================================================================

#==============
size_ticks=14
size_title=16                                                        
size_eje=16
size_letter=20
#==============

Letters=['a', 'b']
fig=plt.figure(figsize=(19, 7))

gs=gridspec.GridSpec(1, 2)
gs.update(left=0.06,right=0.96,bottom=0.14,top=0.88,wspace=0.2,\
hspace=0.45)

#Ranks vs SEMS
ax1=plt.subplot(gs[0,0])

ax1=sns.regplot(x=df["Mean Rank"],y=df["SEM"],order=2,truncate=True)
ax1.set_xlabel("Mean Rank",fontsize=size_eje)
ax1.set_ylabel("Standard error of the mean",fontsize=size_eje)
ax1.tick_params(labelsize= size_ticks)

#Letter for caption                                               
#----------------------------------------------------------------
Lim_x_up=plt.gca().get_xlim()                                     
Lim_y_up=plt.gca().get_ylim()

ax1.text(Lim_x_up[0] + 0.01*(Lim_x_up[1]-Lim_x_up[0]), Lim_y_up[1] + 0.01*(Lim_y_up[1]-Lim_y_up[0]),Letters[0], fontsize=size_letter,fontweight='bold') 
#----------------------------------------------------------------

#Ranks vs. genomic size
ax2=plt.subplot(gs[0,1])
ax2=sns.regplot(x=df["Mean Rank"], y=df["Size(Mb)"],order=0,truncate=True)
ax2.set_xlabel("Mean Rank",fontsize=size_eje)
ax2.set_ylabel("Genome Size (Mb)",fontsize=size_eje)
ax2.tick_params(labelsize= size_ticks)

#Letter for caption                                               
#----------------------------------------------------------------
Lim_x_up=plt.gca().get_xlim()                                     
Lim_y_up=plt.gca().get_ylim()

ax2.text(Lim_x_up[0] + 0.01*(Lim_x_up[1]-Lim_x_up[0]), Lim_y_up[1] + 0.01*(Lim_y_up[1]-Lim_y_up[0]),Letters[1], fontsize=size_letter,fontweight='bold') 
#----------------------------------------------------------------

plt.savefig('../../Plots_Paper/'+ 'Figure_5_Test.pdf',dpi=300)
plt.show()
#====================================================================
