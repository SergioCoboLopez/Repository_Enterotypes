#5/07/2019. Codigo para generar la Nestedness de los cinco datasets considerados. Copia del codigo "Nestedeness"
#original para generar graficos en version paper

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
import pandas as pd
import scipy.cluster
from scipy import stats



#0. FUNCTIONS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#0.1. Read Data
#===============================================
def read_file(inFileName):
    lines = open(inFileName).readlines()
    d = [line.strip().split() for line in lines]
    return d
#===============================================

#0.2. Compute Nestedness
#================================================================================================
def Nestedness_fun(Matrix):

    #0.2.1. Reorder matrix by generalist-to-specialist criteria
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #Dictionary of rows: sum of columns
    #-------------------------------------------
    sum_cols={}
    for row in range(K):
        sum_cols[row]=sum(Matrix[row])
    #-------------------------------------------

    #Order the dictionary by ascending sum of columns
    #-------------------------------------------
    sorted_sum_cols = sorted(sum_cols.items(), key=operator.itemgetter(1),reverse=True)
    #-------------------------------------------

    #Reorder matrix rows by the sum of the columns - key: new; value: old
    #----------------------------------------------------------------------
    Dict_row_membership_fun={}
    Matrix_row_ordered=np.zeros(( K,L ))
    for new_row in range(K):
        Matrix_row_ordered[new_row] = Matrix[ sorted_sum_cols[new_row][0] ]
        Dict_row_membership_fun[new_row]=sorted_sum_cols[new_row][0]
    #----------------------------------------------------------------------
        
    #Dictionary of cols: sum of rows
    #-------------------------------------------
    sum_rows={}
    for col in range(L):
        sum_rows[col]=sum( np.transpose(Matrix_row_ordered)[col] )
    #-------------------------------------------

    #Order the dictionary by ascending sum of rows
    #------------------------------------------------
    sorted_sum_rows = sorted(sum_rows.items(), key=operator.itemgetter(1),reverse=True)
    #------------------------------------------------

    #Reorder matrix rows by the sum of the columns
    #------------------------------------------------
    Dict_column_membership_fun={}
    Matrix_ordered=np.zeros(( K,L ))
    for new_col in range(L):
        np.transpose(Matrix_ordered)[new_col] = \
        np.transpose(Matrix_row_ordered)[ sorted_sum_rows[new_col][0] ]
        Dict_column_membership_fun[new_col]=sorted_sum_rows[new_col][0]
    #------------------------------------------------
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #0.2.2. Compute Nestedness from ordered matrix
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    d_ij_microbes=0;d_ij_patients=0;
    Min_Interactions_microbes=0;Min_Interactions_patients=0;

    for row_i in range(L-1):
        for row_j in range(row_i+1,L):

            #Numerator
            #---------------------------------------------------------------------------------------
            #Microbes
            #.......................................................................
            if row_j<K:
                product_microbes=np.multiply(Matrix_ordered[row_i], Matrix_ordered[row_j])
                shared_interactions_microbes=np.sum( product_microbes)
                d_ij_microbes+=shared_interactions_microbes
            #.......................................................................

            #Patients
            #.......................................................................
            product_patients=np.multiply(np.transpose(Matrix_ordered)[row_i],\
                                         np.transpose(Matrix_ordered)[row_j])
            shared_interactions_patients=np.sum( product_patients)
            d_ij_patients+=shared_interactions_patients
            #.......................................................................

            #---------------------------------------------------------------------------------------

            #Denominator
            #---------------------------------------------------------------------------------------
            if row_j<K:
                Interactions_microbes_i=np.sum(Matrix_ordered[row_i])
                Interactions_microbes_j=np.sum(Matrix_ordered[row_j])
                Min_Interactions_microbes+=min(Interactions_microbes_i, Interactions_microbes_j)

            Interactions_patients_i=np.sum(np.transpose(Matrix_ordered)[row_i])
            Interactions_patients_j=np.sum(np.transpose(Matrix_ordered)[row_j])
            Min_Interactions_patients+=min(Interactions_patients_i, Interactions_patients_j)
            #---------------------------------------------------------------------------------------
    Nestedness=float(d_ij_microbes+d_ij_patients)/float(Min_Interactions_microbes+Min_Interactions_patients)

    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
    return Matrix_ordered, Nestedness, Dict_row_membership_fun, Dict_column_membership_fun
#================================================================================================

#0.3. Randomized Matrix (Hard Null Model)
#==============================================================================
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
#==============================================================================

#0.4. Shannon Entropy
#===========================================
def Shannon(v):
    Entropy=0
    for i in v:
        if i==0:
            Entropy+=0
        else:
            Entropy+=-i*math.log(i,len(v))

    return Entropy
#===========================================

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#1. PARAMETERS FOR INPUT/OUTPUT DATA
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#1.1. EXTERNAL PARAMETERS
#=====================
number_of_links=3  #==
rnd.seed(1)        #==
#=====================

#1.2. INTERNAL PARAMETERS
#===========================================================================================================

#::::::::::::::::::::::::::::::::::::
K=10;L=20
Microbe_Entropies={}
pMatrices={};All_Nestedness={}
Triangles={};Squares={}
Null_Models={}
#::::::::::::::::::::::::::::::::::::


Datasets=['S-8_Tot', 'V-10_Tot', 'V-22_Tot', 'V-23_24_Tot', 'V-25_Tot']

#New names for title-lables
New_Names={'S-8_Tot': 'Liu-2016 (S-8)', 'V-10_Tot': 'Qin-2014 (V-10)', 'V-22_Tot': 'Schirmer-2016 (V-22)',\
        'V-23_24_Tot': 'C. Huttenhower-2012 & Lloyd-Price J-2017 (V-23_24)','V-25_Tot': 'Zeevi-2015 (V-25)'}

Patients_per_Dataset=\
{'S-8_Tot':107, 'V-10_Tot':92, 'V-22_Tot':467, 'V-23_24_Tot':222,'V-25_Tot':883}

Microbes_per_Dataset=\
{'S-8_Tot':128, 'V-10_Tot':137, 'V-22_Tot':134, 'V-23_24_Tot':118,'V-25_Tot':144}
#===========================================================================================================

for Dataset in Datasets:

    Patients=Patients_per_Dataset[Dataset]
    Microbes=Microbes_per_Dataset[Dataset]

    #1.3. PATHS FOR FILES
    #==============================================================
    Path_out='../../Output_Data/'+ 'Leave_One_Out_' + Dataset + '/'
    Path_in= '../../Input_Data/' + 'Leave_One_Out_' + Dataset + '/'
    path_plot='/export/home/shared/Projects/Microbiome/Plots_Paper/'
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

    filename='Dataset_' + Dataset + '_Seed_' + str(Seed) + '_Parameters.txt'
    Data=read_file(Path_out+filename)    
    #======================================================================

    #1.5. READ PMATRICES AND FORMAT INTO MATRICES
    #======================================================================
    counter=0
    pMatrix=np.zeros((K,L))
    Nodes_Links_Dict={};Links_Nodes_Dict={}
    for line in Data:
        if 'slice' and '0]' in line:
            for row in range(K):
                for col in range(L):

                    Prob_of_zero=float(Data[counter+row+1][col])
                    Weighted_Link=1-Prob_of_zero

                    Node_Pair=(row,col)

                    #Binarize matrix elements
                    #::::::::::::::::::::::::::::::::: 
                    if Weighted_Link>=0.5:
                        Link=1
                    else:
                        Link=0
                    #:::::::::::::::::::::::::::::::::

                    Nodes_Links_Dict[Node_Pair]=Link

                    try:
                        Links_Nodes_Dict[Link].append(Node_Pair)
                    except KeyError:
                            Links_Nodes_Dict[Link]=[Node_Pair]

                    pMatrix[row][col]=Link

        counter+=1

    pMatrix_ordered,Nestedness_pMatrix,Dict_row_membership,Dict_column_membership=Nestedness_fun(pMatrix)

    
    pMatrices[Dataset]=pMatrix_ordered
    Triangles[Dataset]=Dict_column_membership
    Squares[Dataset]=Dict_row_membership
    All_Nestedness[Dataset]=Nestedness_pMatrix
    #=======================================================================


    #1.1.1. Taxonomic Information
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    File_Taxon='Names_Microbes_' + Dataset + '.csv'
    Taxon_Data_Raw=pd.read_csv('../../Input_Data/' + File_Taxon,header=None)
    Species_Names=Taxon_Data_Raw[6]
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #Synthetic ranks for latent space. From zero (generalist) to 1 (specialist)
    #::::::::::::::::::::::::::::::::::::::::::
    Latent_Ranks=np.linspace(0, 1, num=20)
    #::::::::::::::::::::::::::::::::::::::::::


    #Eta, membership vectors of microbes
    #========================================
    etas=Data[Patients+1:Patients+1+Microbes]
    for tupla in zip(Species_Names,etas):
        name=tupla[0];eta=tupla[1]
        eta=[float(element) for element in eta]

        #Reorder eta according to nested version
        #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        ordered_eta=[None]*L
        for col_element in Dict_column_membership:
            ordered_eta[ col_element ]=eta[ Dict_column_membership[col_element] ]
        #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


        #Compute rank of vector according to membership weight multiplied by synthetic ranks
        #::::::::::::::::::::::::::::::::::::::::::::::::::::
        print ordered_eta, "\n"
        Rank=0
        for position in zip(ordered_eta, Latent_Ranks):
            Rank+=position[0]*position[1]
        #::::::::::::::::::::::::::::::::::::::::::::::::::::
            
        try:
            Microbe_Entropies[name][Dataset]=Rank
        except KeyError:
            Microbe_Entropies[name]={Dataset:Rank}
    #========================================

    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
print Microbe_Entropies

Microbes_sorted_means={}
Rankings=[]
for microbe in Microbe_Entropies:

    if len(Microbe_Entropies[microbe])>4:
        Rankings.append( (np.mean(Microbe_Entropies[microbe].values()),\
                          scipy.stats.sem(Microbe_Entropies[microbe].values(), ddof=0),\
                          np.std(Microbe_Entropies[microbe].values())) ) 

        Microbes_sorted_means[microbe]={'mean':np.mean(Microbe_Entropies[microbe].values()),\
                                        'SEM': scipy.stats.sem(Microbe_Entropies[microbe].values(),ddof=0),\
                                        'std':np.std(Microbe_Entropies[microbe].values())}

sorted_by_first = sorted(Rankings, key=lambda tup: tup[0])

Means=[]; SEMs=[]; stds=[]
for pair in sorted_by_first:
    print pair
    Means.append(pair[0])
    SEMs.append(pair[1])
    stds.append(pair[2])


#3. PLOTS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
plt.xlabel("Means", fontsize=14)
plt.ylabel("SEMs", fontsize=14)

ax = sns.regplot(Means, SEMs, order=2,truncate=True)#, height=5)

plt.savefig('/export/home/shared/Projects/Microbiome/Plots_Paper/Comparison_Ranks_Latent_5.png',dpi=300)
plt.show()
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
