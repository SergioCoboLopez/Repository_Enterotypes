#1/02/2019. Convertimos las abundancias en 0 o 1 y analizamos luego la nestedness del sistema

import pandas as pd
import time
import random as rnd
from random import shuffle
import math
import numpy as np
from sklearn.model_selection import KFold
import pickle
import sys
import operator
import seaborn as sns
import matplotlib.pyplot as plt
import bisect
import matplotlib.gridspec as gridspec


rnd.seed(int(3333))

#0. FUNCTIONS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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
    for row in range(rows):
        sum_cols[row]=sum(Matrix[row])
    #-------------------------------------------

    #Order the dictionary by ascending sum of columns
    #-------------------------------------------
    sorted_sum_cols = sorted(sum_cols.items(), key=operator.itemgetter(1),reverse=True)
    #-------------------------------------------

    #Reorder matrix rows by the sum of the columns
    #-------------------------------------------
    Dict_row_membership_fun={}
    Matrix_row_ordered=np.zeros(( rows,columns ))
    for new_row in range(rows):
        Matrix_row_ordered[new_row] = Matrix[ sorted_sum_cols[new_row][0] ]
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
    sorted_sum_rows = sorted(sum_rows.items(), key=operator.itemgetter(1),reverse=True)
    #------------------------------------------------ 

    #Reorder matrix rows by the sum of the columns
    #------------------------------------------------
    Dict_column_membership_fun={}
    Matrix_ordered=np.zeros(( rows,columns ))
    for new_col in range(columns):
        np.transpose(Matrix_ordered)[new_col] = \
        np.transpose(Matrix_row_ordered)[ sorted_sum_rows[new_col][0] ]
        Dict_column_membership_fun[new_col]=sorted_sum_rows[new_col][0]
    #------------------------------------------------
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #0.2.2. Compute Nestedness from ordered matrix
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    d_ij_microbes=0;d_ij_patients=0;
    Min_Interactions_microbes=0;Min_Interactions_patients=0;

    for row_i in range(columns-1):
        for row_j in range(row_i+1,columns):

            #Numerator
            #---------------------------------------------------------------------------------------

            #Microbes
            #.......................................................................
            if row_j<rows:
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
            #......................................................................

            #---------------------------------------------------------------------------------------
            
            #Denominator
            #---------------------------------------------------------------------------------------
            if row_j<rows:
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
            
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#1. MAIN
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#1.1. Read Dataset, Crop it and define global variables for matrices
#===================================================================================================

#1.1.1. EXTERNAL PARAMETERS
#---------------------
Dataset=sys.argv[1]#--
#---------------------

AdjacencyMatrix_str=read_file('../../Input_Data/Binarized_Datasets/Binarized_' + Dataset + '.txt')


AdjacencyMatrix=[]
for line in AdjacencyMatrix_str:
    new_line=[float(block) for block in line]
    AdjacencyMatrix.append(new_line)
    
AdjacencyMatrix=np.asarray(AdjacencyMatrix)

global rows;   rows=AdjacencyMatrix.shape[0]
global columns;columns=AdjacencyMatrix.shape[1]
print columns, rows
#===================================================================================================

#1.2. Get Nestedness
#==================================================================
NestedMatrix,Nestedness,Dict_row_membership,Dict_column_membership=\
                                    Nestedness_fun(AdjacencyMatrix)

print NestedMatrix
print "Nestedness",Nestedness, "\n"
#=======================================================================


#1.3. Degree-Distribution-Probabilistic
#================================================================================ 

#1.3.1. Generate matrix of probabilities
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Degree_Probability_Matrix=np.zeros((rows,columns))

for row in range(rows):
    row_prob=sum(NestedMatrix[row])/columns
    for column in range(columns):
        col_prob=sum(np.transpose(NestedMatrix)[column])/rows
        Average_prob=0.5*(row_prob + col_prob)
        Degree_Probability_Matrix[row,column]=Average_prob

print "Degree probability matrix"
print Degree_Probability_Matrix  , '\n'
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#1.3.2. Generate randomized networks
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#Generate many randomized networks and save nestedness
#-------------------------------------------------------------------------------
Nestedness_Vector=[]
for i in range(500):
    Degree_Null_Model,Nestedness_Degree_Null_Model=Null_Model_fun(Degree_Probability_Matrix)
    Nestedness_Vector.append(Nestedness_Degree_Null_Model)
#-------------------------------------------------------------------------------

#Significance of randomized nestedness
#--------------------------------------------------------------------------
counter=0
for element_2 in Nestedness_Vector:
    if element_2 >= Nestedness:
        counter+=1

print "Soft degree distribution randomization (Canonical Link Swap)"
print "pvalue", counter/float(len(Nestedness_Vector))
print "maximum nestedness randomization: ", max(Nestedness_Vector)
print "minimum nestedness randomization: ", min(Nestedness_Vector)
#--------------------------------------------------------------------------

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#================================================================================


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#2. PLOT
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

path_plot= '/export/home/shared/Projects/Microbiome/Plots/Second_Path/Aggregated_Datasets/'


#ax_00=sns.heatmap(NestedMatrix,cmap="RdBu_r",vmin=0,vmax=1)

ax_00=sns.heatmap(NestedMatrix,cmap="Greys",vmin=0,vmax=1)

ax_00.text(10, 12, 'n= %f' %Nestedness,color='white',fontsize=12,
           bbox={'facecolor':'red', 'alpha':0.5, 'pad':2.6})  

#yticks, xticks, ylabels, xlabels                                                                              
#---------------------------------------------
ax_00.set_xlabel("Microbes", fontsize=15)
ax_00.set_ylabel("Patients", fontsize=15)

ax_00.set_xticks([])
ax_00.set_yticks([])
plt.show()



ax_01=sns.distplot(Nestedness_Vector)
ax_01.axvline(x=Nestedness, linewidth=4, color='k')
#ax_01.set_title(Dataset)
ax_01.set_xlabel('Nestedness', fontsize=15)
ax_01.set_ylabel('Density', fontsize=15)
ax_01.text(0.6, 200,"pvalue = " + str(counter/float(len(Nestedness_Vector))), \
         bbox={'facecolor':'red', 'alpha':0.5, 'pad':2.6}, fontsize=12 )
plt.show()


#2.1. Gridspec Parameters                                                                                    
#==============================================================================
fig=plt.figure(figsize=(12, 4.5))
gs = gridspec.GridSpec( 1,2)
gs.update(left=0.07, right=0.95,top=0.95,bottom=0.1, wspace=0.15, hspace=0.2)
#==============================================================================
  

#2.2. Plot Nestedness
#==============================================================================
ax_00 = plt.subplot(gs[0, 0])
ax_00=sns.heatmap(NestedMatrix,cmap="RdBu_r",vmin=0,vmax=1)

ax_00.text(10, 12, 'n= %f' %Nestedness,color='white',fontsize=12,
           bbox={'facecolor':'red', 'alpha':0.5, 'pad':2.6})  

#yticks, xticks, ylabels, xlabels                                                                              
#---------------------------------------------
ax_00.set_xlabel("Microbes", fontsize=15)
ax_00.set_ylabel("Patients", fontsize=15)

ax_00.set_xticks([])
ax_00.set_yticks([])

#---------------------------------------------                                                                 
#==============================================================================


#2.3. Plot Randomization                                                                                     
#==============================================================================                             
ax_01 = plt.subplot(gs[0,1])
ax_01=sns.distplot(Nestedness_Vector)
ax_01.axvline(x=Nestedness, linewidth=4, color='k')
#ax_01.set_title(Dataset)
ax_01.set_xlabel('Nestedness', fontsize=15)
ax_01.set_ylabel('Density', fontsize=15)
ax_01.text(0, 0,"pvalue = " + str(counter/float(len(Nestedness_Vector))), \
         bbox={'facecolor':'red', 'alpha':0.5, 'pad':2.6} )
#==============================================================================


#2.4. Save file
#=======================================================================================
plt.savefig(\
            path_plot + 'Mutualistic_P_Matrix_w_Randomization_'+Dataset+'.png',dpi=300)
plt.show()
#=======================================================================================

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
