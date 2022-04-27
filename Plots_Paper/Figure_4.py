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

#0. FUNCTIONS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    d_ij_microbes=0;d_ij_patients=0;
    Min_Interactions_microbes=0;Min_Interactions_patients=0;

    Index=max(rows, columns)

    for row_i in range(Index-1):
        for row_j in range(row_i+1,Index):

            #Numerator
            #-------------------------------------------------------

            #Patients 
            #.......................................................
            if row_j<rows:
                product_patients=np.multiply(Matrix_ordered[row_i], Matrix_ordered[row_j])
                shared_interactions_patients=np.sum( product_patients)
                d_ij_patients+=shared_interactions_patients
            #.......................................................

            #Microbes
            #.......................................................
            if row_j<columns:
                product_microbes=np.multiply(np.transpose(Matrix_ordered)[row_i],\
                                             np.transpose(Matrix_ordered)[row_j])
                shared_interactions_microbes=np.sum( product_microbes)
                d_ij_microbes+=shared_interactions_microbes
            #.......................................................

            #-------------------------------------------------------

            #Denominator
            #-------------------------------------------------------

            #Patients
            #.......................................................
            if row_j<rows:
                Interactions_patients_i=np.sum(Matrix_ordered[row_i])
                Interactions_patients_j=np.sum(Matrix_ordered[row_j])
                Min_Interactions_patients+=min(Interactions_patients_i, Interactions_patients_j)
            #.......................................................

            #Microbes
            #.......................................................
            if row_j<columns:
                Interactions_microbes_i=np.sum(np.transpose(Matrix_ordered)[row_i])
                Interactions_microbes_j=np.sum(np.transpose(Matrix_ordered)[row_j])
                Min_Interactions_microbes+=min(Interactions_microbes_i, Interactions_microbes_j)
            #.......................................................

            #-------------------------------------------------------

    Nestedness=float(d_ij_microbes+d_ij_patients)/float(Min_Interactions_microbes+Min_Interactions_patients)
    
    return Matrix_ordered, Nestedness, Dict_row_membership_fun, Dict_column_membership_fun
#===================================================================

#0.3. Randomized Matrix (Hard Null Model)
#===================================================================
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
#===================================================================

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#1. PARAMETERS FOR INPUT/OUTPUT DATA
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#1.1. EXTERNAL PARAMETERS
#=====================
number_of_links=3  #==
rnd.seed(1)        #==
#=====================

#1.2. INTERNAL PARAMETERS
#===================================================================

#::::::::::::::::::::::::::::::::::::
K=10;L=20
Nested_Matrices={};All_Nestedness={}
Triangles={};Squares={}
Null_Models={}
#::::::::::::::::::::::::::::::::::::


Datasets=['S-8_Tot','V-10_Tot','V-22_Tot','V-23_24_Tot','V-25_Tot']

Patients_per_Dataset={'S-8_Tot':107,'V-10_Tot':92,\
    'V-22_Tot':467,'V-23_24_Tot':222,'V-25_Tot':883}

Microbes_per_Dataset={'S-8_Tot':128,'V-10_Tot':137,\
    'V-22_Tot':134, 'V-23_24_Tot':118,'V-25_Tot':144}
#===================================================================

for Dataset in Datasets:
    Patients=Patients_per_Dataset[Dataset]
    Microbes=Microbes_per_Dataset[Dataset]

    AdjacencyMatrix_str=read_file('../../Input_Data/Binarized_Datasets/Binarized_' + Dataset + '_Threshold_5e-6.txt')

    
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

    print Dataset
    print rows, columns
    #===============================================================

    #1.2. Get Nestedness 
    #===============================================================
    NestedMatrix,Nestedness,Dict_row_membership,Dict_column_membership=\
                                    Nestedness_fun(AdjacencyMatrix)

    print NestedMatrix
    print "Nestedness",Nestedness, "\n"

    Nested_Matrices[Dataset]=NestedMatrix
    All_Nestedness[Dataset]=Nestedness
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    #2. NULL MODELS
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    #:::::::::::::::::::::::::::::::::::::::::::::::::::::
    np.random.seed(seed=1)
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::

    #===============================================================

    #2.3. Degree-Distribution-Probabilistic
    #===============================================================

    #2.3.1. Generate matrix of probabilities
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    Degree_Probability_Matrix=np.zeros((rows,columns))

    for row in range(rows):
        row_prob=sum(NestedMatrix[row])/columns
        for column in range(columns):
            col_prob=sum(np.transpose(NestedMatrix)[column])/rows
            Average_prob=0.5*(row_prob + col_prob)
            Degree_Probability_Matrix[row,column]=Average_prob
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #2.3.2. Generate randomized networks
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #Generate many randomized networks and save nestedness
    #---------------------------------------------------------------
    Nestedness_Vector=[]
    for i in range(119):
#        print(i)
        Degree_Null_Model,Nestedness_Degree_Null_Model=Null_Model_fun(Degree_Probability_Matrix)
        Nestedness_Vector.append(Nestedness_Degree_Null_Model)
    #---------------------------------------------------------------

    #Significance of randomized nestedness
    #------------------------------------------
    counter=0
    for element_2 in Nestedness_Vector:
        if element_2 >= Nestedness:
            counter+=1
            print("HOLA!")
            print(counter)

    print(len(Nestedness_Vector))

    Null_Models[Dataset]={'vector': Nestedness_Vector, 'p-value': (counter+1)/float(len(Nestedness_Vector))}
    #------------------------------------------

    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #================================================================

    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#3. PLOTS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#4.1.1. Fontsizes and Parameters
#=================
size_eje=14
size_ticks=11
size_title=12
size_letter=15
Letters=['a','b','c','d','e']
#=================


#3.2.1. Gridspec Parameters
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#Monitor Grande
#--------------------------------------------------------------------
#fig=plt.figure(figsize=(19, 7))
#--------------------------------------------------------------------

#TXIKIMONITOR
#--------------------------------------------------------------------
fig=plt.figure(figsize=(13.5, 5))

#Monitor Grande
#--------------------------------------------------------------------
# gs = gridspec.GridSpec(2, 6, width_ratios=[1,1,1,1,1,0.09])
# gs.update(left=0.05,right=0.96,bottom=0.14,top=0.88,wspace=0.3,hspace=0.25)
#--------------------------------------------------------------------

#TXIKIMONITOR
#--------------------------------------------------------------------
gs = gridspec.GridSpec(2, 6, width_ratios=[1,1,1,1,1,0.09])
gs.update(left=0.05,right=0.94,bottom=0.14,top=0.88,wspace=0.3,hspace=0.25)
#--------------------------------------------------------------------

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#(0,6) PLOT COLORBAR
#====================================================================
ax_6 = plt.subplot(gs[0, 5])
colbar_ticks=[0,0.5,1]
cb1 = matplotlib.colorbar.ColorbarBase(ax_6, cmap="Blues",ticks=colbar_ticks)
cb1.set_ticks(colbar_ticks)
cb1.set_ticklabels(colbar_ticks)
cb1.ax.tick_params(labelsize=size_ticks)
cb1.outline.set_visible(False)
cb1.set_label("p(Non-negligible)",fontsize=12)
#====================================================================

#2. PLOT                
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
path_plot= '../../Plots_Paper/'

Counter_Plot=0

for Dataset in Datasets:
    #2.2. Plot Nestedness
    #================================================================
    ax_00 = plt.subplot(gs[0, Counter_Plot])
    ax_00=sns.heatmap(Nested_Matrices[Dataset],cmap="Blues",vmin=0,vmax=1,cbar=False)

    Lim_x_up=plt.gca().get_xlim()
    Lim_y_up=plt.gca().get_ylim()

#Monitor grande
#--------------------------------------------------------------------
    # ax_00.text(0.67*Lim_x_up[1], 0.9*Lim_y_up[0], 'n= %.3f' %All_Nestedness[Dataset],color='black',fontsize=12, weight='bold')
#--------------------------------------------------------------------

#TXIKIMONITOR
#--------------------------------------------------------------------
    ax_00.text(0.47*Lim_x_up[1], 0.9*Lim_y_up[0], 'n= %.3f' %All_Nestedness[Dataset],color='black',fontsize=12, weight='bold')
#--------------------------------------------------------------------


    #yticks, xticks, ylabels, xlabels
    #---------------------------------------------
    ax_00.set_xlabel("Microbes", fontsize=15)
    if Counter_Plot==0:
        ax_00.set_ylabel("Hosts", fontsize=15)

    ax_00.set_xticks([])
    ax_00.set_yticks([])
    #---------------------------------------------

    #Letter for caption                                          
    #----------------------------------------------------------------
    ax_00.text(Lim_x_up[0]-0.1*(Lim_x_up[1]-Lim_x_up[0]),Lim_y_up[1] + \
    0.1*(Lim_y_up[1]-Lim_y_up[0]), Letters[Counter_Plot], fontsize=size_letter,fontweight='bold')
    #----------------------------------------------------------------

    #Dibujar cuadrado/marco                                        
    #----------------------------------------------
    autoAxis0 = ax_00.axis()
    # rec0 = Rectangle(\
    # (autoAxis0[0]-0,autoAxis0[2]-10),\
    # (autoAxis0[1]-autoAxis0[0])+19,(autoAxis0[3]-autoAxis0[2])+9,\
    #         fill=False,lw=1.5,color='gray')

    rec0 = Rectangle(\
    (autoAxis0[0]-0,autoAxis0[2]),\
    (autoAxis0[1]-autoAxis0[0]),(autoAxis0[3]-autoAxis0[2]),\
            fill=False,lw=1.5,color='gray')
        
    rec0 = ax_00.add_patch(rec0)
    rec0.set_clip_on(False)                                         
    #----------------------------------------------
    #================================================================

    #2.3. Plot Randomization
    #================================================================
    ax_01 = plt.subplot(gs[1,Counter_Plot])
    ax_01=sns.distplot(Null_Models[Dataset]['vector'])
    ax_01.axvline(x=All_Nestedness[Dataset], linewidth=4, color='k')
    ax_01.set_xlabel('Nestedness',fontsize=15)

    Lim_x_up_1=plt.gca().get_xlim()
    Lim_y_up_1=plt.gca().get_ylim()

    print Lim_x_up_1
    print 0.8*Lim_x_up_1[1]


    if Counter_Plot==0:
         ax_01.set_ylabel('Density',fontsize=15)
        


    ax_01.text(1.01*Lim_x_up_1[0], 0.9*Lim_y_up_1[1],"pvalue < 0.001", fontsize=12 )

    #ax_01.text(1.2*Lim_x_up_1[0], 0.9*Lim_y_up_1[1],"pvalue < "+ str(Null_Models[Dataset]['p-value']), fontsize=12)
    #================================================================
    Counter_Plot+=1

plt.savefig(path_plot + 'Figure_4_Threshold_5e-6.pdf',dpi=300)
plt.show()
