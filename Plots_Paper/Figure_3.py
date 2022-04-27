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
#====================================================================
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

    #Reorder matrix rows by the sum of the columns
    #-------------------------------------------
    Dict_row_membership_fun={}
    Matrix_row_ordered=np.zeros(( K,L ))
    for new_row in range(K):
        Matrix_row_ordered[new_row] = Matrix[ sorted_sum_cols[new_row][0] ]
        Dict_row_membership_fun[new_row]=sorted_sum_cols[new_row][0]
    #-------------------------------------------
        
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
            #--------------------------------------------------------
            #Microbes
            #........................................................
            if row_j<K:
                product_microbes=np.multiply(Matrix_ordered[row_i], Matrix_ordered[row_j])
                shared_interactions_microbes=np.sum( product_microbes)
                d_ij_microbes+=shared_interactions_microbes
            #........................................................

            #Patients
            #........................................................
            product_patients=np.multiply(np.transpose(Matrix_ordered)[row_i],\
                                         np.transpose(Matrix_ordered)[row_j])
            shared_interactions_patients=np.sum( product_patients)
            d_ij_patients+=shared_interactions_patients
            #........................................................

            #--------------------------------------------------------

            #Denominator
            #--------------------------------------------------------
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

#1.2. INTERNAL PARAMETERS
#====================================================================

#::::::::::::::::::::::::::::::::::::
K=10;L=20
pMatrices={};All_Nestedness={}
Triangles={};Squares={}
Null_Models={}
#::::::::::::::::::::::::::::::::::::


Datasets=['S-8_Tot', 'V-10_Tot', 'V-22_Tot', 'V-23_24_Tot', 'V-25_Tot']

#New names for title-labels
New_Names={'S-8_Tot': 'Liu-2016 (S-8)', 'V-10_Tot': 'Qin-2014 (V-10)', 'V-22_Tot': 'Schirmer-2016 (V-22)',\
        'V-23_24_Tot': 'C. Huttenhower-2012 & Lloyd-Price J-2017 (V-23_24)','V-25_Tot': 'Zeevi-2015 (V-25)'}

Patients_per_Dataset=\
{'S-8_Tot':107, 'V-10_Tot':92, 'V-22_Tot':467, 'V-23_24_Tot':222,'V-25_Tot':883}

Microbes_per_Dataset=\
{'S-8_Tot':128, 'V-10_Tot':137, 'V-22_Tot':134, 'V-23_24_Tot':118,'V-25_Tot':144}
#====================================================================

for Dataset in Datasets:

    Patients=Patients_per_Dataset[Dataset]
    Microbes=Microbes_per_Dataset[Dataset]

    #1.3. PATHS FOR FILES
    #==============================================================
    Path_out='../../Output_Data/'+ 'Leave_One_Out_' + Dataset + '/'
    Path_in= '../../Input_Data/' + 'Leave_One_Out_' + Dataset + '/'
    path_plot='../../Plots_Paper/'
    #==============================================================

    #1.4. CHOOSE BEST LIKELIHOOD
    #===============================================================
    Likelihood_Seed={}
    for Seed in range(1,10,2):
        likelihood_file=\
        'Dataset_'+Dataset+'_Seed_'+str(Seed) + '_LogLikelihood.txt'
        likelihoodfinal=read_file(Path_out + likelihood_file)[-1]
        Likelihood_Seed[float(likelihoodfinal[0])]=Seed

    print Likelihood_Seed
    print max(Likelihood_Seed)
    print Likelihood_Seed[max(Likelihood_Seed)]
    Seed=Likelihood_Seed[max(Likelihood_Seed)]

    filename='Dataset_' + Dataset + '_Seed_' + str(Seed) + '_Parameters.txt'
    Data=read_file(Path_out+filename)    
    #===============================================================

    #1.5. READ PMATRICES AND FORMAT INTO MATRICES
    #===============================================================
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
    #===============================================================

    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    #2. NULL MODELS
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    rows=np.shape(pMatrix_ordered)[0]
    cols=np.shape(pMatrix_ordered)[1]

    #2.1. Non-probabilistic version (microcanonical Erdos-Renyi)
    #================================================================================

    #2.1.1. Reshape matrix to vector
    #::::::::::::::::::::::::::::::::::::
    pMatrix_1D=np.reshape(pMatrix,rows*cols)
    #::::::::::::::::::::::::::::::::::::

    #:::::::::::::::::::::::::::::::::::::::::::::::::::::
    np.random.seed(seed=1)
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::

    #================================================================================

    #2.3. Degree-Distribution-Probabilistic
    #================================================================================

    #2.3.1. Generate matrix of probabilities
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    Degree_Probability_Matrix=np.zeros((rows,cols))

    for row in range(rows):
        row_prob=sum(pMatrix_ordered[row])/cols
        for column in range(cols):
            col_prob=sum(np.transpose(pMatrix_ordered)[column])/rows
            Average_prob=0.5*(row_prob + col_prob)
            Degree_Probability_Matrix[row,column]=Average_prob
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #2.3.2. Generate randomized networks
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #Generate many randomized networks and save nestedness
    #-------------------------------------------------------------------------------
    Nestedness_Vector=[]
    for i in range(1000):
        Degree_Null_Model,Nestedness_Degree_Null_Model=Null_Model_fun(Degree_Probability_Matrix)
        Nestedness_Vector.append(Nestedness_Degree_Null_Model)
    #-------------------------------------------------------------------------------

    #Significance of randomized nestedness
    #------------------------------------------
    counter=0
    for element_2 in Nestedness_Vector:
        if element_2 >= Nestedness_pMatrix:
            counter+=1

    Null_Models[Dataset]={'vector': Nestedness_Vector, 'p-value': (counter+1)/float(len(Nestedness_Vector))}
    #------------------------------------------
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    #================================================================================

    #2.4. Node swapping
    #================================================================================

    for iteration in range(200):

        #Propose Tuples for swapping
        #========================================================================
        Removed_Tuple_0, Removed_Tuple_1=rnd.sample(Links_Nodes_Dict[1],2)

        Tentative_Tuple_0=( Removed_Tuple_0[0],Removed_Tuple_1[1] )
        Tentative_Tuple_1=( Removed_Tuple_1[0],Removed_Tuple_0[1] )

        counter_while=0
        while Nodes_Links_Dict[Tentative_Tuple_0]==1 \
              or Nodes_Links_Dict[Tentative_Tuple_1]==1:
            Removed_Tuple_0,Removed_Tuple_1=rnd.sample(Links_Nodes_Dict[1],2)
            Tentative_Tuple_0=( Removed_Tuple_0[0],Removed_Tuple_1[1] )
            Tentative_Tuple_1=( Removed_Tuple_1[0],Removed_Tuple_0[1] )
            counter_while+=1

            if counter_while>100:
                print "No encontramos ninguna tupla compatible"
                break

        if counter_while>100:
            break

        #========================================================================

        #Do swapping
        #========================================================================
        New_Tuple_0=Tentative_Tuple_0;New_Tuple_1=Tentative_Tuple_1

        #Removed Tuples: update dictionaries
        #::::::::::::::::::::::::::::::::::::::::::::
        Links_Nodes_Dict[0].append(Removed_Tuple_0)
        Links_Nodes_Dict[0].append(Removed_Tuple_1)

        Links_Nodes_Dict[1].remove(Removed_Tuple_0)
        Links_Nodes_Dict[1].remove(Removed_Tuple_1)

        Nodes_Links_Dict[Removed_Tuple_0]=0;Nodes_Links_Dict[Removed_Tuple_1]=0
        #::::::::::::::::::::::::::::::::::::::::::::

        #New Tuples: Update dictionaries
        #::::::::::::::::::::::::::::::::::::::::::::
        Links_Nodes_Dict[0].remove(New_Tuple_0)
        Links_Nodes_Dict[0].remove(New_Tuple_1)

        Links_Nodes_Dict[1].append(New_Tuple_0)
        Links_Nodes_Dict[1].append(New_Tuple_1)

        Nodes_Links_Dict[New_Tuple_0]=1;Nodes_Links_Dict[New_Tuple_1]=1
        #::::::::::::::::::::::::::::::::::::::::::::


    pMatrix_Randomized=np.zeros((K,L))

    for row in range(K):
        for col in range(L):
            Tuple=(row,col)
            pMatrix_Randomized[row][col]=Nodes_Links_Dict[Tuple]

    pMatrix_Rand_ordered,Nestedness_2,Dict_row_membership_Rand,Dict_col_membership_Rand=\
    Nestedness_fun(Degree_Probability_Matrix)

    #================================================================
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#3. PLOTS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Information=["Only_Nestedness", "Nestedness_w_Groups", "Randomization_Nestedness"]
Plot=Information[1]

#Choose Colormap
colmap="Blues" #"PuBu"  #"RdBu_r"


# elif Plot=="Only_Nestedness":

#     #Only Nestedness
#     #================================================================
#     ax=sns.heatmap(pMatrix_ordered,cmap=colmap,linewidths=1,vmin=0,vmax=1,xticklabels=False,yticklabels=False, \
#                    cbar_kws={'label': 'P(Connection)'})

#     ax.text(0.75, 0.75, 'nestedness= %f' %Nestedness_pMatrix,color='white',fontsize=12,
#             bbox={'facecolor':'red', 'alpha':0.5, 'pad':2.6})

#     #---------------------------------------------
# #    ax.set_title("Nestedness,"+Dataset)
#     # ax.set_xlabel("Groups of Microbes",fontsize=12)
#     # ax.set_ylabel("Groups of Patients",fontsize=12)
#     ax.set_xlabel("Groups of Microbes",fontsize=12)
#     ax.set_ylabel("Enterotypes",fontsize=12)
#     #---------------------------------------------
#     plt.savefig(\
#     path_plot + 'Mutualistic_P_Matrix_'+Dataset+'.png',dpi=300)
#     plt.show()
#     #================================================================

if Plot=="Nestedness_w_Groups":
    print "HOLA"
    #================================================================

    #4.1.1. Fontsizes and Parameters
    #=================
    size_eje=14
    size_ticks=11
    size_title=12
    size_letter=15
    Letters=['a','b','c','d','e']
    #=================


    #3.2.0. Shapes and colors for groups of microbes/patients
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    colors_shapes_microbes={}
    shapes=["^", "v"]
    counter=0

    for shape in range(len(shapes)):
        for group in range(L/2):
            colors_shapes_microbes[counter]=( group, shapes[shape] )
            counter+=1
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #3.2.1. Gridspec Parameters
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    fig=plt.figure(figsize=(19, 7))

    gs = gridspec.GridSpec(2, 6, width_ratios=[1,1,1,1,1,0.09])
    gs.update(left=0.012,right=0.96,bottom=0.14,top=0.88,wspace=0.09,hspace=0.25)
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #(0,6) PLOT COLORBAR
    #================================================================
    ax6 = gridspec.GridSpecFromSubplotSpec(\
        2, 2, height_ratios=[0.04, 1], width_ratios=[0.5, 1.1],subplot_spec=gs[5] )

    colorbar=plt.subplot(ax6[1,1])
    colbar_ticks=[0,0.5,1]
    cb1 = matplotlib.colorbar.ColorbarBase(colorbar, cmap=colmap,ticks=colbar_ticks)
    cb1.set_ticks(colbar_ticks)
    cb1.set_ticklabels(colbar_ticks)
    cb1.ax.tick_params(labelsize=size_ticks)
    cb1.outline.set_visible(False)
    cb1.set_label("p(Connection)",fontsize=15)
    #================================================================

    #(0,0:5) NESTEDNESS PLOTS
    #================================================================
    counter_plot=0
    for Dataset in Datasets:
        #-------------------------------------------------------------------------------
        gs_sub =gridspec.GridSpecFromSubplotSpec(3, 3, height_ratios=[0, 0, 1],\
                        width_ratios=[0, 0, 1],subplot_spec=gs[counter_plot] )
        #-------------------------------------------------------------------------------
        
        #Letter for caption
        #------------------------------------------------------------
        Lim_x_up=plt.gca().get_xlim()
        Lim_y_up=plt.gca().get_ylim()

        # squares_left.text(Lim_x_up[0]+0.4*(Lim_x_up[1]-Lim_x_up[0]),Lim_y_up[1]+ 0.2*(Lim_y_up[1]\
        #                 -Lim_y_up[0]),Letters[counter_plot], fontsize=size_letter,fontweight='bold')
        #------------------------------------------------------------


        #3.2.3. Plot Mutualistic Network
        #------------------------------------------------------------
        ax = plt.subplot(gs_sub[2, 2])

        
        #Dibujar cuadrado/marco
        #----------------------------------------------    
        autoAxis0 = ax.axis()

        rec0 = Rectangle(\
                 (autoAxis0[0]-0,autoAxis0[2]-0.0),\
                 (autoAxis0[1]-autoAxis0[0])+19,(autoAxis0[3]-autoAxis0[2])+9,\
                 fill=False,lw=1.5,color='gray')
        
        rec0 = ax.add_patch(rec0)
        rec0.set_clip_on(False)
        #----------------------------------------------

        ax=sns.heatmap(pMatrices[Dataset],cmap=colmap,linewidths=1,vmin=0,vmax=1,cbar=False)
        
	Lim_x_up_main=plt.gca().get_xlim()
	Lim_y_up_main=plt.gca().get_ylim()

        ax.text(Lim_x_up_main[0] - 2,Lim_y_up_main[1] - 1.25 ,Letters[counter_plot], fontsize=size_letter,fontweight='bold')
        # ax.text(Lim_x_up[0]+0.4*(Lim_x_up[1]-Lim_x_up[0]),Lim_y_up[1]+ 0.2*(Lim_y_up[1]\
        #     -Lim_y_up[0]),Letters[counter_plot], fontsize=size_letter,fontweight='bold')
        

        print Lim_y_up_main
        print(Letters[counter_plot])

        ax.text(14,0.9*Lim_y_up_main[0],'n= %.3f' %All_Nestedness[Dataset],\
                color='black',fontsize=12,weight='bold')

        #yticks, xticks, ylabels, xlabels
        #-----------------------------------------------------------------
        # yticks_pos=[row + 0.5 for row in range(K)]
        # yticks_labels=Squares[Dataset].values()
        # ax.set_yticks(yticks_pos)
        # ax.set_yticklabels(yticks_labels,fontsize=size_ticks)

        # xticks_pos=[col + 0.5 for col in range(L)]
        # xticks_labels=Triangles[Dataset].values()
        # ax.set_xticks(xticks_pos)
        # ax.set_xticklabels(xticks_labels,\
        # fontsize=size_ticks,rotation='90')

        # ax.tick_params(axis=u'both', which=u'both',length=0)


        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel("Groups of Microbes",fontsize=size_eje)
        if counter_plot==0:
            ax.set_ylabel("Groups of Hosts",fontsize=size_eje,labelpad=16)
        #------------------------------------------------------------------

        gs_rnd =gridspec.GridSpecFromSubplotSpec(3, 3, height_ratios=[0.11, -0.1, 1], \
                        width_ratios=[0.0, -0.05, 1],subplot_spec=gs[1,counter_plot] )

        #3.2.3. Plot Randomization  
        #-----------------------------------------------------------
        ax_Rnd = plt.subplot(gs_rnd[2, 2])
        ax_Rnd = sns.distplot(Null_Models[Dataset]['vector'])
#        plt.axvline(x=Nestedness_pMatrix, linewidth=4, color='k')
        plt.axvline(x=All_Nestedness[Dataset], linewidth=4, color='k')

        ax_Rnd.set_xlim([0.3,0.9])
        Lim_x_up_1=plt.gca().get_xlim()
        Lim_y_up_1=plt.gca().get_ylim()

        plt.xlabel('Nestedness',fontsize=size_eje)
        if counter_plot==0:
            plt.ylabel('Density',fontsize=size_eje)
        
        plt.text(0.54, 0.9*Lim_y_up_1[1] ,"pvalue < " + str(Null_Models[Dataset]['p-value']), fontsize=12)
        print("hola!")
        print(All_Nestedness[Dataset])

        # randomization=plt.subplot(gs[1, counter_plot])
        # randomization=sns.distplot(Null_Models[Dataset]['vector'])
        # plt.axvline(x=Nestedness_pMatrix, linewidth=4, color='k')
        # randomization.set_xlim([0.3,0.9])
        # Lim_x_up_1=plt.gca().get_xlim()
        # Lim_y_up_1=plt.gca().get_ylim()
        # plt.xlabel('Nestedness',fontsize=size_eje)
        if counter_plot==0:
            plt.ylabel('Density',fontsize=size_eje)


        # plt.text(0.54, 0.9*Lim_y_up_1[1] ,"pvalue < " + str(Null_Models[Dataset]['p-value']), fontsize=12)
        #------------------------------------------------------------

        counter_plot+=1
        #------------------------------------------------------------

    #================================================================

    plt.savefig(path_plot+'Figure_3.pdf',dpi=300)
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    plt.show()
    #================================================================

elif Plot=="Randomization_Nestedness":
    #================================================================
    sns.distplot(Nestedness_Vector)
    plt.axvline(x=Nestedness_pMatrix, linewidth=4, color='k')

    plt.xlabel('Nestedness')
    plt.ylabel('Density')
    pvalue=(counter+1)/float(len(Nestedness_Vector))
    plt.text(min(Nestedness_Vector), 8,"pvalue = %.3f" %pvalue, \
             bbox={'facecolor':'red', 'alpha':0.5, 'pad':2.6} )

    plt.show()
    #================================================================


elif Plot=="Only_Nestedness":

    #Only Nestedness
    #================================================================
    ax=sns.heatmap(pMatrix_ordered,cmap=colmap,linewidths=1,vmin=0,vmax=1,xticklabels=False,yticklabels=False, \
                   cbar_kws={'label': 'P(Connection)'})

    ax.text(0.75, 0.75, 'nestedness= %f' %Nestedness_pMatrix,color='white',fontsize=12,
            bbox={'facecolor':'red', 'alpha':0.5, 'pad':2.6})

    #---------------------------------------------
#    ax.set_title("Nestedness,"+Dataset)
    # ax.set_xlabel("Groups of Microbes",fontsize=12)
    # ax.set_ylabel("Groups of Patients",fontsize=12)
    ax.set_xlabel("Groups of Microbes",fontsize=12)
    ax.set_ylabel("Enterotypes",fontsize=12)
    #---------------------------------------------
    plt.savefig(\
    path_plot + 'Mutualistic_P_Matrix_'+Dataset+'.png',dpi=300)
    plt.show()
     #================================================================
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

