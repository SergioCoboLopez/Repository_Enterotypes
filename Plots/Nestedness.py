#8/02/2019. Con este codigo calculamos la nestedness de acuerdo a la formula del paper de
#Saavedra_et_al-Ecology_and_Evolution. Ademas, planteamos varios null models. En particular un Erdos-Renyi y un
#link swap en su version canonica/probabilistica. Como output, te devuelve una distribucion de las randomizaciones
#de la nestedness con la nestedness "de verdad".

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
    conteiro=0
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
            print shared_interactions_patients, conteiro
            conteiro+=1
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

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#1. PARAMETERS FOR INPUT/OUTPUT DATA
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#1.1. EXTERNAL PARAMETERS
#=====================
Dataset=sys.argv[1]#==
number_of_links=3  #==
rnd.seed(1)        #==
Diseases=0         #==
#=====================

#1.2. INTERNAL PARAMETERS
#===========================================================================================================

#::::::::::::::::::::
K=10;L=20

Dataset_Pieces=Dataset.split("_")
if len(Dataset_Pieces)>3:
    K=int(Dataset_Pieces[2][1:])
    L=int(Dataset_Pieces[3][1:])
    print K, L
    

if "L13" in Dataset:
    L=13
#::::::::::::::::::::
    
Patients_per_Dataset=\
{'S-8_Tot':107, 'V-10_Tot':92, 'V-22_Tot':467, 'V-23_24_Tot':222,'V-25_Tot':883,\
 'He_CD':47, 'He_Healthy':39, 'Nielsen_IBD':78, 'Nielsen_Healthy':58, \
 'Vogtmann_CRC': 47 , 'Vogtmann_Healthy':30 , 'Zeller_CRC': 135 , 'Zeller_Healthy':42 ,\
 'Zhang_RA': 92,'Zhang_Healthy':55,\
 'He_CD_L13':47, 'He_Healthy_L13':39, 'Nielsen_IBD_L13':78, 'Nielsen_Healthy_L13':58, \
 'Vogtmann_CRC_L13': 47 , 'Vogtmann_Healthy_L13':30 , 'Zeller_CRC_L13': 135 , 'Zeller_Healthy_L13':42 ,\
 'Zhang_RA_L13': 92,'Zhang_Healthy_L13':55,\
 'Feng_CRC':85 ,'Feng_CRC_K10_L13':85, 'Feng_CRC_K5_L13':85, 'Feng_CRC_K5_L20':85,  \
 'Feng_Healthy':21 ,'Feng_Healthy_K10_L13':21 , 'Feng_Healthy_K5_L13':21 , 'Feng_Healthy_K5_L20':21,  \
 'Jie_CVD':152 ,'Jie_CVD_K10_L13':152 , 'Jie_CVD_K5_L13':152 , 'Jie_CVD_K5_L20':152,  \
 'Jie_Healthy': 75 ,'Jie_Healthy_K10_L13': 75 , 'Jie_Healthy_K5_L13':75 , 'Jie_Healthy_K5_L20':75,  \
 'Schirmer_IBD': 68 ,'Schirmer_IBD_K10_L13': 68 , 'Schirmer_IBD_K5_L13': 68 , 'Schirmer_IBD_K5_L20': 68,  \
 'Schirmer_Healthy':17,'Schirmer_Healthy_K10_L13':17,'Schirmer_Healthy_K5_L13':17,'Schirmer_Healthy_K5_L20': 17}


Microbes_per_Dataset=\
{'S-8_Tot':128, 'V-10_Tot':137, 'V-22_Tot':134, 'V-23_24_Tot':118,'V-25_Tot':144,\
 'He_CD':313, 'He_Healthy':313, 'Nielsen_IBD':313 , 'Nielsen_Healthy':313,\
 'Vogtmann_CRC':313 , 'Vogtmann_Healthy':313 , 'Zeller_CRC':313 , 'Zeller_Healthy':313 ,\
 'Zhang_RA':313,'Zhang_Healthy':313,\
 'He_CD_L13':313, 'He_Healthy_L13':313, 'Nielsen_IBD_L13':313 , 'Nielsen_Healthy_L13':313,\
 'Vogtmann_CRC_L13':313 , 'Vogtmann_Healthy_L13':313 , 'Zeller_CRC_L13':313 , 'Zeller_Healthy_L13':313 ,\
 'Zhang_RA_L13':313,'Zhang_Healthy_L13':313, \
 'Feng_CRC':296 ,'Feng_CRC_K10_L13':296 ,'Feng_CRC_K5_L13':296 ,'Feng_CRC_K5_L20':296,  \
 'Feng_Healthy':215 ,'Feng_Healthy_K10_L13':215 , 'Feng_Healthy_K5_L13':215 , 'Feng_Healthy_K5_L20':215,  \
 'Jie_CVD':304 ,'Jie_CVD_K10_L13':304 , 'Jie_CVD_K5_L13':304 , 'Jie_CVD_K5_L20':304,  \
 'Jie_Healthy':270 ,'Jie_Healthy_K10_L13':270 , 'Jie_Healthy_K5_L13':270 , 'Jie_Healthy_K5_L20':270,  \
 'Schirmer_IBD':206 ,'Schirmer_IBD_K10_L13':206 , 'Schirmer_IBD_K5_L13':206 , 'Schirmer_IBD_K5_L20':206,  \
 'Schirmer_Healthy':145,'Schirmer_Healthy_K10_L13':145 ,'Schirmer_Healthy_K5_L13':145, \
 'Schirmer_Healthy_K5_L20':145}

Patients=Patients_per_Dataset[Dataset]
Microbes=Microbes_per_Dataset[Dataset]
#===========================================================================================================

if Diseases==1:
    #1.3. PATHS FOR FILES
    #=================================================================================
    Path_out='../../Output_Data/'+ 'Datasets_Diseases/' + Dataset + '/'
    Path_in= '../../Input_Data/' + 'Datasets_Diseases/'
    path_plot='/export/home/shared/Projects/Microbiome/Plots/Second_Path/Nestedness/'
    #=================================================================================

    #1.4. CHOOSE BEST LIKELIHOOD
    #=====================================================
    Likelihood_Seed={}
    for Seed in range(1,50):
        likelihood_file=\
        Dataset + '_Seed_' + str(Seed) + '_LogLikelihood.txt'
        likelihoodfinal=read_file(Path_out + likelihood_file)[-1]
        Likelihood_Seed[float(likelihoodfinal[0])]=Seed

    print Likelihood_Seed
    print max(Likelihood_Seed)
    print Likelihood_Seed[max(Likelihood_Seed)]
    Seed=Likelihood_Seed[max(Likelihood_Seed)]
    
    filename=Dataset + '_Seed_' + str(Seed) + '_Parameters.txt'
    Data=read_file(Path_out+filename)
    #=====================================================
    
else:
    #1.3. PATHS FOR FILES
    #==============================================================
    Path_out='../../Output_Data/'+ 'Leave_One_Out_' + Dataset + '/'
    Path_in= '../../Input_Data/' + 'Leave_One_Out_' + Dataset + '/'
    path_plot='/export/home/shared/Projects/Microbiome/Plots/First_Path/Nestedness/'
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
print pMatrix_ordered
print "Nestedness",Nestedness_pMatrix, "\n"
#=======================================================================

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#2. NULL MODELS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rows=np.shape(pMatrix_ordered)[0]
cols=np.shape(pMatrix_ordered)[1]

# #2.1. Non-probabilistic version (microcanonical Erdos-Renyi)
# #================================================================================

# #2.1.1. Reshape matrix to vector
# #::::::::::::::::::::::::::::::::::::
# pMatrix_1D=np.reshape(pMatrix,rows*cols)
# #::::::::::::::::::::::::::::::::::::

# #2.1.2.Copy vector and shuffle multiple times
# #:::::::::::::::::::::::::::::::::::::::::::::::::::::

# np.random.seed(seed=1)
# Nestedness_Vector_E_R=[]
# for iteration_0 in range(1000):
    
#     #Copy and shuffle
#     #.......................................................
#     pMatrix_1D_shuffled=cp.copy(pMatrix_1D)
#     np.random.shuffle(pMatrix_1D_shuffled)
#     pMatrix_E_R=np.reshape(pMatrix_1D_shuffled,(rows,cols))
#     #.......................................................

#     #Get and save nestedness
#     #......................................................................
#     pMatrix_E_R_ordered,Nestedness_pMatrix_E_R,Dict_row_membership_E_R,Dict_col_membership_E_R=\
#     Nestedness_fun(pMatrix_E_R)
#     Nestedness_Vector_E_R.append(Nestedness_pMatrix_E_R)
#     #......................................................................
    
# #Significance of randomized nestedness
# #------------------------------------------
# counter_E_R=0
# for element_0 in Nestedness_Vector_E_R:
#     if element_0 >= Nestedness_pMatrix:
#         counter_E_R+=1

# print "Nestedness Soft Erdos-Renyi"        
# print "pvalue", counter_E_R/float(len(Nestedness_Vector_E_R))
# #------------------------------------------
# #:::::::::::::::::::::::::::::::::::::::::::::::::::::

# #================================================================================

# #2.2. Probabilistic version (Canonical Erdos-Renyi)
# #================================================================================

# #2.2.1. Generate matrix of probabilities
# #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Total_Degree=sum(sum(pMatrix))
# Probability=Total_Degree/(rows*cols)
# Probability_Matrix=np.full((rows,cols), Probability)
# #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# #2.2.2. Generate many randomized networks 
# #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# #Generate many randomized networks and save nestedness
# #-------------------------------------------------------------------------------
# Nestedness_Vector_E_R_Canonical=[]
# for i in range(1000):
#     E_R_Canonical_ordered,Nestedness_E_R_Canonical=Null_Model_fun(Probability_Matrix)
#     Nestedness_Vector_E_R_Canonical.append(Nestedness_E_R_Canonical)
# #-------------------------------------------------------------------------------

# #Significance of randomized nestedness
# #------------------------------------------
# counter_E_R_Canonical=0
# for element_1 in Nestedness_Vector_E_R_Canonical:
#     if element_1 >= Nestedness_pMatrix:
#         counter_E_R_Canonical+=1

# print "Nestedness Soft Erdos-Renyi"        
# print "pvalue", counter_E_R_Canonical/float(len(Nestedness_Vector_E_R_Canonical)), '\n'
# #------------------------------------------
                                           
# #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# #================================================================================

# #2.3. Degree-Distribution-Probabilistic
# #================================================================================

# #2.3.1. Generate matrix of probabilities
# #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Degree_Probability_Matrix=np.zeros((rows,cols))

# for row in range(rows):
#     row_prob=sum(pMatrix_ordered[row])/cols
#     for column in range(cols):
#         col_prob=sum(np.transpose(pMatrix_ordered)[column])/rows
#         Average_prob=0.5*(row_prob + col_prob)
#         Degree_Probability_Matrix[row,column]=Average_prob

# print "Degree probability matrix"
# print Degree_Probability_Matrix  , '\n'      
# #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# #2.3.2. Generate randomized networks
# #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# #Generate many randomized networks and save nestedness
# #-------------------------------------------------------------------------------
# Nestedness_Vector=[]
# for i in range(1000):
#     Degree_Null_Model,Nestedness_Degree_Null_Model=Null_Model_fun(Degree_Probability_Matrix)
#     Nestedness_Vector.append(Nestedness_Degree_Null_Model)
# #-------------------------------------------------------------------------------

# #Significance of randomized nestedness
# #------------------------------------------
# counter=0
# for element_2 in Nestedness_Vector:
#     if element_2 >= Nestedness_pMatrix:
#         counter+=1

# print "Soft degree distribution randomization (Canonical Link Swap)"
# print "pvalue", counter/float(len(Nestedness_Vector))
# print "maximum nestedness randomization: ", max(Nestedness_Vector)
# print "minimum nestedness randomization: ", min(Nestedness_Vector)
# #------------------------------------------

# #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# #================================================================================

# #2.4. Node swapping
# #================================================================================

# for iteration in range(200):
#     #Propose Tuples for swapping
#     #========================================================================
#     Removed_Tuple_0, Removed_Tuple_1=rnd.sample(Links_Nodes_Dict[1],2)

#     Tentative_Tuple_0=( Removed_Tuple_0[0],Removed_Tuple_1[1] )
#     Tentative_Tuple_1=( Removed_Tuple_1[0],Removed_Tuple_0[1] )

#     counter_while=0
#     while Nodes_Links_Dict[Tentative_Tuple_0]==1 \
#           or Nodes_Links_Dict[Tentative_Tuple_1]==1:
#         Removed_Tuple_0,Removed_Tuple_1=rnd.sample(Links_Nodes_Dict[1],2)
#         Tentative_Tuple_0=( Removed_Tuple_0[0],Removed_Tuple_1[1] )
#         Tentative_Tuple_1=( Removed_Tuple_1[0],Removed_Tuple_0[1] )
#         counter_while+=1
        
#         if counter_while>100:
#             print "No encontramos ninguna tupla compatible"
#             break

#     if counter_while>100:
#         break
    
#     #========================================================================

#     #Do swapping
#     #========================================================================
#     New_Tuple_0=Tentative_Tuple_0;New_Tuple_1=Tentative_Tuple_1

#     #Removed Tuples: update dictionaries
#     #::::::::::::::::::::::::::::::::::::::::::::
#     Links_Nodes_Dict[0].append(Removed_Tuple_0)
#     Links_Nodes_Dict[0].append(Removed_Tuple_1)

#     Links_Nodes_Dict[1].remove(Removed_Tuple_0)
#     Links_Nodes_Dict[1].remove(Removed_Tuple_1)

#     Nodes_Links_Dict[Removed_Tuple_0]=0;Nodes_Links_Dict[Removed_Tuple_1]=0
#     #::::::::::::::::::::::::::::::::::::::::::::

#     #New Tuples: Update dictionaries
#     #::::::::::::::::::::::::::::::::::::::::::::
#     Links_Nodes_Dict[0].remove(New_Tuple_0)
#     Links_Nodes_Dict[0].remove(New_Tuple_1)

#     Links_Nodes_Dict[1].append(New_Tuple_0)
#     Links_Nodes_Dict[1].append(New_Tuple_1)

#     Nodes_Links_Dict[New_Tuple_0]=1;Nodes_Links_Dict[New_Tuple_1]=1
#     #::::::::::::::::::::::::::::::::::::::::::::
    

# pMatrix_Randomized=np.zeros((K,L))

# for row in range(K):
#     for col in range(L):
#         Tuple=(row,col)
#         pMatrix_Randomized[row][col]=Nodes_Links_Dict[Tuple]

# pMatrix_Rand_ordered,Nestedness_2,Dict_row_membership_Rand,Dict_col_membership_Rand=\
# Nestedness_fun(Degree_Probability_Matrix)

#========================================================================
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#3. PLOTS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Information=["Only_Nestedness", "Nestedness_w_Groups", "Nestedness_w_Randomization", "Randomization_Nestedness"]
Plot=Information[-1]

if Plot=="Only_Nestedness":
    
    #Only Nestedness
    #=====================================================================================================
    ax=sns.heatmap(pMatrix_ordered,cmap="RdBu_r",linewidths=1,vmin=0,vmax=1,\
                   xticklabels=False,yticklabels=False, cbar_kws={'label': 'P(Connection)'})
    ax.figure.axes[-1].yaxis.label.set_size(15)

    ax.text(0.75, 0.75, 'nestedness= %f' %Nestedness_pMatrix,color='white',fontsize=12,
            bbox={'facecolor':'red', 'alpha':0.5, 'pad':2.6})

    #---------------------------------------------
#    ax.set_title("Nestedness,"+Dataset)
    ax.set_xlabel("Communities of Microbes",fontsize=15)
    ax.set_ylabel("Communities of Patients",fontsize=15)
    #---------------------------------------------
    plt.savefig(\
    path_plot + 'Mutualistic_P_Matrix_'+Dataset+'.png',dpi=300)
    plt.show()
    #=====================================================================================================

elif Plot=="Nestedness_w_Groups":
    #=====================================================================================================

    #3.2.0. Shapes and colors for groups of microbes/patients
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    colors_shapes_microbes={}
    shapes=["^", "v"]
    counter=0

    for shape in range(len(shapes)):
        for group in range(L/2):
            colors_shapes_microbes[counter]=( group, shapes[shape] )
            counter+=1
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
    #3.2.1. Gridspec Parameters
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    fig=plt.figure(figsize=(14,10))
    gs = gridspec.GridSpec(2, 2, height_ratios=[0.11,1],width_ratios=[0.2,1])
    gs.update(left=0.05, right=0.999,top=0.95,bottom=0.1, wspace=0.05,hspace=0.1)
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #3.2.2. Plot Triangles
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    triangles_up = plt.subplot(gs[0, 1])
    sns.despine(offset=0, ax=triangles_up, top=True, left=True,bottom=True)#Poner fondo de la grafica en blanco

    #3.2.2.1. Triangle parameters and colormap
    #--------------------------------------------------
    Base=0.033; Height=0.42
    Buffer=0.0035
    Right_Vertex=[-Buffer, 0]
    colors_microbes = cm.Paired(np.linspace(0, 1, K))
    #--------------------------------------------------

    #3.2.2.2. Loop to draw triangles
    #---------------------------------------------------------------------------------------------------
    for triangle in range(L):

        #Get membership corresponding to column
        #....................................................
        Membership_of_column=Dict_column_membership[triangle]
        #....................................................

        #Set up,down triangles according to membership: shape, color dictionary
        #.....................................................................................
        if colors_shapes_microbes[Membership_of_column][1]=="^":
            Origin=[Right_Vertex[0]+2*Buffer,0.0*Height];Right_Vertex=[Origin[0]+Base, Origin[1]]
            Top_Vertex=[Right_Vertex[0]-0.5*Base,Origin[1]+ Height]

        else:
            Origin=[Right_Vertex[0]+2*Buffer,Height];Right_Vertex=[Origin[0]+Base, Origin[1]]
            Top_Vertex=[Right_Vertex[0]-0.5*Base,Origin[1]- Height]
        #.....................................................................................

        #Draw and color triangle
        #.............................................................................
        Vertices=[ Origin, Right_Vertex, Top_Vertex, Origin]
        triangles_up.add_patch(patches.Polygon( Vertices,\
        facecolor=colors_microbes[ colors_shapes_microbes[Membership_of_column][0]] ))
        #.............................................................................

        #---------------------------------------------------------------------------------------------------

    #3.2.2.3. Put title and remove ticks and axis
    #-----------------------------------------------------------------
    plt.title(Dataset+"            ",fontsize=18)

    #Quitar xticks,yticks y ejes
    triangles_up.tick_params(
        axis='both',         # changes apply to the x-axis
        which='both',        # both major and minor ticks are affected
        bottom='off',        # ticks along the bottom edge are off
        top='off',           # ticks along the top edge are off
        labelbottom='off',
        left='off',
        labelleft='off'     )
    #-----------------------------------------------------------------

    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #3.2.3. PLOT SQUARES
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    squares_left=plt.subplot(gs[1, 0])
    sns.despine(offset=0, ax=squares_left, top=True, left=True,bottom=True)#Poner fondo de la grafica en blanco

    #3.2.3.1. Square parameters and colormap
    #----------------------------------------------
    Origin_x_left=0.55;
    Width_left=0.18;Height_left=0.06;
    Separation=0.02
    Origin_y_left=1-Separation
    colors_patients = cm.Set3(np.linspace(0, 1, K))
    #---------------------------------------------- 

    #3.1.3.2. Loop to draw squares
    #-----------------------------------------------------------------
    for square in range(K):
        print square, Dict_row_membership[square]

        squares_left.add_patch(patches.Rectangle(\
        (Origin_x_left,Origin_y_left),Width_left,-Height_left,\
        facecolor=colors_patients[ Dict_row_membership[square] ]))

        Origin_y_left=Origin_y_left - Height_left - 2*Separation
    #-----------------------------------------------------------------


    #3.2.3.3. Put title and remove ticks and axis
    #-----------------------------------------------------------------
    #Quitar xticks,yticks y ejes
    squares_left.tick_params(
        axis='both',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom='off',      # ticks along the bottom edge are off
        top='off',         # ticks along the top edge are off
        labelbottom='off',
        left='off',
        labelleft='off'     )
    #-----------------------------------------------------------------

    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


    #3.2.3. Plot Mutualistic Network
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ax = plt.subplot(gs[1, 1])
    ax=sns.heatmap(pMatrix_ordered,cmap="RdBu_r",linewidths=1,vmin=0,vmax=1)

    ax.text(0.75, 0.75, 'n= %f' %Nestedness_pMatrix,color='white',fontsize=20,
            bbox={'facecolor':'red', 'alpha':0.5, 'pad':2.6})

    #yticks, xticks, ylabels, xlabels
    #---------------------------------------------
    yticks_pos=[row + 0.5 for row in range(K)]
    yticks_labels=Dict_row_membership.values()
    ax.set_yticks(yticks_pos)
    ax.set_yticklabels(yticks_labels,fontsize=15)

    xticks_pos=[col + 0.5 for col in range(L)]
    xticks_labels=Dict_column_membership.values()
    ax.set_xticks(xticks_pos)
    ax.set_xticklabels(xticks_labels,fontsize=15)

    ax.set_xlabel("Groups of Microbes",fontsize=18)
    ax.set_ylabel("Groups of Patients",fontsize=18)
    #---------------------------------------------

    plt.savefig(\
    path_plot + 'Mutualistic_P_Matrix_w_Memberships_'+Dataset+'.pdf',dpi=300)

    plt.savefig(\
    path_plot + 'Mutualistic_P_Matrix_w_Memberships_'+Dataset+'.png',dpi=300)

    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    plt.show()
    #============================================================================================================
    
elif Plot=="Nestedness_w_Randomization":
    #======================================================================================================

    #3.3.1. Gridspec Parameters
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    fig=plt.figure(figsize=(12,5))
    gs = gridspec.GridSpec(1, 2)
    gs.update(left=0.07, right=0.95,top=0.95,bottom=0.1, wspace=0.15, hspace=0.2)
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #3.3.2. Plot Nestedness
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ax_00 = plt.subplot(gs[0, 0])
    ax_00=sns.heatmap(pMatrix_ordered,cmap="RdBu_r",linewidths=1,vmin=0,vmax=1)

    ax_00.text(0.75, 0.75, 'n= %f' %Nestedness_pMatrix,color='white',fontsize=20,
            bbox={'facecolor':'red', 'alpha':0.5, 'pad':2.6})

    #yticks, xticks, ylabels, xlabels
    #---------------------------------------------
    yticks_pos=[row + 0.5 for row in range(K)]
    yticks_labels=Dict_row_membership.values()
    ax_00.set_yticks(yticks_pos)
    ax_00.set_yticklabels(yticks_labels)
    
    xticks_pos=[col + 0.5 for col in range(L)]
    xticks_labels=Dict_column_membership.values()
    ax_00.set_xticks(xticks_pos)
    ax_00.set_xticklabels(xticks_labels)
    
    ax_00.set_xlabel("Groups of Microbes")
    ax_00.set_ylabel("Groups of Patients")
    #---------------------------------------------

    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #3.3.3. Plot Randomization
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ax_01 = plt.subplot(gs[0,1])
    ax_01=sns.distplot(Nestedness_Vector)
    ax_01.axvline(x=Nestedness_pMatrix, linewidth=4, color='k')
    ax_01.set_title(Dataset)
    ax_01.set_xlabel('Nestedness')
    ax_01.set_ylabel('Density')
    ax_01.text(min(Nestedness_Vector), 8,"pvalue = " + str(counter/float(len(Nestedness_Vector))), \
             bbox={'facecolor':'red', 'alpha':0.5, 'pad':2.6} )
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #Save file
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    plt.show()
    plt.savefig(\
                path_plot + 'Mutualistic_P_Matrix_w_Randomization_'+Dataset+'.png',dpi=300)
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
            
    #======================================================================================================
    
    
elif Plot=="Randomization_Nestedness":
    #============================================================================
    sns.distplot(Nestedness_Vector)
    plt.axvline(x=Nestedness_pMatrix, linewidth=4, color='k')
    plt.title(Dataset)
    plt.xlabel('Nestedness')
    plt.ylabel('Density')
    plt.text(min(Nestedness_Vector), 8,"pvalue = " + str(counter/float(len(Nestedness_Vector))), \
             bbox={'facecolor':'red', 'alpha':0.5, 'pad':2.6} )
    plt.savefig(path_plot + Dataset + '_Significance_Nestedness.pdf',dpi=300)
    plt.savefig(path_plot + Dataset + '_Significance_Nestedness.png',dpi=300)
    plt.show()
    #============================================================================
    
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
