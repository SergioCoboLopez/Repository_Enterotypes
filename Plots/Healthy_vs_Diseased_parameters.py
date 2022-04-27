#17/07/2019. Entrenamos el MMSBM con los datasets healthy y diseased juntos y miramos a ver que pasa con los
#parametros del modelo

import numpy as np
import math
import operator
import matplotlib.gridspec as gridspec
import matplotlib
from pylab import *
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import seaborn as sns


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

#0.2. Shannon Entropy
#===========================================                                                                      
def Shannon(v,base):
    Entropy=0
    for i in v:
        if i==0:
            Entropy+=0
        else:
            Entropy+=-i*math.log(i,base)

    return Entropy
#===========================================

#0.3. Nestedness 
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
    Nestedness=float(d_ij_microbes+d_ij_patients)\
    /float(Min_Interactions_microbes+Min_Interactions_patients)
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::         

    return Matrix_ordered, Nestedness, Dict_row_membership_fun, Dict_column_membership_fun
#================================================================================================

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++


#1. MAIN
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#1.1. Input Parameters
#===================================================================================
Path_Input='/export/home/shared/Projects/Microbiome/Output_Data/Datasets_Diseases/'
Path_Plot='/export/home/shared/Projects/Microbiome/Plots/Second_Path/Aggregated_Datasets/'

Datasets=['All_He', 'All_Nielsen', 'All_Vogtmann', 'All_Zeller', 'All_Zhang']

K=10;L=20

Patients_per_Dataset=\
{'All_He':{'Healthy':39, 'Diseased':47}, 'All_Nielsen':{'Healthy':58, 'Diseased':78 } , \
'All_Vogtmann': {'Healthy':30, 'Diseased':49 } , 'All_Zeller': {'Healthy': 42, 'Diseased':136 }, \
'All_Zhang': {'Healthy': 55, 'Diseased':92 } }

Microbes_per_Dataset=\
{'All_He':263, 'All_Nielsen':280 , 'All_Vogtmann':247 , 'All_Zeller':291 ,\
 'All_Zhang': 298}
#===================================================================================

#1.2. Compute means and entropies of membership vectors for all datasets
#=============================================================================================================
Healthy_Metrics={}; Diseased_Metrics={};Nested_Metrics={}

for Data_read in Datasets:
    
    #3.1. Choose best likelihood                                                                                  
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::        
    path_Out='/export/home/shared/Projects/Microbiome/Plots/Second_Path/Aggregated_Datasets/'

    Likelihood_Seed={}
    for Seed in range(1,10,2):
        Likelihood_File=\
        Data_read + '/' + Data_read + '_Seed_' + str(Seed) + '_LogLikelihood.txt'
        LikelihoodFinal=read_file(Path_Input + Likelihood_File)[-1]
        Likelihood_Seed[float(LikelihoodFinal[0])]=Seed
        Seed=Likelihood_Seed[max(Likelihood_Seed)]
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #3.2. Parameters of Data                                           
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    Patients=Patients_per_Dataset[Data_read]['Healthy'] + Patients_per_Dataset[Data_read]['Diseased'] 
    Microbes=Microbes_per_Dataset[Data_read]

    File_Data=Data_read + '/' + Data_read + '_Seed_' + str(Seed) + '_Parameters.txt'
    Data=read_file(Path_Input+File_Data)
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


    #3.3. Read and save membership vectors                                                                       
    #:::::::::::::::::::::::::::::::::::::::::::::
    theta=[];eta=[]
    theta=Data[:Patients]
    eta=Data[Patients+1:Patients+1+Microbes]
    #:::::::::::::::::::::::::::::::::::::::::::::

    #3.4. Format microbe and patient membership vectors to float
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    Metrics_Microbe={}
    for microbe in range(len(eta)):
        eta[microbe] = [float(membership) for membership in eta[microbe] ]
        Metrics_Microbe[microbe]={'Membership': eta[microbe], 'Shannon': Shannon(eta[microbe], L) }

    Metrics_Patient={}
    for patient in range(len(theta)):
        theta[patient] = [float(membership) for membership in theta[patient] ]
        Metrics_Patient[patient]={'Membership': theta[patient], 'Shannon': Shannon(theta[patient], K) }
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #3.5. Split healthy and diseased patients
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    Metrics_Patient_Healthy = dict(Metrics_Patient.items()[:Patients_per_Dataset[Data_read]['Healthy']])
    Metrics_Patient_Diseased = dict(Metrics_Patient.items()\
    [Patients_per_Dataset[Data_read]['Healthy']:\
     Patients_per_Dataset[Data_read]['Healthy']+1 + Patients_per_Dataset[Data_read]['Diseased']])

    Healthy_Metrics[Data_read]=Metrics_Patient_Healthy.values()
    Diseased_Metrics[Data_read]=Metrics_Patient_Diseased.values()
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


    #3.6. NESTEDNESS
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
    #3.6.1. Read Pmatrices and format into matrices
    #-------------------------------------------------------------------
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

    Nested_Matrix,Nestedness,Dict_row_membership,Dict_column_membership=Nestedness_fun(pMatrix)
    print Nested_Matrix
    print "Nestedness",Nestedness, "\n"

    Nested_Metrics[Data_read]={'Matrix':Nested_Matrix, 'Nestedness_value':Nestedness}

    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#=============================================================================================================

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#2. PLOTS. GRIDSPEC 5X3: 5Datasets; means, shannons and nestedness of joint datasets
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#2.1. Single Panels
#==============================================================================
Data_to_Plot='All_Zhang' #'All_He' 'All_Nielsen', 'All_Vogtmann', 'All_Zeller', 'All_Zhang'

#2.1.1. Gridspec
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
fig=plt.figure(figsize=(19,5))
gs = gridspec.GridSpec(1, 3)
gs.update(left=0.05, right=0.999,top=0.95,bottom=0.1, wspace=0.2,hspace=0.6)
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


#2.1.2. Compute Aggregated Memberships
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#2.1.2.1. Healthy Aggregated Memberships
#------------------------------------------------------------------------------------------------------------

#Extract Memberships
Healthy_Memberships=[element_M['Membership'] for element_M in  Healthy_Metrics[Data_to_Plot]]

#Sum
Aggregated_Healthy_Membership_raw=[0]*K
for Memberships in Healthy_Memberships:
    for k in range(K):
        Aggregated_Healthy_Membership_raw[k]+=Memberships[k]

#Normalize
Aggregated_Healthy_Membership=[0]*K
for k_norm_h in range(len(Aggregated_Healthy_Membership)):
    Aggregated_Healthy_Membership[k_norm_h]=\
    Aggregated_Healthy_Membership_raw[k_norm_h]/sum(Aggregated_Healthy_Membership_raw)

print Aggregated_Healthy_Membership
#------------------------------------------------------------------------------------------------------------

#2.1.2.2. Diseased Aggregated Memberships
#------------------------------------------------------------------------------------------------------------

#Extract Memberships
Diseased_Memberships=[element_M['Membership'] for element_M in  Diseased_Metrics[Data_to_Plot]]

#Sum
Aggregated_Diseased_Membership_raw=[0]*K
for Memberships in Diseased_Memberships:
    for k in range(K):
        Aggregated_Diseased_Membership_raw[k]+=Memberships[k]

#Normalize
Aggregated_Diseased_Membership=[0]*K
for k_norm_d in range(len(Aggregated_Diseased_Membership)):
    Aggregated_Diseased_Membership[k_norm_d]=\
    Aggregated_Diseased_Membership_raw[k_norm_d]/sum(Aggregated_Diseased_Membership_raw)

print Aggregated_Diseased_Membership, sum(Aggregated_Diseased_Membership)
#------------------------------------------------------------------------------------------------------------

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#2.1.3. Compute Shannon Entropies
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Healthy_Shannons=[element_S['Shannon'] for element_S in Healthy_Metrics[Data_to_Plot]]
Diseased_Shannons=[element_S['Shannon'] for element_S in  Diseased_Metrics[Data_to_Plot]]
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#2.1.4. Plot Figures
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#Memberships
#-----------------------------------------------------------------------------------------------
Memberships=plt.subplot(gs[0, 0])

#Barplots
ind = np.arange(K)
width=0.4
Healthy=plt.bar(ind - 0.5*width, Aggregated_Healthy_Membership, width )
Diseased=plt.bar(ind + 0.5*width, Aggregated_Diseased_Membership, width)

#Axes and ticks
plt.title("Membership Distribution of Healthy vs Diseased Patients: " + Data_to_Plot)
plt.legend((Healthy[0], Diseased[0]), ('Healthy', 'Diseased'))
plt.xlabel("Membership")
plt.ylabel("Membership ratio")
Memberships.set_xticks(range(K))
#------------------------------------------------------------------------------------


#Shannon Entropies
#----------------------------------------------------------------------------------------------------
Shannons=plt.subplot(gs[0, 1])

#Distplots
sns.distplot(Healthy_Shannons, kde=False, norm_hist=True, bins=10, label='Healthy')
sns.distplot(Diseased_Shannons, kde=False, norm_hist=True, bins=10, label='Diseased')


#Axes and ticks
plt.title("Shannon Entropy Healthy vs Diseased Patients: " + Data_to_Plot)
plt.xlabel("H")
plt.ylabel("Density")
plt.legend()
Memberships.set_xticks(range(K))
#----------------------------------------------------------------------------------------------------


#Nestedness
#----------------------------------------------------------------------------------------------------
Nestedness_Plot=plt.subplot(gs[0,2])

Nestedness_Plot=sns.heatmap(Nested_Metrics[Data_to_Plot]['Matrix'],cmap="RdBu_r",linewidths=1,vmin=0,vmax=1)
plt.title("Nestedness")

#yticks, xticks, ylabels, xlabels                                                                              
#---------------------------------------------                                                                 
yticks_pos=[row + 0.5 for row in range(K)]
yticks_labels=Dict_row_membership.values()
Nestedness_Plot.set_yticks(yticks_pos)
Nestedness_Plot.set_yticklabels(yticks_labels,fontsize=10)

xticks_pos=[col + 0.5 for col in range(L)]
xticks_labels=Dict_column_membership.values()
Nestedness_Plot.set_xticks(xticks_pos)
Nestedness_Plot.set_xticklabels(xticks_labels,fontsize=10)
#---------------------------------------------

plt.xlabel("Groups of Microbes")
plt.ylabel("Groups of Patients")

Nestedness_Plot.text(0.75, 0.75, 'n= %f' %Nested_Metrics[Data_to_Plot]['Nestedness_value'],color='white',fontsize=20,
            bbox={'facecolor':'red', 'alpha':0.5, 'pad':2.6})
#----------------------------------------------------------------------------------------------------
plt.savefig(\
    Path_Plot + 'Whole_Datasets_'+Data_to_Plot+'.png',dpi=300)


plt.show()
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#==============================================================================


'''
#2.2. Full Panel
#==============================================================================


#2.1.1. Gridspec
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
fig=plt.figure(figsize=(15,20))
gs = gridspec.GridSpec(5, 3)
gs.update(left=0.05, right=0.999,top=0.95,bottom=0.1, wspace=0.5,hspace=0.6)
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#2.1.2. Plot Figures
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
counter_plot=0
for Data_Plot in Datasets:

    Healthy_Means=[element_M['Mean'] for element_M in  Healthy_Metrics[Data_Plot]]
    Healthy_Shannons=[element_S['Shannon'] for element_S in Healthy_Metrics[Data_Plot]]

    Diseased_Means=[element_M['Mean'] for element_M in  Diseased_Metrics[Data_Plot]]
    Diseased_Shannons=[element_S['Shannon'] for element_S in  Diseased_Metrics[Data_Plot]]

    Means=plt.subplot(gs[counter_plot, 0])
    sns.distplot(Healthy_Means, norm_hist=True)
    sns.distplot(Diseased_Means, norm_hist=True)


    Shannons=plt.subplot(gs[counter_plot, 1])
    sns.distplot(Healthy_Shannons, norm_hist=True)
    sns.distplot(Diseased_Shannons, norm_hist=True)


    Nestedness=plt.subplot(gs[counter_plot,2])
    Nestedness=sns.heatmap(Nested_Metrics[Data_read]['Matrix'],cmap="RdBu_r",linewidths=1,vmin=0,vmax=1)

    counter_plot+=1
    
plt.show()
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#==============================================================================    





#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

'''
