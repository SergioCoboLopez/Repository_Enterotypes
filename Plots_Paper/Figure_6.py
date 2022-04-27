#20/04/2020. Calculamos el rango de microbios y de hosts para cada
#dataset y en total.
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
    sorted_sum_cols = \
    sorted(sum_cols.items(), key=operator.itemgetter(1),reverse=True)
    #-------------------------------------------

    #Reorder matrix rows by the sum of the columns
    #-------------------------------------------
    Dict_row_membership_fun={}
    Matrix_row_ordered=np.zeros((rows,columns))
    for new_row in range(rows):
        old_row=sorted_sum_cols[new_row][0]
        
        Matrix_row_ordered[new_row]=\
        Matrix[old_row]
        
        Dict_row_membership_fun[new_row]=old_row
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

        old_col=sorted_sum_rows[new_col][0]
        
        np.transpose(Matrix_ordered)[new_col] = \
        np.transpose(Matrix_row_ordered)[old_col]
        Dict_column_membership_fun[new_col]=old_col
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
                min(Interactions_microbes_i,Interactions_microbes_j)
            #.......................................................

            #-------------------------------------------------------

    Nestedness=float(d_ij_microbes+d_ij_patients)/float(Min_Interactions_microbes+Min_Interactions_patients)
    
    return Matrix_ordered, Nestedness, Dict_row_membership_fun, Dict_column_membership_fun
#===================================================================

#0.3. Compute Nestedness
#===================================================================
def Genomic_Size_fun(microbial_species,Dataset,Missing_Microbes):

    #Cogemos informacion del dataframe                            
    #------------------------------------------------------------ 
    Matching_Microbes=Genomic_Data_Raw.loc\
    [Genomic_Data_Raw['#Organism Name']==microbe]
    
    Sizes_Raw=Matching_Microbes['Size(Mb)_y']
    #------------------------------------------------------------

    #Si no encuentra la especie miramos a ver si se ha equivocado 
    #............................................................
    if Matching_Microbes.empty==True:

        #We start looking at Vinod's File                         
        #------------------------------------------------------   
        Alter_Info=Missing_Microbes\
        [Missing_Microbes['Total']==microbe]               
        
        #If we find something on Vinod's file we enter there
        #-------------------------------------------------------  
        if Alter_Info.empty==False:                              
            Alter_Info=Alter_Info.reset_index(drop=True)         
            Alter_Name=Alter_Info['Alternative species names'][0]
            In_Dataset=Alter_Info['Found in dataset'][0]         
                                                                     
            #If the microbe exists in our dataset with other name
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
                    if ("_"+microbe_splits[3]+"_" in Name_Genomic\
                    and microbe_splits[3]!="bacterium"):
                        Matching_Microbes=Genomic_Data_Raw.loc\
                        [Genomic_Data_Raw['#Organism Name']==\
                         Name_Genomic]
                
                        Sizes_Raw=Matching_Microbes['Size(Mb)_y']

    #Compute mean of all microbes genetic info found              
    #-----------------------------------------------              
    Sizes=[Size for Size in Sizes_Raw]                            
    MeanSize=np.mean(Sizes)                                       
    #-----------------------------------------------

    return MeanSize
#===========================================                         

#0.4. Shannon Diversity
#===========================================                         
def Shannon(v):            
    Shannon_Diversity=0
    for i in v:                                                      
        if i==0:                                                     
            Shannon_Diversity+=0
        else:                                                        
            Shannon_Diversity+=-i*math.log(i,math.e)
#            Shannon_Diversity+=-i*math.log(i,cosa)
                                                                     
    return Shannon_Diversity
#===========================================

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#1. PARAMETERS FOR INPUT/OUTPUT DATA
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
K=10;L=20
Microbe_Presence=4
Microbes_sorted={};Host_Rank_Shannon={}
Hosts_sorted={};Microbes_size={}

path_plot= '../../Plots_Paper/'

Datasets=\
['S-8_Tot', 'V-10_Tot', 'V-22_Tot', 'V-23_24_Tot', 'V-25_Tot']

Patients_per_Dataset= {'S-8_Tot':107,\
'V-10_Tot':92,'V-22_Tot':467,'V-23_24_Tot':222,'V-25_Tot':883} 
                                                                     
Microbes_per_Dataset={'S-8_Tot':128,\
'V-10_Tot':137,'V-22_Tot':134,'V-23_24_Tot':118,'V-25_Tot':144}
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



#2. DICTIONARIES FOR RANKS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
for Dataset in Datasets:

    #1.1. Adjacency Matrix and taxonomies                            
    #=============================================================== 
    Patients=Patients_per_Dataset[Dataset]                           
    Microbes=Microbes_per_Dataset[Dataset]                           
    path_Out='../../Output_Data/'+ 'Leave_One_Out_' + Dataset + '/'
    path_In= '../../Input_Data/' + 'Leave_One_Out_' + Dataset + '/'
    
    #1.1.1. Taxonomic Information                                    
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    File_Taxon='Names_Microbes_' + Dataset + '.csv'                  
    Taxon_Data_Raw=pd.read_csv('../../Input_Data/'                  
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

    #1.1.4. Choose best likelihood
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    Likelihood_Seed={}                                              
    for Seed in range(1,10,2): 
        Likelihood_File='Dataset_' + Dataset + \
        '_Seed_' + str(Seed) + '_LogLikelihood.txt'
        
        LikelihoodFinal=read_file(path_Out + Likelihood_File)[-1]
        Likelihood_Seed[float(LikelihoodFinal[0])]=Seed              

    Seed=Likelihood_Seed[max(Likelihood_Seed)]                       
    
    File='Dataset_'+Dataset+'_Seed_'+str(Seed)+'_Parameters.txt'
    Data=read_file(path_Out+File)
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #1.1.5. Nestedness for finding ranks
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #1.1.5.1. Adjacency Matrix
    #----------------------------------------------------------------
    AdjacencyMatrix_str=read_file(                              
    '../../Input_Data/Binarized_Datasets/Binarized_'+Dataset+'.txt') 
    AdjacencyMatrix=[]                                               
    for line in AdjacencyMatrix_str:                                 
        new_line=[float(block) for block in line]                    
        AdjacencyMatrix.append(new_line)                             
        
    AdjacencyMatrix=np.asarray(AdjacencyMatrix)
    #----------------------------------------------------------------

    #1.1.5.2. Trasponemos por consistencia con los latent enterotypes
    #---------------------------------------------                   
    AdjacencyMatrix=np.transpose(AdjacencyMatrix)                    
    #---------------------------------------------                   

    #1.1.5.3. Get Nestedness and ranks of microbes/hosts
    #----------------------------------------------------------------
    global rows;   rows=AdjacencyMatrix.shape[0]                     
    global columns;columns=AdjacencyMatrix.shape[1]
    
    NestedMatrix,Nestedness,Dict_row_membership\
    ,Dict_column_membership=Nestedness_fun(AdjacencyMatrix)
    
    #Dict -> Dict[new]=old;Inverted_Dict -> Inverted_Dict[old]=new
    #old=name(as number); new= rank
    
    Inverted_Dict_column_membership =dict(\
    (v_c,k_c) for k_c,v_c in Dict_column_membership.iteritems())

    Inverted_Dict_row_membership =dict(\
    (v_r,k_r) for k_r,v_r in Dict_row_membership.iteritems())        
    #----------------------------------------------------------------


    #2.2.2 Compute Shannon Diversity and relate to hosts rank
    #========================================================
    AdjacencyMatrix_real=read_file('../../Input_Data/'
    +Dataset+'.txt')

    counter_host_shannon=0
    for Microbe_Profile in np.transpose(AdjacencyMatrix_real):
        
        #Shannon diversity
        #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        Microbe_Profile=[float(abund) for abund in Microbe_Profile]
        Shannon_Div_Host=Shannon(Microbe_Profile)
        #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        #Rank of host
        #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        Rank_Host=Inverted_Dict_row_membership[counter_host_shannon]
        #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        #Information into dict
        #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        try:
            Host_Rank_Shannon[Dataset][Rank_Host]=Shannon_Div_Host
            
        except KeyError:
            Host_Rank_Shannon[Dataset]={Rank_Host:Shannon_Div_Host}
        #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        counter_host_shannon+=1
    #========================================================
        
        
    #1.1.3.4. Dictionary species_name-rank
    #----------------------------------------------------------------
    counter_microbes=0                                               
    for microbe in Species_Names:
        #........................................................
        Genomic_Size=Genomic_Size_fun\
        (microbe, Genomic_Data_Raw, Missing_Microbes_Record)
        if np.isnan(Genomic_Size)==True:
            print "holi!"
            print Microbe
        Microbes_size[microbe]=Genomic_Size
        #........................................................
        
        try:                                                         
            Microbes_sorted[microbe][Dataset]=\
            Inverted_Dict_column_membership[counter_microbes]        
            
        except KeyError:                                             
            Microbes_sorted[microbe]={Dataset:\
            Inverted_Dict_column_membership[counter_microbes]}
            
        counter_microbes+=1
    #----------------------------------------------------------------

    #1.1.3.5. Dictionary host_name-rank
    #----------------------------------------------------------------
    counter_hosts=0                                               
    for host in range(Patients_per_Dataset[Dataset]):
        
        try:                                                         
            Hosts_sorted[host][Dataset]=\
                Inverted_Dict_row_membership[counter_hosts]        
            
        except KeyError:                                             
            Hosts_sorted[host]={Dataset:\
            Inverted_Dict_row_membership[counter_hosts]}
            
        counter_hosts+=1
    #----------------------------------------------------------------
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#3. PLOTS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#3.1. Extract lists from dictionary of dictionaries
#=================================================================
Ranks_Hosts=[ [] for i in range(5)]
Shannon_Hosts=[ [] for j in range(5)]

counter_datasets=0
for Dataset in Datasets:
    print Dataset
    for rank in Host_Rank_Shannon[Dataset]:
        Shannon_Diversity=Host_Rank_Shannon[Dataset][rank]
        Ranks_Hosts[counter_datasets].append(rank)
        Shannon_Hosts[counter_datasets].append(Shannon_Diversity)

    counter_datasets+=1
#=================================================================


#3.2. Make Plots
#=================================================================

#3.2.1. Gridspec
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
fig=plt.figure(figsize=(13, 4))                                      
gs = gridspec.GridSpec(1, 5)          
gs.update(left=0.05,right=0.99,bottom=0.14,top=0.88,wspace=0.0)

size_letter=15                                                       
Letters=['a','b','c','d','e']
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



#3.2.2.Run over datasets
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Counter_Plt=0
for Dataset in Datasets:
    
    ax=plt.subplot(gs[0, Counter_Plt])
    ax=plt.scatter(Ranks_Hosts[Counter_Plt],Shannon_Hosts[Counter_Plt], s=6)
    #title,titcks,axes
    #---------------------------------------------
    plt.title(Dataset[:-4])
    plt.ylim(0.45, 3.7)

    Lim_x_up=plt.gca().get_xlim()                                    
    Lim_y_up=plt.gca().get_ylim()
    
    plt.xlabel("Rank", fontsize=15)                        
    if Counter_Plt==0:                                              
        plt.ylabel("Shannon Diversity", fontsize=15)

    else:

        plt.yticks([])
    #---------------------------------------------

    #Letter for caption                                              
    #----------------------------------------------------------------
    plt.text(Lim_x_up[0] + 0.01*(Lim_x_up[1]-Lim_x_up[0]),
        Lim_y_up[1] + 0.01*(Lim_y_up[1]-Lim_y_up[0]), 
        Letters[Counter_Plt],fontsize=size_letter,fontweight='bold')
    #----------------------------------------------------------------
        
    Counter_Plt+=1
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

plt.savefig(path_plot + 'Figure_6.pdf',dpi=300)
plt.show()
#=================================================================

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#4. CSVs
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#CSV for patients
#====================================================================
# Dataframe=pd.DataFrame( {'Name':Microbes_sorted.keys(),
#    'S-8':Ranks_Microbes[0], 'V-10':Ranks_Microbes[1],\
#    'V-22':Ranks_Microbes[2],'V-23-24':Ranks_Microbes[3],\
#    'V-25':Ranks_Microbes[4],'Mean':Ranks_Microbes[5],\
#    'SEM':Ranks_Microbes[6], 'Size(Mb)':Ranks_Microbes[7]})


# Ranks_Hosts_1=[ [] for i in range(5)]
# #for element in Hosts_sorted:
# for element in Hosts_sorted:
#     print element, Hosts_sorted[element]

#CSV for microbes
#====================================================================
Ranks_Microbes=[[] for i in range(8)]

for element in Microbes_sorted:
    
    Mean_Rank=np.mean(Microbes_sorted[element].values())
    Ranks_Microbes[5].append(Mean_Rank)
    SEM=scipy.stats.sem(Microbes_sorted[element].values() ,ddof=0)
    Ranks_Microbes[6].append(SEM)

    Ranks_Microbes[7].append(Microbes_size[element])
    
    index_counter=0
    for Dataset in Datasets:

        try:
            Ranks_Microbes[index_counter].append(\
            Microbes_sorted[element][Dataset])
	    
        except KeyError:
           Ranks_Microbes[index_counter].append(np.nan)
           
        index_counter+=1

Dataframe=pd.DataFrame( {'Name':Microbes_sorted.keys(),
   'S-8':Ranks_Microbes[0], 'V-10':Ranks_Microbes[1],\
   'V-22':Ranks_Microbes[2],'V-23-24':Ranks_Microbes[3],\
   'V-25':Ranks_Microbes[4],'Mean':Ranks_Microbes[5],\
   'SEM':Ranks_Microbes[6], 'Size(Mb)':Ranks_Microbes[7]})

Dataframe=Dataframe.sort_values(by='Mean', axis=0, ascending=True)

Dataframe.to_csv("../../Plots_Paper/Ranks_and_Genome_Size_Choricitos.csv",\
index=False, columns=['Name','S-8','V-10','V-22','V-23-24','V-25',\
                      'Mean','SEM','Size(Mb)'])

#===================================================================

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
