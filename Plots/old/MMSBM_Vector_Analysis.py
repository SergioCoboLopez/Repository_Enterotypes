#15/05/2018. Este codigo hace plots de varias cosas:
#1. De las distribuciones de entropias de Shannon de personas and microbios
#2. Las abundancias de \theta_i y \eta_i
#3. Usuarios y microbios segregados por entropias de Shannon y sus niveles de
#pertenencia a \theta_i y \eta_i

#Aparte de eso, naturalmente, el codigo calcula entropias de Shannon y segrega usuarios y microbios por niveles
#de entropia de shannon

#.-.. --- ... -... --- .-. -... --- -. . ... ... --- -. ..- -. --- ... .-.. .- -.. .-. --- -. . ...
import numpy as np
import math
from scipy.stats import kendalltau
import scipy.cluster
import bisect
import seaborn as sns
from matplotlib import colors
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib
from itertools import groupby
from operator import itemgetter

from collections import Counter
#.-.. --- ... -... --- .-. -... --- -. . ... ... --- -. ..- -. --- ... .-.. .- -.. .-. --- -. . ...

#0.FUNCTIONS                                                                  
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#0.1. Shannon Entropy                                  
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

#0.2. Reading data
#=====================================================
def read_file(inFileName):
        lines = open(inFileName).readlines()
        d = [line.strip().split() for line in lines]
        return d
#=====================================================

#0.2. Generate File Name string
#=====================================================
def get_name(Seed,K,L):
        Name_of_File='0_Seed_' + str(Seed) \
            + '_K_' + str(K) + '_L_' + str(L) + '_Parameters.txt'
        return Name_of_File
#=====================================================

#0.3. Generate Entropy Vectors
#=====================================================                      
def List_Entropy(Data,NGroups):
    Entropy_Data=[];Players_Data={}

    for node in Data:
        #format data to floats
        node=[float(parameter_weight) for parameter_weight in node] 
        Entropy_Data.append( Shannon( node[1:len(node)], NGroups ) )
        Players_Data[node[0]]=\
        [ node[1:len(node)], Shannon( node[1:len(node) ], NGroups ) ]
    return Entropy_Data, Players_Data
#=====================================================

#0.4. Separate theta vectors by ranges in entropies
#=========================================================================
def Split_by_Entropy_Range(Patients_Dict, Entropy_Categories):

    #0.4.1. Take input dictionary and split it by ranges of entropies defined
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    PatientsEntropy={}#New dictionary by ranges of entropy
    for line in Patients_Dict:
        Entropy=Patients_Dict[line][1]
        i = bisect.bisect(Entropy_Categories,Entropy ) #cut by categories
        Patients_Dict[line].append(i-1)	#append category to input dictionary
        
        #0.4.1.1. Build new dictionary
        #--------------------------------------------------------------------
        try:
            PatientsEntropy[ Patients_Dict[line][2] ].\
                append( Patients_Dict[line][0:2]  )
        except KeyError:
            PatientsEntropy[ Patients_Dict[line][2] ]= \
                        [ Patients_Dict[line][0:2] ]
        #--------------------------------------------------------------------
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    
    #0.4.2. Sort vectors of thetas by decreasing Shannon entropy
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #0.4.2.1. Declare arrays of patient thetas sorted by entropy
    #----------------------------------------
    Patients=[ [] for i in range(len(PatientsEntropy))]
    #----------------------------------------

    #0.4.2.2. Create arrays of thetas
    #--------------------------------------------------------------------------
    counter=0
    for EntropyRange in PatientsEntropy:

        Patients[counter]=[thetas[0] for thetas in PatientsEntropy[EntropyRange]  ]
        Patients[counter]=np.array( Patients[counter] )
        
        counter+=1
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    return Patients
#==============================================================================

def Cluster_Matrix(Matrix):

    #Clustering vector
    #--------------------------------------------------------------------------
    Clustering=scipy.cluster.hierarchy.fclusterdata\
                (Matrix,0.0,method='average')
    #--------------------------------------------------------------------------

    #Dictionary cluster_number: node
    #----------------------------------
    Clustering_Dict={}
    for line in range(0,len(Clustering)):
        try:
            Clustering_Dict[ Clustering[line]  ].append( line  )
        except KeyError:
            Clustering_Dict[ Clustering[line]  ]=[line]
    #----------------------------------
    
    #Rearrange matrix with Dictionary
    #----------------------------------
    Matrix_Clustered=np.zeros( Matrix.shape )

    counter=0;New_Order=[];
    Clusters=[]#-->Each element is the cluster to which the analogous element in Matrix_Clustered belongs
    for key in Clustering_Dict:
        for node in Clustering_Dict[key]:
            
            Clusters.append(key)
            
            New_Order.append(node)
            
            Matrix_Clustered[counter]=Matrix[node]
            counter+=1
    #----------------------------------
    
    return Matrix_Clustered,Clusters,New_Order
#==============================================================================


#0.5. Cluster Patients (Microbes) and reorder matrices
#==============================================================================
def Cluster_Magic(Patients_or_Microbes, OnlyRows = 0):

    #0.5.1. Declare variables
    #----------------------------------
    Entropy_Rank=len(Patients_or_Microbes)
    New_Order_Rows=[0]*Entropy_Rank
    New_Order_Cols=[0]*Entropy_Rank
    Clusters_Rows=[0]*Entropy_Rank
    Clusters_Cols=[0]*Entropy_Rank
    Patients_or_Microbes_Clustered=[0]*Entropy_Rank
    #----------------------------------
    
    for Rank in range(Entropy_Rank):

        #Clustering by rows
        #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        
        #Excepction if there's only one Patient/Microbe
        #----------------------------------
        if len(Patients_or_Microbes[Rank])==1:
            Transposed=np.transpose(Patients_or_Microbes[Rank])
            Clustered_by_cols,Clusters_Cols[Rank],New_Order_Cols[Rank]=\
                                            Cluster_Matrix( Transposed )
            Patients_or_Microbes_Clustered[Rank]=\
                    np.transpose(Clustered_by_cols)
            continue
        #----------------------------------

        #Clustering
        #----------------------------------------------------------------------
        Clustered_by_rows,Clusters_Rows[Rank],New_Order_Rows[Rank]=\
                        Cluster_Matrix( Patients_or_Microbes[Rank] )
        #----------------------------------------------------------------------
        
        #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        
        #Clustering by cols
        #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        #Transpose clustered-by-rows matrix
        #-------------------------------------------
        Transposed=np.transpose( Clustered_by_rows )
        #-------------------------------------------

        #Clustering
        #----------------------------------------------------------------------
        Clustered_by_cols,Clusters_Cols[Rank],New_Order_Cols[Rank]=\
                                        Cluster_Matrix( Transposed )
        #----------------------------------------------------------------------

        #Detranspose
        #-------------------------------------------
        Patients_or_Microbes_Clustered[Rank]=np.transpose(Clustered_by_cols)
        #-------------------------------------------
        #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    return Patients_or_Microbes_Clustered,Clusters_Rows,Clusters_Cols,\
        New_Order_Rows, New_Order_Cols
#==============================================================================



#0.5. Cluster Patients (Microbes) and reorder matrices
#=============================================================================
# def Cluster_Magic_Old(Patients_or_Microbes):

#     #0.5.1. Declare variables
#     #----------------------------------
#     Entropy_Rank=len(Patients_or_Microbes)
#     Clustering_Rows=[0]*Entropy_Rank
#     Clustering_Cols=[0]*Entropy_Rank
#     Patients_or_Microbes_Clustered_by_row=[0]*Entropy_Rank
#     Patients_or_Microbes_Clustered=[0]*Entropy_Rank
#     #----------------------------------
    
#     for Rank in range(Entropy_Rank):

#         #Clustering by rows
#         #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#         #Excepction if there's only one Patient/Microbe
#         #----------------------------------
#         if len(Patients_or_Microbes[Rank])==1:
# #            Patients_or_Microbes_Clustered_by_row[Rank]=Patients_or_Microbes[Rank]
#             Patients_or_Microbes_Clustered[Rank]=Patients_or_Microbes[Rank]
#             continue
#         #----------------------------------

#         #Clustering vector
#         #-------------------------------------------------------------------
#         Clustering_Rows[Rank]=scipy.cluster.hierarchy.fclusterdata\
#                                (Patients_or_Microbes[Rank],0.0,method='average')
#         #-------------------------------------------------------------------

#         #Dictionary cluster_number: node
#         #----------------------------------
#         Clustering_RowsDict={}
#         for line in range(0,len(Clustering_Rows[Rank])):
#             try:
#                 Clustering_RowsDict[ Clustering_Rows[Rank][line]  ].append( line  )
#             except KeyError:
#                 Clustering_RowsDict[ Clustering_Rows[Rank][line]  ]=[line]
#         #----------------------------------

#         #Rearrange matrix with Dictionary
#         #----------------------------------
#         Patients_or_Microbes_Clustered_by_row[Rank]=np.zeros(Patients_or_Microbes[Rank].shape)
        
#         counter_row=0
#         for key in Clustering_RowsDict:
#             for node in Clustering_RowsDict[key]:
#                 Patients_or_Microbes_Clustered_by_row[Rank][counter_row]=Patients_or_Microbes[Rank][node]
#                 counter_row+=1
#         #----------------------------------
#         #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        
#         #Clustering by cols
#         #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#         Transposed=np.transpose( Patients_or_Microbes_Clustered_by_row[Rank] )

#         Clustering_Cols[Rank]=scipy.cluster.hierarchy.fclusterdata(Transposed,0.0,method='average')

        
#         #Dictionary cluster_number: node
#         #----------------------------------
#         Clustering_Cols_Dict={}
#         for line in range(0,len(Clustering_Cols[Rank])):
#             try:
#                 Clustering_Cols_Dict[ Clustering_Cols[Rank][line]  ].append( line  )
#             except KeyError:
#                 Clustering_Cols_Dict[ Clustering_Cols[Rank][line]  ]=[line]
#         #----------------------------------

#         #Rearrange matrix with Dictionary
#         #----------------------------------
#         Patients_or_Microbes_Clustered[Rank]=np.zeros( Transposed.shape )
#         counter_col=0
        
#         for key in Clustering_Cols_Dict:
#             for node in Clustering_Cols_Dict[key]:
#                 Patients_or_Microbes_Clustered[Rank][counter_col]=Transposed[node]
#                 counter_col+=1
#         #----------------------------------

#         Patients_or_Microbes_Clustered[Rank]=np.transpose(Patients_or_Microbes_Clustered[Rank])
#         #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        
#     return Patients_or_Microbes_Clustered_by_row
#=============================================================================

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#1. READ DATA AND GET MATRICES FOR HEATMAP
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#1.1. File
#=================================================
path='../../Output_Data/ToyDataset/ModelParameters/'
K=10;L=20;
Patients=370;Microbes=1220

#Look for the file with the right seed
#-------------------------------------
for Seed in range(1111,9999+1,2222):
    File=get_name(Seed,K,L)
    try:
        Data=read_file(path+File)
        break
    except IOError:
        continue
#-------------------------------------

#=================================================

#1.2. Compute thetas, etas and aggregated
#=================================================
theta=[];contador_theta=0
eta=[];contador_eta=0

for line in Data[:Patients]:
    line.insert(0, contador_theta )
    contador_theta+=1
    
theta=Data[:Patients]
    
for line in Data[Patients+1:Patients+1+Microbes]:
    line.insert(0, contador_eta )
    contador_eta+=1
    
eta=Data[Patients+1:Patients+1+Microbes]

theta_aggregate=[0]*K
for lineT in theta:
    for i in range(K):
        theta_aggregate[i]+=float(lineT[i+1])

eta_aggregate=[0]*L
for lineE in eta:    
    for j in range(L):
        eta_aggregate[j]+=float(lineE[j+1])    
#=================================================

#1.3. Get Entropy Vectors                     
#======================================================
EntropyPatients,Patients_Dict=List_Entropy( theta, K );
EntropyMicrobes,Microbes_Dict=List_Entropy( eta  , L );
#======================================================

#1.4. Get parameters for vertical lines
#=================================================
Entropy_Categories_K=[]
for k in range(1,K+1):
    EntropyBox_K=-math.log(1/float(k),K)
    Entropy_Categories_K.append(EntropyBox_K)

Entropy_Categories_L=[] 
for l in range(1,L+1):
    EntropyBox_L=-math.log(1/float(l),L)
    Entropy_Categories_L.append(EntropyBox_L)
#=================================================

#1.5. Split Dictionaries by entropies
#==================================================================
Patients=Split_by_Entropy_Range(Patients_Dict,Entropy_Categories_K)
Microbes=Split_by_Entropy_Range(Microbes_Dict,Entropy_Categories_L)
#==================================================================

#1.6. Clustering patients and microbes
#==============================================================================
Patients_Clustered,Clusters_Patients_Row,Clusters_Patients_Col,\
New_Order_Rows_Patients, New_Order_Cols_Patients=Cluster_Magic(Patients)
Microbes_Clustered,Clusters_Microbes_Row,Clusters_Microbes_Col, \
New_Order_Rows_Microbes, New_Order_Cols_Microbes=Cluster_Magic(Microbes)
#==============================================================================



#1.7. METADATA
#==============================================================================

#1.7.1. Read Metadata
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
path_meta= '../../Input_Data/taxonomies_for_toy_dataset.txt'
Metadata=read_file(path_meta)
contador_microbes=0
Microbes_Genus={}
Genus_Microbes={}
for line in Metadata:
    #Get genus and remove quotes
    #-----------------------------------------
    genus_with_quote=line[0].split('_')[0] #--
    genus=genus_with_quote[1:]             #--
    #-----------------------------------------
    
    Microbes_Genus[contador_microbes]=genus
    try:
        Genus_Microbes[genus].append(contador_microbes)
    except KeyError:
        Genus_Microbes[genus]=[contador_microbes]

    contador_microbes+=1
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#1.7.2.Build Reciprocal Dictionary. key: Shannon. Value: eta vector, Microbe_Id
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Microbes_Dict_Reciproco={}
for key in Microbes_Dict:
    Microbe_Id=key;Membership_vector=Microbes_Dict[key][0];
    Shannon_Microbe=Microbes_Dict[key][1]
    Microbes_Dict_Reciproco[Shannon_Microbe]=(Membership_vector, Microbe_Id)
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#1.7.3. Given the clustered microbe matrices, create a list (of lists) with the original ids of the microbes in the same order.
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Id_Original= [ [] for i in range(0, len(Microbes) ) ]

#Run over clustered microbial matrices
for entropy_id in range(0,len(Microbes_Clustered)):
    #Run over rows in each matrix
    for microbe in Microbes_Clustered[entropy_id]: 
        
        #Coges cada fila de la matriz, calculas la entropia de Shannon y la \
        #buscas en el diccionario reciproco
        try:
            (Id_Original[entropy_id]).append( Microbes_Dict_Reciproco[Shannon(microbe,L)][1] )
        except KeyError:
        #Diferencias numericas en entropia de Shannon clusterizada y original  
        #Coges shannon del vector y buscas el valor mas parecido en diccionario
        #reciproco
            #------------------------------------------------------------------
            idx = (np.abs(np.asarray(Microbes_Dict_Reciproco.keys()) - Shannon(microbe,L))).argmin()
            (Id_Original[entropy_id]).append( Microbes_Dict_Reciproco[Microbes_Dict_Reciproco.keys()[idx]][1] )
            #------------------------------------------------------------------
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#1.7.4. Build a dictionary such that - key: genus; value: integer representing category of the genus 
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Genus_Categories={}
for key in Genus_Microbes:#Correr por diccionario con key=Genus, value=microbe

    Species_per_Genus=len(Genus_Microbes[key])
        
    if Species_per_Genus==1:
        Genus_Categories[key]=0

    elif Species_per_Genus>1 and Species_per_Genus<=4:
        Genus_Categories[key]=1
    elif Species_per_Genus>4 and Species_per_Genus<=20:
        Genus_Categories[key]=2        
    
Genus_Categories['Mycoplasma']=3;Genus_Categories['Streptococcus']=4;
Genus_Categories['Lactobacillus']=5;Genus_Categories['Clostridium']=6
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#1.7.5. Given the clustered microbe matrices, generate a list of lists with
       #the category of genuses to which the corresponding microbe belongs
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Genus_Numbers=[ [] for _ in range(len(Id_Original)) ]
yticks=[ [] for _ in range(len(Id_Original)) ]

for heatmap in range(len(Id_Original)):
    for microbe in Id_Original[heatmap]:
        yticks[heatmap].append(Microbes_Genus[microbe])
        #Para cada microbio en la id_list, busco primero el genus y, con el
        #genus, busco el entero que le toca
        Genus_Numbers[heatmap].append(\
        Genus_Categories[ Microbes_Genus[microbe] ] )
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


#1.7.6. Statistics on the metadata: percentage of categories present in each
#entropy rank and similar
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
counter=[0,0,0,0,0,0,0]
for line in Genus_Numbers:
    for microbe in line:
        counter[microbe]+=1
print counter

for line in Genus_Numbers:
    print Counter(line)
    for key in Counter(line):
        print key, Counter(line)[key],len(line), float(Counter(line)[key])/float(len(line))\
        ,counter[key], float(Counter(line)[key])/float(counter[key])
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
#==============================================================================

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#2. PLOTS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Plots=['ShannonEntropy', 'theta+eta_Occupation', 'Heatmaps_Patients','Heatmaps_Microbes',\
       'Heatmap_Microbes_with_Metadata','nothing']
Information=Plots[4]

#Fontsizes
#:::::::::::::
size_eje=15
size_ticks=13
size_title=15
size_letter=17
#:::::::::::::

pathFig='/export/home/shared/Projects/Microbiome/Plots/'

color1=(float(141)/float(255),float(160)/float(255),float(153)/float(203),1)

if Information==Plots[0]:
    
    fig=plt.figure(figsize=(10,5.3))
    gs = gridspec.GridSpec(1, 2)
    gs.update(left=0.08, right=0.99, wspace=0.25,hspace=0.5)

    Bins=10
    color=(float(28)/float(255),float(144)/float(255),float(153)/float(255),1)

    ticks_posx=[0,0.25,0.5,0.75,1]
    ticks_posy=[tick for tick in range(0,5,1)]

    #(0,0) PLOT
    #===================================================
    ax_t0 = plt.subplot(gs[0, 0])
    sns.distplot(EntropyPatients,bins=K, kde=False, norm_hist=True,color=color)
    print sum(EntropyPatients)
    #Lineas verticales
    #-----------------------------------------------
    for k_line in Entropy_Categories_K:
        plt.axvline(x=k_line, color='gray', linestyle='--')
        if k_line > max(EntropyPatients):
            break
    #-----------------------------------------------
    ax_t0.set_xticks(ticks_posx)
    ax_t0.set_xticklabels(ticks_posx,fontsize=size_ticks)
    ax_t0.set_yticks(ticks_posy)
    ax_t0.set_yticklabels(ticks_posy,fontsize=size_ticks)
    ax_t0.set_ylim(0,5)
    ax_t0.set_title(r'$\theta$ (Patients), K='+str(K),fontsize=size_title)
    ax_t0.set_xlabel('H',fontsize=size_eje)
    ax_t0.set_ylabel('Probability Density Function',fontsize=size_eje)
    #===================================================

    #(1,0) PLOT                                                                                                  
    #===================================================
    ax_t1 = plt.subplot(gs[0, 1])
    sns.distplot(EntropyMicrobes,bins=L, kde=False, norm_hist=True, color=color)
    #Lineas verticales
    #-----------------------------------------------
    for l_line in Entropy_Categories_L:
        plt.axvline(x=l_line, color='gray', linestyle='--')
        if l_line > max(EntropyMicrobes):
            break
    #-----------------------------------------------
    ax_t1.set_xticks(ticks_posx)
    ax_t1.set_xticklabels(ticks_posx,fontsize=size_ticks)
    ax_t1.set_yticks([])
    ax_t1.set_ylim(0,5)
    ax_t1.set_title(r'$\eta$ Microbes L='+str(L),fontsize=size_title)
    ax_t1.set_xlabel('H',fontsize=size_eje)
    #===================================================

    #===================================================
    path_out='/export/home/shared/Projects/Microbiome/Plots/'
    # plt.savefig(path_out + 'Shannon_K=' + str(K) + '_L=' + str(L) +'.pdf')
    # plt.savefig(path_out + 'ShannonEntropy_pngs/Shannon_K=' + str(K) + '_L=' + str(L) +'.png')

    plt.show()
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++                        

elif Information==Plots[1]:

    fig=plt.figure(figsize=(10,5.3))
    gs = gridspec.GridSpec(1, 2)
    gs.update(left=0.08, right=0.99, wspace=0.25,hspace=0.5)
    
    #(0,0) PLOT
    #===================================
    ax0 = plt.subplot(gs[0, 0])
    x0 = [index0 for index0 in range(K)]
    y0 = theta_aggregate
    ax0.set_title(r'$\theta_i$ distributions',fontsize=size_title)
        
    sns.barplot(x0, y0, color=color1)
    #===================================

    #(0,1) PLOT
    #===================================
    ax1 = plt.subplot(gs[0, 1])
    x1 = [index1 for index1 in range(L)]
    y1 = eta_aggregate
    ax1.set_title(r'$\eta_i$ distributions',fontsize=size_title)
    sns.barplot(x1, y1, color=color1)
    #===================================

    plt.show()
#==============================================================================
elif Information==Plots[2]:

    columnas=len(Patients) #Numero de heatmaps que vamos a tener
    fig=plt.figure(figsize=(5*(columnas),10))

    #Definir los width ratios de acuerdo al numero de heatmaps
    #.........................................................
    Ancho_Plots=[1]*columnas
    Ancho_Plots.append( 1/float(columnas+10 ) )
    #.........................................................
    
    gs = gridspec.GridSpec(2, columnas+1, width_ratios=Ancho_Plots, height_ratios=[0.2,0.8]   )
    gs.update(left=0.038, right=0.96, wspace=0.13,hspace=0.05, bottom=0.05, top=0.96)


    #Bucle para pintar Heatmaps
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    for plot in range( columnas ):

        #(:,0) PLOT (Histograma)
        #----------------------------------------------------------------------
        ax0 = plt.subplot(gs[0, plot])
        x0 = [index0 for index0 in range(K)]

        #Change column order
        #----------------------------------------------
        y0 = sum(Patients[plot])
        y1=[]
        for New_Index in New_Order_Cols_Patients[plot]:
            y1.append( y0[New_Index] ) 
        #----------------------------------------------

        sns.barplot(x0, y1, color=color1)

        #Axis
        #----------------------------------------------------------------------
        ax0.set_title(str("%.2f"%Entropy_Categories_K[plot])+\
                      '< H < '+str("%.2f" % Entropy_Categories_K[plot+1]))
        ax0.set(xticklabels=New_Order_Cols_Patients[plot])
        #----------------------------------------------------------------------
        
        #(:,1) PLOT (Heatmap)
        #----------------------------------------------------------------------
        ax = plt.subplot(gs[1, plot])
        
        Heatmap=Patients_Clustered[plot]
        
        sns.heatmap(Heatmap,cmap="RdBu_r",vmin=0,vmax=1,cbar=False,xticklabels=New_Order_Cols_Patients[plot],\
                    yticklabels=False)
        ax.set_xlabel('K',fontsize=size_eje);

        #Ejes
        #..................................................
        if plot==0:
            ax.set_ylabel('Patients(i)',fontsize=size_eje);
        #..................................................

        #----------------------------------------------------------------------
        
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
            
        #PLOT: COLORBAR
        #======================================================================
	ax1=plt.subplot(gs[1,-1 ]  )
	colbar_ticks=[0,0.5,1]
	cb1 = matplotlib.colorbar.ColorbarBase(ax1, cmap="RdBu_r",\
                                               ticks=colbar_ticks)
        cb1.set_ticks(colbar_ticks)
        cb1.set_ticklabels(colbar_ticks)
        cb1.ax.tick_params(labelsize=size_ticks)
        cb1.outline.set_visible(False)
        cb1.set_label(r"$\theta^K_i$",fontsize=size_ticks)
        #======================================================================
    # plt.savefig(pathFig+\
    # 'Heatmap_Thetas'+'_K_'+str(K)+'_L_'+ str(L)+'.png', dpi=300)
    plt.show()
#==============================================================================

#==============================================================================
elif Information==Plots[3]:

    columnas=len(Microbes) #Numero de heatmaps que vamos a tener
    fig=plt.figure(figsize=(5*(columnas),10))

    #Definir los width ratios de acuerdo al numero de heatmaps
    #.........................................................
    Ancho_Plots=[1]*columnas
    Ancho_Plots.append( 1/float(columnas+10 ) )
    #.........................................................
    
    gs = gridspec.GridSpec(2, columnas+1, width_ratios=Ancho_Plots, height_ratios=[0.2,0.8]   )
    gs.update(left=0.025, right=0.96, wspace=0.1,hspace=0.05, bottom=0.05, top=0.96)

    #Bucle para pintar Heatmaps
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    for plot in range( columnas ):

        #(:,0) PLOT (Histograma)
        #----------------------------------------------------------------------
        ax0 = plt.subplot(gs[0, plot])
        x0 = [index0 for index0 in range(L)]

        #Change column order
        #----------------------------------------------
        y0 = sum(Microbes[plot])
        y1=[]
        for New_Index in New_Order_Cols_Microbes[plot]:
            y1.append( y0[New_Index] ) 
        #----------------------------------------------
            
        sns.barplot(x0, y1, color=color1)

        #Ejes
        #----------------------------------------------
        ax0.set_title(str("%.2f"%Entropy_Categories_L[plot])+'< H < '+str("%.2f" % Entropy_Categories_L[plot+1]))
        ax0.set(xticklabels=New_Order_Cols_Microbes[plot])
        #----------------------------------------------

        #----------------------------------------------------------------------


        #(:,1) PLOT (Heatmap)
        #----------------------------------------------------------------------
        ax = plt.subplot(gs[1, plot])
        
        Heatmap=Microbes_Clustered[plot]

        sns.heatmap(Heatmap,cmap="RdBu_r",cbar=False,vmin=0,vmax=1,xticklabels=New_Order_Cols_Microbes[plot], \
                    yticklabels=False)

        #Ejes
        #..................................................
        ax.set_xlabel('L',fontsize=size_eje);
        
        if plot==0:
            ax.set_ylabel('Microbes(j)',fontsize=size_eje);

        # ytick_pos=[ n+0.5 for n in range(len(yticks[plot])) ]
        # ax.set_yticks( ytick_pos )
        # ax.set_yticklabels(yticks[plot],fontsize=8)
        #..................................................

        #----------------------------------------------------------------------
        
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
            
        #PLOT: COLORBAR
        #======================================================================
	ax1=plt.subplot(gs[1,-1 ]  )
	colbar_ticks=[0,0.5,1]
	cb1 = matplotlib.colorbar.ColorbarBase(\
              ax1, cmap="RdBu_r",ticks=colbar_ticks)
        cb1.set_ticks(colbar_ticks)
        cb1.set_ticklabels(colbar_ticks)
        cb1.ax.tick_params(labelsize=size_ticks)
        cb1.outline.set_visible(False)
        cb1.set_label(r"$\eta^L_j$",fontsize=size_ticks)
        #======================================================================

    # plt.savefig(pathFig+'Heatmap_Etas'+ '_K_' + str(K) + '_L_' + str(L) + '.png', dpi=300)
    plt.show()    
#==============================================================================

#==============================================================================
elif Information==Plots[4]:

    #Parameters on the gridspec plot
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    columnas=2*len(Microbes) #2 columnas por cada rango de entropia considerado
    fig=plt.figure(figsize=(5*len(Microbes),10))

    #Definir los width ratios de acuerdo al numero de heatmaps
    #.........................................................
    #Ancho del metadata, del heatmap y multiplicado por rangos de entropia
    Ancho_Plots=[0.08,0.92]*len(Microbes) 
    Ancho_Plots.append( 1/float(columnas+10 ) ) #Ancho de la colorbar
    #.........................................................
    
    gs = gridspec.GridSpec(2, columnas+1, width_ratios=Ancho_Plots, height_ratios=[0.15,0.85]   )
    gs.update(left=0.025, right=0.96, wspace=0.1,hspace=0.05, bottom=0.05, top=0.96)
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    

    #Bucle para pintar Heatmaps
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    for plot in range( 0,columnas,2 ):

        Index=int(0.5*plot)#An index for plot position,another for entropy rank

        #(:,0) PLOT (Histograma)
        #----------------------------------------------------------------------
        ax0 = plt.subplot(gs[0, plot+1])
        x0 = [index0 for index0 in range(L)]

        #Change column order
        #______________________________________________
        y0 = sum(Microbes[Index])
        y1=[]
        for New_Index in New_Order_Cols_Microbes[Index]:
            y1.append( y0[New_Index] )
        #______________________________________________

        #Plot
        #________________________________
        sns.barplot(x0, y1, color=color1)
        #________________________________
        
        #Ejes
        #___________________________________________
        ax0.set_title(str("%.2f"%Entropy_Categories_L[Index])+'< H <'+str("%.2f" % Entropy_Categories_L[Index+1]))
        ax0.set(xticklabels=New_Order_Cols_Microbes[Index])
        #___________________________________________

        #----------------------------------------------------------------------


        #(:,1) PLOT (Heatmaps)
        #----------------------------------------------------------------------

        #Heatmap de metadata
        #______________________________________________________________________
        ax1 = plt.subplot(gs[1, plot])
        
        Heatmap1=np.transpose( np.matrix(Genus_Numbers[Index]) )

        sns.heatmap(Heatmap1,cmap="Set1",cbar=False,vmin=0,vmax=8, annot=False,xticklabels=False,yticklabels=False)
        #Ejes
        #..................................................
        if plot==0:
            ax1.set_ylabel('Microbes(j)',fontsize=size_eje);
        #..................................................
        #______________________________________________________________________

        
        #Heatmap principal
        #______________________________________________________________________
        ax2 = plt.subplot(gs[1, plot+1])

        Heatmap2=Microbes_Clustered[Index]

        sns.heatmap(Heatmap2,cmap="RdBu_r",cbar=False,vmin=0,vmax=1,xticklabels=New_Order_Cols_Microbes[Index], \
                    yticklabels=False)

        ax2.set_xlabel('L',fontsize=size_eje);
        #______________________________________________________________________

        #----------------------------------------------------------------------

        
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

            
        #PLOT: COLORBAR
        #======================================================================
	ax1=plt.subplot(gs[1,-1 ]  )
	colbar_ticks=[0,0.5,1]
	cb1 = matplotlib.colorbar.ColorbarBase(ax1, cmap="RdBu_r",ticks=colbar_ticks)
        cb1.set_ticks(colbar_ticks)
        cb1.set_ticklabels(colbar_ticks)
        cb1.ax.tick_params(labelsize=size_ticks)
        cb1.outline.set_visible(False)
        cb1.set_label(r"$\eta^L_j$",fontsize=size_ticks)
        #======================================================================

    # plt.savefig(pathFig+\
    # 'Heatmap_Etas'+ '_K_' + str(K) + '_L_' + str(L) + 'metadata.png', dpi=300)
    plt.show()    
#==============================================================================

#====================
else:
    print "No plots"
#====================


plt.show()
