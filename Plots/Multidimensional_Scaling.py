#7/05/2019. Codigo para hacer un multidimensional scaling de los datos de los vectores de membership
#de los microbios

import sys
import os
import numpy as np
import pandas as pd
import random as rnd
import copy as cp
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
from matplotlib import colors
import seaborn as sns
import scipy.cluster
from scipy import stats
from scipy.spatial import distance
import math
import operator
from collections import Counter
import time
from sklearn.manifold import MDS

#0. FUNCTIONS
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#0.1. Read Data
#===============================================
def read_file(inFileName):
    lines = open(inFileName).readlines()
    d = [line.strip().split() for line in lines]
    return d
#=============================================== 

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



#1. PARAMETERS FOR INPUT/OUTPUT DATA
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#1.1. EXTERNAL PARAMETERS
#========================
Dataset=sys.argv[1]   #==
#========================

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
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Taxonomic_Orders=['Kingdom','Phylum','Class','Order','Family','Genus','Species']

Taxonomic_level=5
Taxonomic_Order=Taxonomic_Orders[Taxonomic_level]

File_Taxon='Names_Microbes_' + Dataset + '.csv'
Taxon_Data=pd.read_csv(path_In+File_Taxon,header=None)
Taxon_Data=Taxon_Data[Taxon_Data.columns[Taxonomic_level]] 
Taxon_List=Taxon_Data.tolist()
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
for microbe in range(len(eta)):
        eta[microbe] = [float(membership) for membership in eta[microbe] ]
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#4. Do Multidimensional scaling
#+++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++
seed = np.random.RandomState(seed=1)
embedding = MDS(n_components=3, random_state=seed)
print embedding
eta_transformed=embedding.fit_transform(eta)
#+++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++

#5. GROUP MICROBES BY TAXONOMIC ORDER
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Taxonomy_Microbes_Dict={}
for microbe_1 in range( len(Taxon_List) ):
        Id_eta_tuple=( microbe_1, eta_transformed[microbe_1] )

        try:
            Taxonomy_Microbes_Dict[ Taxon_List[microbe_1] ].append( Id_eta_tuple )
        except KeyError:
            Taxonomy_Microbes_Dict[ Taxon_List[microbe_1] ]=[ Id_eta_tuple ]
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#6. CLEAN DICTIONARY OF MONO/BI-CLADE SPECIES
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#We look in the taxonomy dictionary for clades that host more than two
#microbial species and save them to the new dictionary

Purged_Taxonomy_Microbes_Dict={}
for clade in Taxonomy_Microbes_Dict:
    print clade
    print Taxonomy_Microbes_Dict[clade]
    
    if len(Taxonomy_Microbes_Dict[clade])>2:
        Purged_Taxonomy_Microbes_Dict[clade]=Taxonomy_Microbes_Dict[clade]
    
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#7. PLOTS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#7.1. Colors for scatterplot
#========================================================================

#7.1.1. Choose the colormap to chop
#::::::::::::::::::::::::::::::::::::::::::::::
Number_of_Colors=len( Purged_Taxonomy_Microbes_Dict.keys() )

print Number_of_Colors
    
markers=["o", "^", "v", "s", "D", "X", "P"]
colors_microbes = cm.Set1(np.linspace(0, 1, 9))

Plot_Id=[]


for marker in markers:
    for color in colors_microbes:
        Plot_Id.append( (color, marker) )

Names_of_Clades = Purged_Taxonomy_Microbes_Dict.keys()
Dict_marker_color={}

counter_Id=0
for clade in Names_of_Clades:
    Dict_marker_color[clade]=Plot_Id[counter_Id]
    counter_Id+=1

#========================================================================

Number_of_Dimensions=2

if Number_of_Dimensions==2:
    #7.3. 2-dimensional scatterplot
    #============================================================================
    fig = plt.figure(figsize=(12,10))
    ax2d=fig.add_subplot(111)
    for clade in Purged_Taxonomy_Microbes_Dict:
        for microbe in Purged_Taxonomy_Microbes_Dict[clade]:
            plt.scatter(\
            microbe[1][0], microbe[1][1], color=Dict_marker_color[clade][0], marker=Dict_marker_color[clade][1],s=40)

    plt.title("Multidimensional Scaling " + str(Dataset)+ " " + str(Taxonomic_Order) + " 2 dimensions")
    ax2d.set_xlabel("X"); ax2d.set_ylabel("Y")

    plt.savefig(\
    '../../Plots/First_Path/Multdimensional_Scaling_2D_' + Dataset + '_' + Taxonomic_Order + '.pdf',dpi=300)
    plt.show()
    #============================================================================

elif Number_of_Dimensions==3:
    
    #7.4. 3-dimensional scatterplot
    #=========================================================
    fig = plt.figure(figsize=(15,12))
    ax3d = Axes3D(fig)        

    for clade in Purged_Taxonomy_Microbes_Dict:
        for microbe in Purged_Taxonomy_Microbes_Dict[clade]:
            ax3d.scatter(microbe[1][0], microbe[1][1], microbe[1][2],\
                         color=Dict_marker_color[clade][0], marker=Dict_marker_color[clade][1], s=50)

    plt.title("Multidimensional Scaling " + str(Dataset)+ " " + str(Taxonomic_Order) + " 3 dimensions")
    ax3d.set_xlabel("X"); ax3d.set_ylabel("Y");ax3d.set_zlabel("Z")


    plt.savefig(\
    '../../Plots/First_Path/Multdimensional_Scaling_3D_' + Dataset + '_' + Taxonomic_Order + '.pdf',dpi=300)
    plt.show()
    #=========================================================


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
