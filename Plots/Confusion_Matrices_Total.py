#25/10/2018. Este script genera todas las confusion matrices para un dataset
#determinado.

import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import scipy.cluster



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


#0.2. Clean File
#======================================================
def clean_data(data):
    scores=[]
    for line in data:
        if len(line)==0:
            continue
        line=map(float,line)
        scores.append(line)

    return scores
#======================================================

#0.3. Create Confusion Matrix
#======================================================
def Confusion_Matrix_func(Confusion_Matrices_dict,Patient):

    Confusion_Matrix_Instance=np.zeros((cols,rows))
    
    if Patient==None:
        #Do total Confusion Matrix
        
        #Create raw Confusion Matrix
        #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        for Patient in Confusion_Matrices_dict:
            for Fold in Confusion_Matrices_dict[Patient]:
                Confusion_Matrix_Instance+=\
                        Confusion_Matrices_dict[Patient][Fold]
        #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    else:
        for Fold in Confusion_Matrices_dict[Patient]:
            Confusion_Matrix_Instance+=\
                        Confusion_Matrices_dict[Patient][Fold]


    #Normalize
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    Normalization=sum(Confusion_Matrix_Instance)
    for row in range(0,rows):
        for col in range(0,cols):
            Confusion_Matrix_Instance[row][col]=\
            Confusion_Matrix_Instance[row][col]/Normalization[col]    
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
            
    return Confusion_Matrix_Instance
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++


#1.GENERATE CONFUSION MATRICES
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#1.1. Parameters for reading data
#========================================================

Seed=np.random.seed(1111)#--->Para escoger la muestra de 9 pacientes

#1.1.1. EXTERNAL PARAMETERS
#---------------------
Dataset=sys.argv[1]#--
Number_of_Links=3  #--
#---------------------

#Paths for files
#------------------------------------------------------------
path_Out='../../Output_Data/'+ 'Leave_One_Out_' + Dataset + '/'
path_In= '../../Input_Data/' + 'Leave_One_Out_' + Dataset + '/'
#------------------------------------------------------------

#Read best results from best likelihoods
#--------------------------------------------
Results=read_file(path_In + 'Best_Seeds.txt')
#--------------------------------------------

#1.2. Internal parameters and structures
#==========================================================================
Confusion_Matrices_dict={}
rows=Number_of_Links;cols=Number_of_Links;
Link_Counter=[0]*Number_of_Links#->We count links present in the 30 pat. sample
#==========================================================================

#1.3. Generate+write results
#==============================================================================
Confusion_Matrix_Total=np.zeros((rows,cols))
for line in Results:

    #1.3.1. Read Input files
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    Patient=line[0];Fold=line[1];Seed=line[2]
    Parameters = str(Patient) + "_Fold_" + str(Fold) +  "_Seed_" + str(Seed)

    Scores = path_Out +  Parameters +  '_scores.txt'
    raw_scores=read_file(Scores)
    Check=read_file(path_In + 'P_' + Patient + '_F_'+Fold+ '_TestCheck.txt')
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #1.3.2. Read and compute MMSBM scores
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #1.3.2.1. Clean scores file
    #----------------------------
    scores=clean_data(raw_scores)
    #----------------------------

    #1.3.2.2. Collect scores in dictionary
    #------------------------------------------
    scores_dict={}
    Another_scores_dict={}
    
    for line in scores:
        Microbe=line[0];Person=line[1];Bet=line[2:];Winner=Bet.index(max(Bet))

        Pareja_Tupla=(Microbe, Person)

        scores_dict[ Pareja_Tupla ]=Winner#Assign prediction
        Another_scores_dict[ Pareja_Tupla  ]=Bet
    #------------------------------------------
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
    #1.3.3. Build individual (per Fold) Confusion Matrices and Dictionary
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #1.3.3.1. Confusion Matrices
    #----------------------------------------------------------------------
    Confusion_Matrix_Patient_Fold=np.zeros((rows,cols))
    
    for line in Check:
        MicrobeCheck=float(line[0]); PersonCheck=float(line[1])
        Abundance=float(line[2])
        TuplaCheck=(MicrobeCheck,PersonCheck)
        
        Link_Counter[int(Abundance)]+=1#Abundance Counter

        # if Abundance==1:

        #     #Analizar Scores en detalle para link 1
        #     #-------------------------------------------------------
        #     Scores_One.append( Another_scores_dict[TuplaCheck] )#---
        #     #-------------------------------------------------------
            
        #     for rating in range(len(One_Ratings)):
        #         One_Ratings[rating]+=Another_scores_dict[TuplaCheck][rating]


        # elif Abundance==0:
            
        #     for rating in range(len(Zero_Ratings)):
        #         Zero_Ratings[rating]+=Another_scores_dict[TuplaCheck][rating]
                
        # elif Abundance==2:    
        #     for rating in range(len(Two_Ratings)):
        #         Two_Ratings[rating]+=Another_scores_dict[TuplaCheck][rating]

        # elif Abundance==3:    
        #     for rating in range(len(Three_Ratings)):
        #         Three_Ratings[rating]+=Another_scores_dict[TuplaCheck][rating]

        #Individual Confusion Matrix
        #---------------------------------------------------------------
        Confusion_Matrix_Patient_Fold\
            [scores_dict[TuplaCheck], int(Abundance)]+=1
        #---------------------------------------------------------------
        
    #----------------------------------------------------------------------

    #1.3.3.2. Fill dictionary
    #-----------------------------------------------------------------------
    try:
        Confusion_Matrices_dict[Patient][Fold]=Confusion_Matrix_Patient_Fold
    except KeyError:
        try:
            Confusion_Matrices_dict[Patient].append(Fold)
        except KeyError:
            Confusion_Matrices_dict[Patient]=\
            {Fold:Confusion_Matrix_Patient_Fold}
    #-----------------------------------------------------------------------
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#==============================================================================

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++


#3. CREATE CONFUSION MATRICES
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
print Confusion_Matrices_dict.keys(), len(Confusion_Matrices_dict.keys())

#3.1. Choose random Patients
#=======================================================================
Chosen_Patients= np.random.choice(Confusion_Matrices_dict.keys(), size=9)
print Chosen_Patients, len(Chosen_Patients)
#=======================================================================

#3.2. Generate Confusion Matrices
#===================================================================
Confusion_Matrices_Patients=[]

#3.2.1. By Patient
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
for chosen in Chosen_Patients:
    Confusion_Matrices_Patients.append\
        (Confusion_Matrix_func(Confusion_Matrices_dict,chosen))
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#3.2.2. Total 
#::::::::::::::::::::::::::::::::::::::::::::::::::::
Confusion_Matrix_Aggregated=\
Confusion_Matrix_func(Confusion_Matrices_dict,None)
print Confusion_Matrix_Aggregated
#::::::::::::::::::::::::::::::::::::::::::::::::::::
#===================================================================

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#3. PLOTS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#3.1. Gridspec of Confusion Matrices
#============================================================================

#3.1.1. Fontsizes and Parameters
#=============================
size_eje=12
size_ticks=13
size_title=15
size_letter=17
#=============================

Information=['Conf_Matrices','Single_Conf_Matrix','Heatmap_Scores']
Information_Requested=Information[0]

if Information_Requested=='Conf_Matrices':
    #3.1.1 Parameters for gridspec
    #=============================
    Width=Number_of_Links*2
    Height=Width*1.2
    fig=plt.figure(figsize=(Width,Height))
    gs = gridspec.GridSpec(4, 3)
    gs.update(\
    left=0.08,right=0.99,bottom=0.07,top=0.96,wspace=0.25,hspace=0.65) 
    #=============================
    
    #3.1.3. Fill Gridspec
    #=========================================================================
    #0.0. Aggregated Confusion Matrix
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ax_Tot = plt.subplot(gs[0, 0])
    Heatmap=Confusion_Matrix_Aggregated
    sns.heatmap(Confusion_Matrix_Aggregated,linewidths=0.2,\
                annot=True,cmap="RdBu_r",cbar=False,vmin=0,vmax=1)

    ax_Tot.set_title('Aggregated from 30 patients',fontsize=size_eje)
    ax_Tot.set_xlabel('Actual',fontsize=size_eje)
    ax_Tot.set_ylabel('Predicted',fontsize=size_eje)
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #(0,1). Write number of links found in sample
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    axText= plt.subplot(gs[0,1])
    sns.set_style("white")
    sns.despine
    axText.set_yticks([])
    axText.set_xticks([])

    #Loop to write individual link + total
    #--------------------------------------------------------------
    for link_txt in range(Number_of_Links):
        axText.text(0,0.9-link_txt*0.22,\
        'Link ' + str(link_txt)+ ': '+str(Link_Counter[link_txt]) )

    axText.text(0.5,0.05,'Total: '+str(sum(Link_Counter)) )
    #--------------------------------------------------------------
    
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #(1,0) to (2,2). Link Distribution in Sample
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    Counter_Plots=0
    for row_plot in range(1, 1+ len(Chosen_Patients)/3 ):
        for col_plot in range( len(Chosen_Patients)/3 ):
            ax = plt.subplot(gs[row_plot, col_plot])
            Heatmap=Confusion_Matrices_Patients[Counter_Plots]
            sns.heatmap(Heatmap,\
            cmap="RdBu_r", annot=True,cbar=False,vmin=0,vmax=1,linewidths=0.2)

            ax.set_title(\
            'Patient Number '+Chosen_Patients[Counter_Plots],fontsize=size_eje)
            ax.set_xlabel('Actual',fontsize=size_eje)

            if col_plot==0:            
                ax.set_ylabel('Predicted',fontsize=size_eje)

            Counter_Plots+=1
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    plt.savefig('../../Plots/First_Path/Confusion_Matrices_'+Dataset+'.pdf',\
                dpi=300)
    plt.show()
    #=========================================================================


elif Information_Requested=='Single_Conf_Matrix':
    
    #3.2. Confusion Matrix of single patient
    #=========================================================================
    fig=plt.figure(figsize=(7,6))
    
    Heatmap=Confusion_Matrix_Aggregated
    Plot=sns.heatmap(Confusion_Matrix_Aggregated,linewidths=0.2,\
                annot=True,cmap="RdBu_r",cbar=True,vmin=0,vmax=1)
    

    Plot.set_title(\
    'Aggregated ' +Dataset,fontsize=12)

    Plot.set_xlabel('Actual',fontsize=14)
    Plot.set_xticklabels(['Zero', 'Low','Medium', 'High'])
    Plot.set_ylabel('Predicted',fontsize=14)
    Plot.set_yticklabels( ['High', 'Medium', 'Low', 'Zero'] )

    plt.savefig('../../Plots/First_Path/Confusion_Matrix_'+Dataset+'.pdf',\
                dpi=300)
    plt.show()
    #========================================================================= 

elif Information_Requested=='Heatmap_Scores':    
    #3.3. Heatmap of score analysis
    #=========================================================================

    #3.3.1. Convert list of lists to matrix
    #:::::::::::::::::::::::::::::::::::::::::::::
    Matrix=np.asarray(Scores_One)
    print Matrix, len(Matrix), np.shape(Matrix)
    #:::::::::::::::::::::::::::::::::::::::::::::

    #3.3.2. Do clustering on the matrix to improve visualization
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #Clustering
    #--------------------------------------------------------
    Clustering=scipy.cluster.hierarchy.fclusterdata\
                                (Matrix,0.0,method='average')
    #--------------------------------------------------------

    #Dictionary cluster_number: node
    #--------------------------------------------------------
    Clustering_Dict={}
    for line in range(0,len(Clustering)):
        try:
            Clustering_Dict[ Clustering[line] ].append(line)
        except KeyError:
            Clustering_Dict[ Clustering[line] ]=[line]
    #--------------------------------------------------------

    #Rearrange matrix with Dictionary
    #----------------------------------
    Matrix_Clustered=np.zeros( Matrix.shape )

    counter=0;New_Order=[];
    Clusters=[]
    for key in Clustering_Dict:
        for node in Clustering_Dict[key]:

            Clusters.append(key)

            New_Order.append(node)

            Matrix_Clustered[counter]=Matrix[node]
            counter+=1
    #----------------------------------
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
    #3.3.3. Plot heatmap
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    fig=plt.figure(figsize=(7,12))
    Plot=sns.heatmap(Matrix_Clustered,cmap="Greys")
    Plot.set_title(Dataset)
    Plot.set_yticklabels( [] )
    plt.savefig('../../Plots/First_Path/Scores_On_One_'+Dataset+'.pdf',\
                dpi=300)
    plt.show()
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    #=========================================================================

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
