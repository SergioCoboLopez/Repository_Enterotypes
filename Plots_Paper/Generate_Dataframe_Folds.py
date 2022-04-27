#5/09/2019. Calculamos accuracy, precision y recall del algoritmo de clustering a partir de los scores generados

import sys
import os
import numpy as np
import pandas as pd
import scipy.cluster
from scipy import stats
import math
import pickle

import seaborn as sns
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt

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

#0.2. Precision
#===============================================
def Precision_fun(Matrix):
    Predicted=Matrix[1]
    Precision=Matrix[1][1]/(sum(Predicted) + 1e-8)
    return Precision
#===============================================

#0.3. Recall
#===============================================
def Recall_fun(Matrix):
    Actual=np.transpose(Matrix)[1]
    Recall=Matrix[1][1]/sum(Actual)
    return Recall
#===============================================

#0.4. Accuracy
#===============================================
def Accuracy_fun(Matrix):
    Accuracy=0
    for Link in range(Matrix.shape[1]):
        Accuracy+=Matrix[Link][Link]/\
        (sum(sum(Matrix)))
    return Accuracy
#===============================================

#0.5. Clean File                                             
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

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++


#1. DISTANCES INTRA/INTER
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#1.1. Load Data
#==================================================================
Distance_Data = pickle.load( open( "Pairs_of_Distances.p", "rb" ) )
#==================================================================


#1.2. Compute averages
#===================================================================
Datasets=['S-8_Tot','V-10_Tot','V-22_Tot','V-23_24_Tot','V-25_Tot' ]
Taxonomic_Orders=['Phylum','Class','Order','Family','Genus']
#===================================================================

#2. DIFFERENCES MMSBM/CLUSTER
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Number_of_Links=2

K=10;L=20

rows=Number_of_Links;cols=Number_of_Links
Link_Counter=[0]*Number_of_Links#->We count links present in the 30 pat. sample
Folds=5;Seeds=1

Patients_per_Dataset=\
{'Test':10, 'S-8_Tot':107, 'V-10_Tot':92, 'V-22_Tot':467, 'V-23_24_Tot':222,\
 'V-25_Tot':883}

Microbes_per_Dataset=\
{'Test':10, 'S-8_Tot':128, 'V-10_Tot':137, 'V-22_Tot':134, 'V-23_24_Tot':118,\
 'V-25_Tot':144}

Differences={'S-8_Tot':{}, 'V-10_Tot':{}, 'V-22_Tot':{}, 'V-23_24_Tot':{}, 'V-25_Tot':{}}

Datasets=['S-8_Tot', 'V-10_Tot', 'V-22_Tot', 'V-23_24_Tot', 'V-25_Tot']
Datasets_label=['S-8', 'V-10', 'V-22', 'V-23_24', 'V-25']

Confusion_Matrices_per_Patient={'S-8_Tot':{}, 'V-10_Tot':{}, 'V-22_Tot':{}, 'V-23_24_Tot':{}, 'V-25_Tot':{} }
Patient_wise_Metrics_Clustering={'S-8_Tot':{}, 'V-10_Tot':{}, 'V-22_Tot':{}, 'V-23_24_Tot':{}, 'V-25_Tot':{} }

Confusion_Matrices_per_Patient_MMSBM={'S-8_Tot':{},'V-10_Tot':{},'V-22_Tot':{},'V-23_24_Tot':{},'V-25_Tot':{} }
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Confusion_Matrices_per_Patient_MMSBM_Accuracy=\
{'S-8_Tot':{},'V-10_Tot':{},'V-22_Tot':{},'V-23_24_Tot':{},'V-25_Tot':{} }
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Patient_wise_Metrics_MMSBM={'S-8_Tot':{}, 'V-10_Tot':{}, 'V-22_Tot':{}, 'V-23_24_Tot':{}, 'V-25_Tot':{} }
Patient_wise_Metrics_Cluster={'S-8_Tot':{}, 'V-10_Tot':{}, 'V-22_Tot':{}, 'V-23_24_Tot':{}, 'V-25_Tot':{} }

#BUILD CONFUSION MATRICES
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#2.1. Paths Clustering
#====================================================================
path_common_Clustering='../../RClustering/'
#====================================================================

#2.2. Diccionario de confusion matrices de paciente -> Seed -> Fold
#====================================================================
for Dataset in Datasets:

    Patients=Patients_per_Dataset[Dataset]
    Microbes=Microbes_per_Dataset[Dataset]

    path_I_Clustering='../../Input_Data/Leave_One_Out_' + Dataset + '/'
    path_O_Clustering=path_common_Clustering+'robust-clustering-metagenomics/Comparison/Output_Data/'+Dataset+'/'

    Confusion_Matrices_dict={}

    for Patient_Id in range(Patients):
        for Fold_Id in range(Folds):
            for Seed_Id in range(1,Seeds+1):

                #2.2.1. Correr sobre scores
                #::::::::::::::::::::::::::::::::::::::::::::::::::::
                Name='P_' + str(Patient_Id) + '_F_' + str(Fold_Id)
                Full_Name_Seed=Name + '_Seed_' + str(Seed_Id)

                TestCheck=read_file(path_I_Clustering + Name + '_TestCheck.txt')
                Scores=read_file(path_O_Clustering + Full_Name_Seed + '_Scores.txt')

                scores_dict={}
                for line in Scores:
                    line=[float(element) for element in line]
                    Microbe=float(line[0]);Patient=float(line[1]);Bet=line[2:]
                    Pareja_Tupla=(Microbe, Patient)

                    if Number_of_Links==3:
                        #............................................
                        Winner=Bet.index(max(Bet))
                        scores_dict[ Pareja_Tupla ]=Winner#Assign prediction
                        #............................................

                    #Reducir a dos scores para precision y recall
                    elif Number_of_Links==2:
                        #...........................................
                        New_line=[Microbe,Patient,Bet[0],\
                                  Bet[1]+Bet[2]]
                        Winner=New_line[2:].index(max(New_line[2:]))
                        scores_dict[ Pareja_Tupla ]=Winner
                        #...........................................

                #:::::::::::::::::::::::::::::::::::::::::::::::::::

                #2.2.2. Crear matrices de confusion
                #:::::::::::::::::::::::::::::::::::::::::::::::::::
                Confusion_Matrix_Patient_Fold_Seed=np.zeros((rows,cols))

                for line in TestCheck:
                    MicrobeCheck=float(line[0]); PersonCheck=float(line[1])
                    Abundance=float(line[2])
                    TuplaCheck=(MicrobeCheck,PersonCheck)

                    #DOS LINKS
                    #...............................................
                    if Number_of_Links==2:
                        if Abundance==2:
                            Abundance=1.0
                    #...............................................

                    Link_Counter[int(Abundance)]+=1#Abundance Counter

                    #Create individual Confusion Matrix
                    #-----------------------------------------------
                    Confusion_Matrix_Patient_Fold_Seed\
                        [scores_dict[TuplaCheck], int(Abundance)]+=1
                    #-----------------------------------------------

                #1.3.3.2. Fill dictionary
                #---------------------------------------------------
                try:
                    Confusion_Matrices_dict[Patient_Id][Seed_Id]\
                        [Fold_Id]=Confusion_Matrix_Patient_Fold_Seed
                except KeyError:
                    try:
                        Confusion_Matrices_dict[Patient_Id][Seed_Id].append(Fold_Id)
                    except KeyError:
                        try:
                            Confusion_Matrices_dict[Patient_Id]\
                            [Seed_Id]={Fold_Id:Confusion_Matrix_Patient_Fold_Seed}
                        except KeyError:
                            try:
                                Confusion_Matrices_dict[Patient_Id].append(Seed_Id)
                            except KeyError:
                                Confusion_Matrices_dict[Patient_Id]=\
                                {Seed_Id:{Fold_Id:Confusion_Matrix_Patient_Fold_Seed}}
                #---------------------------------------------------
                #:::::::::::::::::::::::::::::::::::::::::::::::::::
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    #3. MMSBMS
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    #3.1. Paths MMSBM
    #================================================================
    path_Out='../../Output_Data/'+ 'Leave_One_Out_' + Dataset + '/'
    path_In= '../../Input_Data/' + 'Leave_One_Out_' + Dataset + '/'
    #================================================================

    #3.2. Choose best seeds
    #================================================================

    #Read best results from best likelihoods
    #--------------------------------------------                    
    Results=read_file(path_In + 'Best_Seeds.txt')
    #--------------------------------------------

    Confusion_Matrices_MMSBM_dict={}    
    Confusion_Matrices_MMSBM_dict_Accuracy={}    

    for line_file in Results:

        #2.2.1. Read Input files
        #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        Patient_R=line_file[0];Fold=line_file[1];Seed=line_file[2]
        Parameters = str(Patient_R) + "_Fold_" + str(Fold) +  "_Seed_" + str(Seed)

        Scores = path_Out + Parameters +  '_scores.txt'
        raw_scores=read_file(Scores)
        Check=read_file(path_In + 'P_' + Patient_R + '_F_'+Fold+ '_TestCheck.txt')
        #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


        #2.2.2. Read and compute MMSBM scores
        #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        #2.2.2.1. Clean scores file
        #----------------------------
        scores=clean_data(raw_scores)
        #----------------------------

        #2.2.2.2. Collect scores in dictionary
        #-----------------------------------------------------------
        #~~~~~~~~~~~~~~~~~~~~~~~~~
        scores_dict_Accuracy={} #~  
        #~~~~~~~~~~~~~~~~~~~~~~~~~
        scores_dict={}
        for line in scores:

            Microbe_R=line[0];Person=line[1];Bet=line[2:]
            Max_score=max(Bet);Pareja_Tupla=(Microbe_R, Person)

            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            #TRES LINKS: Accuracy
            #.......................................................
            Winner_Accuracy=Bet.index(max(Bet))
            #Assign prediction
            scores_dict_Accuracy[ Pareja_Tupla ]=Winner_Accuracy
            #.......................................................

            #DOS LINKS: Precision y Recall
            #.......................................................
            New_line=[ Microbe_R, Person, Bet[0], Bet[1]+Bet[2] ]
            Winner=New_line[2:].index(max(New_line[2:]))
            scores_dict[ Pareja_Tupla ]=Winner
            #.......................................................
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        Matrix_Counter_perFold=np.zeros((rows,cols))

        Matrix_Counter_perFold_Accuracy=np.zeros((3,3))

        for line in Check:

            MicrobeCheck=float(line[0]); PersonCheck=float(line[1])
            Abundance=float(line[2])
            TuplaCheck=(MicrobeCheck,PersonCheck)

            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            #Individual Confusion Matrix 3x3, Accuracy
            #-------------------------------------------------------
            Matrix_Counter_perFold_Accuracy[scores_dict_Accuracy[TuplaCheck],int(Abundance)]+=1
            #-------------------------------------------------------
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            #DOS LINKS 
            #.......................................................
            if Number_of_Links==2:
                if Abundance==2:
                    Abundance=1.0
            #.......................................................

            #Individual Confusion Matrix 
            #-------------------------------------------------------
            Matrix_Counter_perFold[scores_dict[TuplaCheck],int(Abundance)]+=1
            #-------------------------------------------------------

        #2.2.3.3. Fill dictionaries
        #-----------------------------------------------------------
        try:
            Confusion_Matrices_MMSBM_dict[Patient_R][Fold]=Matrix_Counter_perFold
        except KeyError:

            try:
                Confusion_Matrices_MMSBM_dict[Patient_R].append(Fold)
            except KeyError:
                Confusion_Matrices_MMSBM_dict[Patient_R]=\
                {Fold:Matrix_Counter_perFold}

        #Dictionaries 3x3, Accuracy
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        try:
            Confusion_Matrices_MMSBM_dict_Accuracy[Patient_R][Fold]=Matrix_Counter_perFold_Accuracy
        except KeyError:

            try:
                Confusion_Matrices_MMSBM_dict_Accuracy[Patient_R].append(Fold)
            except KeyError:
                Confusion_Matrices_MMSBM_dict_Accuracy[Patient_R]=\
                {Fold:Matrix_Counter_perFold_Accuracy}
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        #-----------------------------------------------------------

        #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#BUILD DICTIONARIES OF METRICS
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    for Host_R in Confusion_Matrices_MMSBM_dict:
        for Fold in Confusion_Matrices_MMSBM_dict[Host_R]:
            #Compute metrics
            #---------------------------------------------------
            Accuracy_Cluster=Accuracy_fun(Confusion_Matrices_dict[int(Host_R)][1][int(Fold)])
            Accuracy_MMSBM=Accuracy_fun(Confusion_Matrices_MMSBM_dict[Host_R][Fold])
            #---------------------------------------------------

            #---------------------------------------------------
            Precision_Cluster=Precision_fun(Confusion_Matrices_dict[int(Host_R)][1][int(Fold)])
            Precision_MMSBM=Precision_fun(Confusion_Matrices_MMSBM_dict[Host_R][Fold])
            #---------------------------------------------------

            #---------------------------------------------------
            Recall_Cluster=Recall_fun(Confusion_Matrices_dict[int(Host_R)][1][int(Fold)])
            Recall_MMSBM=Recall_fun(Confusion_Matrices_MMSBM_dict[Host_R][Fold])
            #---------------------------------------------------

            #Save to dictionaries 
            #---------------------------------------------------
            try:
                Patient_wise_Metrics_MMSBM[Dataset][Host_R][Fold]=\
                {"Accuracy":Accuracy_MMSBM, "Precision":Precision_MMSBM, "Recall":Recall_MMSBM}

                Patient_wise_Metrics_Cluster[Dataset][Host_R][Fold]=\
                {"Accuracy":Accuracy_Cluster, "Precision":Precision_Cluster,"Recall":Recall_Cluster}
                
            except KeyError:
                Patient_wise_Metrics_MMSBM[Dataset][Host_R]=\
                {Fold:{"Accuracy":Accuracy_MMSBM,"Precision":Precision_MMSBM,"Recall":Recall_MMSBM}}
                Patient_wise_Metrics_Cluster[Dataset][Host_R]=\
                {Fold:{"Accuracy":Accuracy_Cluster, "Precision":Precision_Cluster, "Recall":Recall_Cluster}}
            #---------------------------------------------------
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#SAVE DICTIONARIES TO DATAFRAMES
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Dict_Dataframe={'Performance':[],'Patient':[],'Fold':[],'Dataset':[],'Metric':[],'Method':[]}

#MMSBM Dictionary                                                                                                          
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::              
for Dataset in Patient_wise_Metrics_MMSBM:  
    for Patient in Patient_wise_Metrics_MMSBM[Dataset]:
        for Fold in Patient_wise_Metrics_MMSBM[Dataset][Patient]:
            for metric in Patient_wise_Metrics_MMSBM[Dataset][Patient][Fold]:
                
                Dict_Dataframe['Performance'].append(Patient_wise_Metrics_MMSBM[Dataset][Patient][Fold][metric])
                Dict_Dataframe['Fold'].append(Fold)
                Dict_Dataframe['Patient'].append(Patient)
                #Recortar el "Tot" del dataset
                #---------------------------------------------------
                Dataset_Formal=Dataset.split("_")[0]
                if len(Dataset.split("_"))>2:
                    Dataset_Formal="_".join(Dataset.split("_")[0:2])
                #---------------------------------------------------
                Dict_Dataframe['Dataset'].append(Dataset_Formal)
                Dict_Dataframe['Metric'].append(metric)
                Dict_Dataframe['Method'].append('MMSBM')
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


#Clustering Dictionary
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::              
for Dataset in Patient_wise_Metrics_Cluster:  
    for Patient in Patient_wise_Metrics_Cluster[Dataset]:
        for Fold in Patient_wise_Metrics_Cluster[Dataset][Patient]:
            for metric in Patient_wise_Metrics_Cluster[Dataset][Patient][Fold]:
                
                Dict_Dataframe['Performance'].append(Patient_wise_Metrics_Cluster[Dataset][Patient][Fold][metric])
                Dict_Dataframe['Fold'].append(Fold)
                Dict_Dataframe['Patient'].append(Patient)
                #Recortar el "Tot" del dataset
                #---------------------------------------------------
                Dataset_Formal=Dataset.split("_")[0]
                if len(Dataset.split("_"))>2:
                    Dataset_Formal="_".join(Dataset.split("_")[0:2])
                #---------------------------------------------------
                Dict_Dataframe['Dataset'].append(Dataset_Formal)
                Dict_Dataframe['Metric'].append(metric)
                Dict_Dataframe['Method'].append('Cluster')
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


Dataframe_Test=pd.DataFrame(Dict_Dataframe)                                                                                               
Dataframe_Test.to_csv("Metrics_by_Folds.csv")
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
