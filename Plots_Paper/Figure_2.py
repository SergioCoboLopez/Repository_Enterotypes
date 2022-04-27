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

#CALCULATIONS
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Distances={}                                                     
for Dataset in Distance_Data:
    for taxonomy in Distance_Data[Dataset]:                      
            
        Intra=Distance_Data[Dataset][taxonomy]['Intra']          

        Inter=Distance_Data[Dataset][taxonomy]['Inter']          

        #ERROR                                                   
        #--------------------------------------------------------
        #We do a propagation of errors to compute the actual error:
        #sqrt[ SEM(Intra)^2/Intra^2 + SEM(Inter)^2/Inter^2 ]

        Inter_Mean=np.mean(Inter)
        SEM_Inter=scipy.stats.sem(Inter,ddof=0)
        
        Intra_Mean=np.mean(Intra)
        SEM_Intra=scipy.stats.sem(Intra,ddof=0)

        #Error                                                    
        #........................................................ 
        Error=math.sqrt( \
        (SEM_Intra/Intra_Mean)**2 + (SEM_Inter/Inter_Mean)**2 )
        #.......................................................

        #-------------------------------------------------------
        try:                                                     
            Distances[Dataset][taxonomy]={"Intra_Mean":Intra_Mean,"Inter_Mean":Inter_Mean,"Error":Error }                             
        except KeyError:    
            Distances[Dataset]={taxonomy:{"Intra_Mean":Intra_Mean,"Inter_Mean":Inter_Mean,"Error":Error } }
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#2. DIFFERENCES MMSBM/CLUSTER
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Number_of_Links=3
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Differences_test={'S-8_Tot':{}, 'V-10_Tot':{}, 'V-22_Tot':{}, 'V-23_24_Tot':{}, 'V-25_Tot':{}}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

#2. CLUSTERING
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
#===================================================================


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

#4.Difference between patients averaging by folds with the 2 methods
#===================================================================

    for Host_R in Confusion_Matrices_MMSBM_dict:
        Differences_Prec=[];Differences_Recall=[];Differences_Acc=[]
        Relative_Differences_Prec=[];Relative_Differences_Recall=[];
        Relative_Differences_Acc=[]

        Accuracy_Cluster_List=[]
        Accuracy_MMSBM_List=[]
        
        Precision_Cluster_List=[]
        Precision_MMSBM_List=[]
        
        Recall_Cluster_List=[]
        Recall_MMSBM_List=[]


        for Fold in Confusion_Matrices_MMSBM_dict[Host_R]:

            #---------------------------------------------------
            Accuracy_Cluster=Accuracy_fun(Confusion_Matrices_dict[int(Host_R)][1][int(Fold)])
            Accuracy_MMSBM=Accuracy_fun(Confusion_Matrices_MMSBM_dict[Host_R][Fold])
                
            Difference_Accuracy=Accuracy_MMSBM - Accuracy_Cluster
            Differences_Acc.append(Difference_Accuracy)

            Relative_Difference_Accuracy=\
            Difference_Accuracy/Accuracy_Cluster
            Relative_Differences_Acc.append(\
            Relative_Difference_Accuracy)

            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            Accuracy_Cluster_List.append(Accuracy_Cluster)
            Accuracy_MMSBM_List.append(Accuracy_MMSBM)
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
            #---------------------------------------------------

            #---------------------------------------------------
            Precision_Cluster=Precision_fun(\
            Confusion_Matrices_dict[int(Host_R)][1][int(Fold)])
            Precision_MMSBM=Precision_fun(\
            Confusion_Matrices_MMSBM_dict[Host_R][Fold])

            Difference_Precision=Precision_MMSBM - Precision_Cluster
            Differences_Prec.append(Difference_Precision)

            Relative_Difference_Precision=\
            Difference_Precision/Precision_Cluster
            
            if math.isnan(Relative_Difference_Precision)==True:
                    Relative_Difference_Precision=0

            if math.isinf(Relative_Difference_Precision)==True:
                Relative_Difference_Precision=1
                    
            Relative_Differences_Prec.append(\
            Relative_Difference_Precision)

            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            Precision_Cluster_List.append(Precision_Cluster)
            Precision_MMSBM_List.append(Precision_MMSBM)
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
            #---------------------------------------------------

            #---------------------------------------------------
            Recall_Cluster=Recall_fun(\
            Confusion_Matrices_dict[int(Host_R)][1][int(Fold)])
            Recall_MMSBM=Recall_fun(\
            Confusion_Matrices_MMSBM_dict[Host_R][Fold])

            Difference_Recall=Recall_MMSBM - Recall_Cluster
            Differences_Recall.append(Difference_Recall)

            Relative_Difference_Recall=\
            Difference_Recall/Recall_Cluster
            
            if math.isnan(Relative_Difference_Recall)==True:
                Relative_Difference_Recall=0
            if math.isinf(Relative_Difference_Recall)==True:
                Relative_Difference_Recall=1
            Relative_Differences_Recall.append(\
                    Relative_Difference_Recall)

            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            Recall_Cluster_List.append(Recall_Cluster)
            Recall_MMSBM_List.append(Recall_MMSBM)
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
            #---------------------------------------------------

        #-------------------------------------------------------
        Mean_Difference_Accuracy=np.mean(Differences_Acc)
        Mean_Difference_Prec=np.mean(Differences_Prec)
        Mean_Difference_Recall=np.mean(Differences_Recall)

        Mean_Relative_Difference_Accuracy=\
        np.mean(Relative_Differences_Acc)
        Mean_Relative_Difference_Prec=\
        np.mean(Relative_Differences_Prec)
        Mean_Relative_Difference_Recall=\
        np.mean(Relative_Differences_Recall)

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Mean_Accuracy_Cluster=np.mean(Accuracy_Cluster_List)
        Mean_Accuracy_MMSBM=np.mean(Accuracy_MMSBM_List)
        print(Dataset)
        print(Host_R)
        print(Mean_Accuracy_MMSBM)

        Mean_Precision_Cluster=np.mean(Precision_Cluster_List)
        Mean_Precision_MMSBM=np.mean(Precision_MMSBM_List)

        Mean_Recall_Cluster=np.mean(Recall_Cluster_List)
        Mean_Recall_MMSBM=np.mean(Recall_MMSBM_List)
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        #-------------------------------------------------------

        #Relative improvement
        #.......................................................
        Differences_test[Dataset][Host_R]={\
        'Accuracy': Mean_Relative_Difference_Accuracy,\
        'Precision': Mean_Relative_Difference_Prec,\
        'Recall': Mean_Relative_Difference_Recall}

        #Media del cociente
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Differences_test[Dataset][Host_R]={\
        # 'Accuracy': Mean_Accuracy_MMSBM/Mean_Accuracy_Cluster,\
        # 'Precision':Mean_Precision_MMSBM/Mean_Precision_Cluster,\
        # 'Recall': Mean_Recall_MMSBM/Mean_Recall_Cluster}

        # Differences_test[Dataset][Host_R]={\
        # 'Accuracy': Mean_Accuracy_Cluster/Mean_Accuracy_MMSBM,\
        # 'Precision':Mean_Precision_Cluster/Mean_Precision_MMSBM,\
        # 'Recall': Mean_Recall_Cluster/Mean_Recall_MMSBM}
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        #Cociente de las medias
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Differences[Dataset][Host_R]={\
        'Accuracy': (Mean_Accuracy_MMSBM, Mean_Accuracy_Cluster),\
        'Precision':(Mean_Precision_MMSBM, Mean_Precision_Cluster),\
        'Recall': (Mean_Recall_MMSBM, Mean_Recall_Cluster)}
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #.......................................................

        #-------------------------------------------------------

        #=======================================================
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#4. DATAFRAMES
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#4.1. Relative Distances
#===================================================================
    
#Cociente de medias
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Mean_Accuracies=[];Mean_Precisions=[];Mean_Recalls=[]
SEM_Accuracies=[];SEM_Precisions=[];SEM_Recalls=[]

Mean_Relative_Accuracies=[];Mean_Relative_Precisions=[]
Mean_Relative_Recalls=[]

SEM_Relative_Accuracies=[];SEM_Relative_Precisions=[]
SEM_Relative_Recalls=[]

for Dataset in Datasets:

    Accuracies_MMSBM=[];Accuracies_Cluster=[]
    Precisions_MMSBM=[];Precisions_Cluster=[]
    Recalls_MMSBM=[];Recalls_Cluster=[]
    
    for Patient in Differences[Dataset]:
        
        Accuracies_MMSBM.append(Differences[Dataset][Patient]['Accuracy'][0])
        Accuracies_Cluster.append(Differences[Dataset][Patient]['Accuracy'][1])

        Precisions_MMSBM.append(Differences[Dataset][Patient]['Precision'][0])
        Precisions_Cluster.append(Differences[Dataset][Patient]['Precision'][1])

        Recalls_MMSBM.append(Differences[Dataset][Patient]['Recall'][0])
        Recalls_Cluster.append(Differences[Dataset][Patient]['Recall'][1])
        
    #------------------------------------------------------------
    Mean_Accuracies.append(math.log(np.mean(Accuracies_MMSBM)/np.mean(Accuracies_Cluster)))
    SEM_Accuracies.append(scipy.stats.sem(Accuracies_MMSBM,ddof=0)/scipy.stats.sem(Accuracies_Cluster,ddof=0))
    
    Mean_Precisions.append(math.log(np.mean(Precisions_MMSBM)/np.mean(Precisions_Cluster)))
    SEM_Precisions.append(scipy.stats.sem(Precisions_MMSBM,ddof=0)/scipy.stats.sem(Precisions_Cluster,ddof=0))

    Mean_Recalls.append(math.log(np.mean(Recalls_MMSBM)/np.mean(Recalls_Cluster)))
    SEM_Recalls.append(scipy.stats.sem(Recalls_MMSBM,ddof=0)/scipy.stats.sem(Recalls_Cluster,ddof=0))
    #-----------------------------------------------------------

print(Mean_Accuracies)
print(Mean_Recalls)
print(Mean_Precisions)
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


#PLOTS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
fig=plt.figure(figsize=(12.5,4.5))
gs=gridspec.GridSpec(1, 2)
gs.update(left=0.07,right=0.99,bottom=0.1,top=0.955,wspace=0.23,\
          hspace=0.25)

#(0,0). Plot
#====================================================================
ax_00=plt.subplot(gs[0,0])

xticks_pos=[tick for tick in Datasets_label]

plt.errorbar(xticks_pos, Mean_Accuracies, yerr=0, marker='o', label='Accuracy')
plt.errorbar(xticks_pos, Mean_Precisions, yerr=0, marker='o', label='Precision')
plt.errorbar(xticks_pos, Mean_Recalls, yerr=0, marker='o', label='Recall')

plt.xlabel('Dataset', fontsize= 12)
plt.ylabel(r'$\log \frac{\langle Latent \rangle}{\langle Conventional \rangle}$',fontsize=12 )

plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.legend()
#====================================================================


#(0,1). Plot
#====================================================================
ax_01=plt.subplot(gs[0,1])

xticks_pos=[tick for tick in Taxonomic_Orders  ]                  
#Plot                                                            
#---------------------------------------------------------------
for Data_plot in Datasets:

    Means=[];Error=[]                                            
    for clade in Taxonomic_Orders:                               
            
        Means.append(math.log(Distances[Data_plot][clade]["Intra_Mean"]/Distances[Data_plot][clade]["Inter_Mean"],math.e) )

        Error.append( (Distances[Data_plot][clade]["Error"]) )
            

    plt.errorbar(xticks_pos, Means, yerr=Error, marker='o')      
    #--------------------------------------------------------------- 
Datasets_Names=['S-8','V-10','V-22','V-23_24','V-25' ]
    
#Cosmetics                                                       
#---------------------------------------------------------------
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)

plt.legend(Datasets_Names)                                       
plt.ylabel(r'log( $\langle$ $d_{\it{Intra}}$ $\rangle$/$\langle d_{\it{Inter}} \rangle$ )',fontsize=12)
plt.xlabel('Clades',fontsize=12)
#---------------------------------------------------------------
#====================================================================


path_plot='../../Plots_Paper/'
plt.savefig(path_plot + 'Figure_2.pdf',dpi=300)

plt.show()
#====================================================================
