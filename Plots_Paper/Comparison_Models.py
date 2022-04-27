#5/09/2019. Calculamos accuracy, precision y recall del algoritmo de clustering a partir de los scores generados

import sys
import os
import numpy as np
import pandas as pd
import scipy.cluster
from scipy import stats
import math

import seaborn as sns
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


#1. PARAMETERS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#NanSaver
#NanSaver=1e-8
path_plot='../../Plots_Paper/'

#Choose measure
Measures=['Absolute', 'Relative']
Measure='Relative'


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
                    Confusion_Matrices_dict[Patient_Id][Seed_Id][Fold_Id]=Confusion_Matrix_Patient_Fold_Seed
                except KeyError:
                    try:
                        Confusion_Matrices_dict[Patient_Id][Seed_Id].append(Fold_Id)
                    except KeyError:
                        try:
                            Confusion_Matrices_dict[Patient_Id][Seed_Id]=\
                            {Fold_Id:Confusion_Matrix_Patient_Fold_Seed}
                        except KeyError:
                            try:
                                Confusion_Matrices_dict[Patient_Id].append(Seed_Id)
                            except KeyError:
                                Confusion_Matrices_dict[Patient_Id]=\
                                {Seed_Id:{Fold_Id:Confusion_Matrix_Patient_Fold_Seed}}
                #---------------------------------------------------
                #:::::::::::::::::::::::::::::::::::::::::::::::::::
#===================================================================

    #===============================================================
    if Measure=='Absolute':
        #2.3. Aggregate Confusion matrices by Folds
        #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        for Patient_Id in Confusion_Matrices_dict:
            for Seed_Id in Confusion_Matrices_dict[Seed_Id]:
                Confusion_raw=sum(Confusion_Matrices_dict[Patient_Id][Seed_Id].values())
                try:
                    Confusion_Matrices_per_Patient[Dataset][Patient_Id][Seed_Id]=Confusion_raw
                except KeyError:
                    Confusion_Matrices_per_Patient[Dataset][Patient_Id]=\
                    {Seed_Id: Confusion_raw}
        #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        #2.4. Compute Patient-wise metrics
        #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        for Patient_Id in Confusion_Matrices_per_Patient[Dataset]:   
            for Seed_Id in Confusion_Matrices_per_Patient[Dataset][Patient_Id]:

                #2.4.5. Precision, Recall, Accuracy
                #---------------------------------------------------
                Precision=Precision_fun(Confusion_Matrices_per_Patient[Dataset][Patient_Id][Seed_Id])
                Recall=Recall_fun(Confusion_Matrices_per_Patient[Dataset][Patient_Id][Seed_Id])
                Accuracy=Accuracy_fun(Confusion_Matrices_per_Patient[Dataset][Patient_Id][Seed_Id])
                #---------------------------------------------------

                Patient_wise_Metrics_Clustering[Dataset][Patient_Id]=\
                                        {"Accuracy":Accuracy,"Precision":Precision,"Recall":Recall}
        #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        
    #===============================================================

    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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

    #===============================================================
    if Measure=='Absolute':
        #2.3. Aggregate Confusion matrices by Folds
        #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        for Patient_Id in Confusion_Matrices_MMSBM_dict:
                Confusion_MMSBM=sum(Confusion_Matrices_MMSBM_dict[Patient_Id].values())
                Confusion_Matrices_per_Patient_MMSBM[Dataset][Patient_Id]=Confusion_MMSBM
        #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        #2.3. Aggregate Confusion matrices 3x3, Accuracy by Folds
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        for Patient_Id in Confusion_Matrices_MMSBM_dict_Accuracy:
                Confusion_MMSBM_Accuracy=sum(Confusion_Matrices_MMSBM_dict_Accuracy[Patient_Id].values())
                Confusion_Matrices_per_Patient_MMSBM_Accuracy[Dataset][Patient_Id]=Confusion_MMSBM_Accuracy
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        #2.4. Compute Patient-wise metrics
        #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        for Patient_Id in Confusion_Matrices_per_Patient_MMSBM[Dataset]:   

                #2.4.5. Precision, Recall, Accuracy
                #---------------------------------------------------
                Precision=Precision_fun(Confusion_Matrices_per_Patient_MMSBM[Dataset][Patient_Id])
                Recall=Recall_fun(Confusion_Matrices_per_Patient_MMSBM[Dataset][Patient_Id])
                Accuracy=Accuracy_fun(Confusion_Matrices_per_Patient_MMSBM_Accuracy[Dataset][Patient_Id])
                #---------------------------------------------------

                Patient_wise_Metrics_MMSBM[Dataset][Patient_Id]\
                    ={"Accuracy":Accuracy, "Precision":Precision, "Recall":Recall}
        #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #===============================================================

    #4. Difference between patients after averaging by folds with the two methods
    #===============================================================
    if Measure=='Relative':
        for Patient_R in Confusion_Matrices_MMSBM_dict:
            Differences_Prec=[]; Differences_Recall=[]; Differences_Acc=[]
            Relative_Differences_Prec=[];Relative_Differences_Recall=[];Relative_Differences_Acc=[]

            for Fold in Confusion_Matrices_MMSBM_dict[Patient_R]:

                #---------------------------------------------------
                Accuracy_Cluster=Accuracy_fun(Confusion_Matrices_dict[int(Patient_R)][1][int(Fold)])
                Accuracy_MMSBM=Accuracy_fun(Confusion_Matrices_MMSBM_dict[Patient_R][Fold])
                
                Difference_Accuracy=Accuracy_MMSBM - Accuracy_Cluster
                Differences_Acc.append(Difference_Accuracy)

                Relative_Difference_Accuracy=Difference_Accuracy/Accuracy_Cluster
                Relative_Differences_Acc.append(Relative_Difference_Accuracy)
                #---------------------------------------------------

                #---------------------------------------------------
                Precision_Cluster=Precision_fun(Confusion_Matrices_dict[int(Patient_R)][1][int(Fold)])
                Precision_MMSBM=Precision_fun(Confusion_Matrices_MMSBM_dict[Patient_R][Fold])

                Difference_Precision=Precision_MMSBM - Precision_Cluster
                Differences_Prec.append(Difference_Precision)

                Relative_Difference_Precision=Difference_Precision/Precision_Cluster
                if math.isnan(Relative_Difference_Precision)==True:
                    Relative_Difference_Precision=0

                if math.isinf(Relative_Difference_Precision)==True:
                    Relative_Difference_Precision=1
                    
                Relative_Differences_Prec.append(Relative_Difference_Precision)
                #---------------------------------------------------

                #---------------------------------------------------
                Recall_Cluster=Recall_fun(Confusion_Matrices_dict[int(Patient_R)][1][int(Fold)])
                Recall_MMSBM=Recall_fun(Confusion_Matrices_MMSBM_dict[Patient_R][Fold])

                Difference_Recall=Recall_MMSBM - Recall_Cluster
                Differences_Recall.append(Difference_Recall)

                Relative_Difference_Recall=Difference_Recall/Recall_Cluster
                if math.isnan(Relative_Difference_Recall)==True:
                    Relative_Difference_Recall=0
                if math.isinf(Relative_Difference_Recall)==True:
                    Relative_Difference_Recall=1
                Relative_Differences_Recall.append(Relative_Difference_Recall)
                #---------------------------------------------------

            #-------------------------------------------------------
            Mean_Difference_Accuracy=np.mean(Differences_Acc)
            Mean_Difference_Prec=np.mean(Differences_Prec)
            Mean_Difference_Recall=np.mean(Differences_Recall)

            Mean_Relative_Difference_Accuracy=np.mean(Relative_Differences_Acc)
            Mean_Relative_Difference_Prec=np.mean(Relative_Differences_Prec)
            Mean_Relative_Difference_Recall=np.mean(Relative_Differences_Recall)
            #-------------------------------------------------------

            #-------------------------------------------------------

            #Relative improvement
            #.......................................................
            Differences[Dataset][Patient_R]={\
            'Accuracy': Mean_Relative_Difference_Accuracy,\
            'Precision': Mean_Relative_Difference_Prec, \
            'Recall': Mean_Relative_Difference_Recall}
            #.......................................................

            #-------------------------------------------------------

            #=======================================================

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
print(Differences)
#4. DATAFRAMES AND PLOTS
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#4.2. Absolute Distances
#===================================================================
if Measure=='Absolute':

    Dict_Dataframe={'Performance':[],'Patient':[],\
                    'Dataset':[],'Metric':[],'Method':[]}

    #MMSBM
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    for Dataset in Patient_wise_Metrics_MMSBM:
        for Patient in Patient_wise_Metrics_MMSBM[Dataset]:
             for metric in Patient_wise_Metrics_MMSBM[Dataset][Patient]:
                 Dict_Dataframe['Performance'].append(Patient_wise_Metrics_MMSBM[Dataset][Patient][metric])
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

    #CLUSTERING
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    for Dataset in Patient_wise_Metrics_Clustering:
        for Patient in Patient_wise_Metrics_Clustering[Dataset]:
             for metric in Patient_wise_Metrics_Clustering[Dataset][Patient]:
                 Dict_Dataframe['Performance'].append(Patient_wise_Metrics_Clustering[Dataset][Patient][metric])
                 Dict_Dataframe['Patient'].append(Patient)
                 #Recortar el "Tot" del dataset
                 #---------------------------------------------------
                 Dataset_Formal=Dataset.split("_")[0]
                 if len(Dataset.split("_"))>2:
                     Dataset_Formal="_".join(Dataset.split("_")[0:2])
                 #---------------------------------------------------
                 Dict_Dataframe['Dataset'].append(Dataset_Formal)
                 Dict_Dataframe['Metric'].append(metric)
                 Dict_Dataframe['Method'].append('Clustering')
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    Dataframe_Test=pd.DataFrame(Dict_Dataframe)
    Dataframe_Test.to_csv("Comparison_MMSBM_RClustering.csv")
#===================================================================


    g=sns.catplot(x='Metric', y='Performance',col='Dataset',hue='Method', col_order=Datasets_label,data=Dataframe_Test, kind='box',height=4.5, aspect=0.9)
    g.set_axis_labels("", "Performance")
    g.set_titles("{col_name}")

    #Boxplot
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    plt.text(-13.75,1.06, 'a', fontsize=12,fontweight='bold')
    plt.text(-10.5, 1.06, 'b', fontsize=12,fontweight='bold')
    plt.text(-7.25, 1.06, 'c', fontsize=12,fontweight='bold')
    plt.text(-4,    1.06, 'd', fontsize=12,fontweight='bold')
    plt.text(-0.75, 1.06, 'e', fontsize=12,fontweight='bold')
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #Violinplot
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    plt.text(-13.75,1.15, 'a', fontsize=12,fontweight='bold')
    plt.text(-10.5, 1.15, 'b', fontsize=12,fontweight='bold')
    plt.text(-7.25, 1.15, 'c', fontsize=12,fontweight='bold')
    plt.text(-4,    1.15, 'd', fontsize=12,fontweight='bold')
    plt.text(-0.75, 1.15, 'e', fontsize=12,fontweight='bold')
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    plt.savefig(\
    path_plot + 'Prediction_MMSBM_vs_Clustering.pdf',dpi=300)
    plt.show()
#=============================================================================================================

#4.1. Relative Distances
#===================================================================
if Measure=='Relative':
    
    #Means of means of pacients (relative accuracies, precisions, recalls)
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    Mean_Accuracies=[];Mean_Precisions=[];Mean_Recalls=[]
    SEM_Accuracies=[];SEM_Precisions=[];SEM_Recalls=[];

    Mean_Relative_Accuracies=[];Mean_Relative_Precisions=[];Mean_Relative_Recalls=[]
    SEM_Relative_Accuracies=[];SEM_Relative_Precisions=[];SEM_Relative_Recalls=[];

    for Dataset in Datasets:

        Relative_Accuracies=[];Relative_Precisions=[];Relative_Recalls=[]        
        for Patient in Differences[Dataset]:
            Relative_Accuracies.append(\
            Differences[Dataset][Patient]['Accuracy'])
            Relative_Precisions.append(\
            Differences[Dataset][Patient]['Precision'])
            Relative_Recalls.append(\
            Differences[Dataset][Patient]['Recall'])
        #------------------------------------------------------------
        Mean_Accuracies.append(np.mean(Relative_Accuracies))
        SEM_Accuracies.append(scipy.stats.sem(Relative_Accuracies,ddof=0))

        Mean_Precisions.append(np.mean(Relative_Precisions))
        SEM_Precisions.append(scipy.stats.sem(Relative_Precisions,ddof=0))

        Mean_Recalls.append(np.mean(Relative_Recalls))
        SEM_Recalls.append(scipy.stats.sem(Relative_Recalls,ddof=0))
        #-----------------------------------------------------------

    xticks_pos=[tick for tick in Datasets_label  ]
    plt.errorbar(xticks_pos, Mean_Accuracies, yerr=SEM_Accuracies, marker='o', label='Accuracy')
    plt.errorbar(xticks_pos, Mean_Precisions, yerr=SEM_Precisions, marker='o', label='Precision')
    plt.errorbar(xticks_pos, Mean_Recalls, yerr=SEM_Recalls, marker='o', label='Recall')
#    plt.title("Relative Difference of Predictions", fontsize=14)
    plt.xlabel('Dataset', fontsize= 14)
    plt.ylabel(r'$\langle$ Latent - Conventional $\rangle$',\
               fontsize=14 )
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend()

    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    plt.savefig(\
    path_plot + 'Prediction_Relative_Improvement.pdf',dpi=300)
    plt.show()



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
