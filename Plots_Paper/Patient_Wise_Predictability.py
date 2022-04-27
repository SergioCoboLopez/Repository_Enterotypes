#1/7/2019. Actualizacion del codigo para generar figuras para el paper. El codigo genera un panel de cada una de
#las metricas (Precision, Recall, Accuracy y Held-out-Loglikelihood)

import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
import seaborn as sns
import scipy.cluster
from scipy import stats
import math
import operator
import copy as cp
import time

NanSaver=1e-10


#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
class Patient:

    def __init__(self, name):
        self.name=name

    def membership(self, membership):
        self.membership=[float(element) for element in membership]

    def top_membership(self, membership):
        self.top_membership= self.membership.index( max(self.membership))

    def raw_conf_matrix(self, matrix):
        self.raw_conf_matrix=matrix


    def Precision(self):
        Predicted=self.raw_conf_matrix[1]
        self.Precision=self.raw_conf_matrix[1][1]/sum(Predicted)


    def Recall(self):
        Actual=np.transpose(self.raw_conf_matrix)[1]
        self.Recall=self.raw_conf_matrix[1][1]/sum(Actual)


    def Accuracy_Patients(self):
        Accuracy=0
        for Link in range(self.raw_conf_matrix.shape[1]):
            Accuracy+=self.raw_conf_matrix[Link][Link]/\
                            (sum(sum(self.raw_conf_matrix)))
        self.Accuracy_Patients=Accuracy
            
    def H_Patients(self,base):
        H=0
        for i in self.membership:
            if i==0:
                H+=0
            else:
                H+=-i*math.log(i,base)

        self.H_Patients=H
#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""


#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
class Microbe:

    def __init__(self, name):
        self.name=name
    
    def membership(self, membership):
        self.membership=[float(element) for element in membership]

    def top_membership(self, membership):
        self.top_membership= self.membership.index( max(self.membership))

    def raw_conf_matrix(self, matrix):
        self.raw_conf_matrix=matrix

    def Precision(self):
        Predicted=self.raw_conf_matrix[1]
        self.Precision=self.raw_conf_matrix[1][1]/sum(Predicted)


    def Recall(self):
        Actual=np.transpose(self.raw_conf_matrix)[1]
        self.Recall=self.raw_conf_matrix[1][1]/sum(Actual)


    def Accuracy_Microbes(self):
        Accuracy=0
        for Link in range(self.raw_conf_matrix.shape[1]):
            Accuracy+=self.raw_conf_matrix[Link][Link]/\
                            (sum(sum(self.raw_conf_matrix)))

        self.Accuracy_Microbes=Accuracy


    def H_Microbes(self,base):
        H=0
        for i in self.membership:
            if i==0:
                H+=0
            else:
                H+=-i*math.log(i,base)

        self.H_Microbes=H
#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

        
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

#0.3. Normalize Matrix
#======================================================
def Normalize_Matrix(Unnormalized_Matrix):
    nan_saver=1e-10
    Normalization=sum(Unnormalized_Matrix) + nan_saver
    Normalized_Matrix=np.zeros((rows,cols))
    for row in range(0,rows):
        for col in range(0,cols):
            Normalized_Matrix[row][col]=\
            Unnormalized_Matrix[row][col]/Normalization[col]

    return Normalized_Matrix
#======================================================

#0.4. Shannon Entropy
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

#0.5. Generate Entropy Vectors
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

#0.6. Total Accuracy, Precision, and recall
#=====================================================
def Total_Validations(Total_Confusion_Matrix):
    Accuracy_Tot=0;Precision_Tot=0;Recall_Tot=0

    #Accuracy
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::
    for Link in range(Total_Confusion_Matrix.shape[1]):
        Accuracy_Tot+=Total_Confusion_Matrix[Link][Link]/\
                   (sum(sum(Total_Confusion_Matrix)))
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::

    #Precision
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::
    Predicted=Total_Confusion_Matrix[1]
    Precision_Tot=Total_Confusion_Matrix[1][1]/sum(Predicted)
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::

    #Recall
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::
    Actual=np.transpose(Total_Confusion_Matrix)[1]
    Recall_Tot=Total_Confusion_Matrix[1][1]/sum(Actual)
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::

    return Accuracy_Tot, Precision_Tot, Recall_Tot
#=====================================================

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++


#1. PARAMETERS FOR INPUT/OUTPUT DATA
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#1.1. EXTERNAL PARAMETERS
#=====================
#Metric=sys.argv[1] #==
Number_of_Links=2  #==
#=====================

#1.2. INTERNAL PARAMETERS
#==================================================================
K=10;L=20

Datasets=['S-8_Tot', 'V-10_Tot', 'V-22_Tot', 'V-23_24_Tot', 'V-25_Tot'] 

#Datasets=['S-8_Tot', 'V-10_Tot']

Patients_per_Dataset=\
{'S-8_Tot':107, 'V-10_Tot':92, 'V-22_Tot':467, 'V-23_24_Tot':222,\
 'V-25_Tot':883}

Microbes_per_Dataset=\
{'S-8_Tot':128, 'V-10_Tot':137, 'V-22_Tot':134, 'V-23_24_Tot':118,\
 'V-25_Tot':144}
#==================================================================

#1.3. INITIALIZE ELEMENTS
#======================================
Patient_Precisions_Datasets={}
Patient_Recalls_Datasets={}
Patient_Accuracies_Datasets={}
Patient_H_O_LogLikelihood_Datasets={}
Patient_top_memberships_Datasets={}
Patient_Entropies_Datasets={}

Microbe_Precisions_Datasets={}
Microbe_Recalls_Datasets={}   
Microbe_Accuracies_Datasets={}
Microbe_top_memberships_Datasets={}
Microbe_Entropies_Datasets={}
Microbe_H_O_LogLikelihood_Datasets={}
#======================================

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#2. CALULATIONS FOR ALL DATASETS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
for Dataset in Datasets:

    print Dataset

    #2.1. PARAMETERS
    #==============================================================
    path_Out='../../Output_Data/'+ 'Leave_One_Out_' + Dataset + '/'
    path_In= '../../Input_Data/' + 'Leave_One_Out_' + Dataset + '/'


    Patients=Patients_per_Dataset[Dataset]
    Microbes=Microbes_per_Dataset[Dataset]
    #==============================================================

    #2.2. SHANNON ENTROPY
    #================================================================

    #2.2.1. Choose best likelihood
    #================================================================
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

    File='Dataset_' + Dataset + '_Seed_' + str(Seed) + '_Parameters.txt'
    Data=read_file(path_Out+File)
    #================================================================

    #2.2.2 Compute thetas, etas and aggregated
    #=================================================
    theta=[];contador_theta=0
    eta=[];contador_eta=0

    Patient_Classes=[]
    for line in Data[:Patients]:
        line.insert(0, contador_theta)

        #Clases
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~
        New_Patient=Patient(contador_theta)
        New_Patient.membership(line[1:])
        New_Patient.top_membership(line[1:])
        New_Patient.H_Patients(K)
        Patient_Classes.append(New_Patient)
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~

        contador_theta+=1

    theta=Data[:Patients]

    #Clases
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Patient_memberships=[pat.membership for pat in Patient_Classes]
    # Patient_top_memberships=[pat.top_membership for pat in Patient_Classes]
    # Patient_Entropies=[pat.H_Patients for pat in Patient_Classes]
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Microbe_Classes=[]
    for line in Data[Patients+1:Patients+1+Microbes]:
        line.insert(0, contador_eta )

        #Clases
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~
        New_Microbe=Microbe(contador_eta)
        New_Microbe.membership(line[1:])
        New_Microbe.top_membership(line[1:])
        New_Microbe.H_Microbes(L)
        Microbe_Classes.append(New_Microbe)
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~

        contador_eta+=1

    eta=Data[Patients+1:Patients+1+Microbes]

    #Clases
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Microbe_memberships=[micro.membership for micro in Microbe_Classes]
    # Microbe_top_memberships=[micro.top_membership for micro in Microbe_Classes]
    Microbe_Entropies=[micro.H_Microbes for micro in Microbe_Classes]
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    theta_aggregate=[0]*K
    for lineT in theta:
        for i in range(K):
            theta_aggregate[i]+=float(lineT[i+1])

    eta_aggregate=[0]*L
    for lineE in eta:
        for j in range(L):
            eta_aggregate[j]+=float(lineE[j+1])
    #=================================================

    #2.2.3. Get Entropy Vectors
    #======================================================
    EntropyPatients,Patients_Dict=List_Entropy( theta, K );
    EntropyMicrobes,Microbes_Dict=List_Entropy( eta  , L );
    #======================================================

    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    #2. GENERATE RESULTS
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    #Read best results from best likelihoods
    #--------------------------------------------
    Results=read_file(path_In + 'Best_Seeds.txt')
    #--------------------------------------------

    #2.1. Internal parameters and structures
    #==================================================================
    Confusion_Matrices_raw_dict={}
    Confusion_Matrices_dict={}

    rows=Number_of_Links;cols=Number_of_Links;

    Held_Out_LogLikelihood_Patients={}

    #2.1.1. MICROBES
    #-------------------------------------------------------------------------
    Held_Out_LogLikelihood_Microbes={}

    Confusion_Matrices_Microbes={}
    for Microbe_Conf in range(0,Microbes):
        Confusion_Matrices_Microbes[Microbe_Conf]=np.zeros((rows,cols))

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Confusion_Matrices_Microbes_3={}
    for Microbe_Conf in range(0,Microbes):
        Confusion_Matrices_Microbes_3[Microbe_Conf]=np.zeros((3,3))
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #-------------------------------------------------------------------------

    #==================================================================

    #2.2. Generate+write results
    #==============================================================================
    Confusion_Matrix_Total=np.zeros((rows,cols))
    for line_file in Results:

        #2.2.1. Read Input files
        #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        Patient_R=line_file[0];Fold=line_file[1];Seed=line_file[2]
        Parameters = str(Patient_R) + "_Fold_" + str(Fold) +  "_Seed_" + str(Seed)

        Scores = path_Out + Parameters +  '_scores.txt'
        raw_scores=read_file(Scores)
        Check=read_file(path_In + 'P_' + Patient_R + '_F_'+Fold+ '_TestCheck.txt')

        #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        #2.2.2. Read and compute MMSBM scores
        #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        #2.2.2.1. Clean scores file
        #----------------------------
        scores=clean_data(raw_scores)
        #----------------------------

        #2.2.2.2. Collect scores in dictionary
        #---------------------------------------------------------------------------------    
        Held_Out_LogLikelihood_per_Patient=0
        scores_dict={};scores_dict_Accuracy={}

        for line in scores:

            Microbe_R=line[0];Person=line[1];Bet=line[2:]
            Max_score=max(Bet);Pareja_Tupla=(Microbe_R, Person)

            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            #TRES LINKS: Accuracy
            #.......................................................
            Winner_Accuracy=Bet.index(max(Bet))
            scores_dict_Accuracy[ Pareja_Tupla ]=Winner_Accuracy#Assign prediction
            #.......................................................

            #DOS LINKS: Precision y Recall
            #.......................................................
            New_line=[ Microbe_R, Person, Bet[0], Bet[1]+Bet[2] ]
            Winner=New_line[2:].index(max(New_line[2:]))
            scores_dict[ Pareja_Tupla ]=Winner
            #.......................................................
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


            #Compute individual Held_Out_LogLikelihood for patient (maximum score)
            #.....................................................................
            Held_Out_LogLikelihood_Pair=math.log(Max_score + NanSaver )
            Held_Out_LogLikelihood_per_Patient+=Held_Out_LogLikelihood_Pair
            #.....................................................................


            #Held_Out_LogLikelihood MICROBES
            #.........................................................................
            try:
                Held_Out_LogLikelihood_Microbes[Microbe_R][Fold]+=Held_Out_LogLikelihood_Pair
            except KeyError:
                try:
                    Held_Out_LogLikelihood_Microbes[Microbe_R][Fold]=Held_Out_LogLikelihood_Pair
                except KeyError:
                    Held_Out_LogLikelihood_Microbes[Microbe_R]={Fold:Held_Out_LogLikelihood_Pair}
            #.........................................................................
            
        #Save Held_Out_LogLikelihood in dictionary with patients and folds
        #.......................................................................................
        try:
            Held_Out_LogLikelihood_Patients[Patient_R][Fold]=Held_Out_LogLikelihood_per_Patient
        except KeyError:
            Held_Out_LogLikelihood_Patients[Patient_R]={Fold:Held_Out_LogLikelihood_per_Patient}
        #........................................................................................

        #---------------------------------------------------------------------------------

        #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


        #2.2.3. Build individual (per Fold) Confusion Matrices and Dictionaries
        #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        #2.2.3.1. Confusion Matrices
        #----------------------------------------------------------------------
        Matrix_Counter_perFold=np.zeros((rows,cols))

        #~~~~~~~~~~~~~~~~~~~~~
        Counter_Anomalies=0
        #~~~~~~~~~~~~~~~~~~~~~

        for line in Check:
            MicrobeCheck=float(line[0]); PersonCheck=float(line[1])
            Abundance=float(line[2])        
            TuplaCheck=(MicrobeCheck,PersonCheck)


            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if scores_dict_Accuracy[TuplaCheck] != scores_dict[TuplaCheck]:
                Counter_Anomalies+=1

            Confusion_Matrices_Microbes_3[MicrobeCheck]\
                [scores_dict_Accuracy[TuplaCheck],int(Abundance)]+=1
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            #DOS LINKS
            #.......................................................
            if Number_of_Links==2:
                if Abundance==2:
                    Abundance=1.0
            #.......................................................

            #Individual Confusion Matrix
            #-------------------------------------------------------------------
            Matrix_Counter_perFold[scores_dict[TuplaCheck],int(Abundance)]+=1
            #-------------------------------------------------------------------

            Confusion_Matrices_Microbes[MicrobeCheck]\
                [scores_dict[TuplaCheck],int(Abundance)]+=1
        #----------------------------------------------------------------------

        #2.2.3.2. Normalize Confusion Matrix
        #----------------------------------------------------------------
        Conf_Mat_Patient_Fold=Normalize_Matrix(Matrix_Counter_perFold)
        #----------------------------------------------------------------

        #2.2.3.3. Fill dictionaries
        #-----------------------------------------------------------------------  
        try:
            Confusion_Matrices_raw_dict[Patient_R][Fold]=Matrix_Counter_perFold
            Confusion_Matrices_dict[Patient_R][Fold]=Conf_Mat_Patient_Fold
        except KeyError:

            try:
                Confusion_Matrices_raw_dict[Patient_R].append(Fold)
                Confusion_Matrices_dict[Patient_R].append(Fold)
            except KeyError:
                Confusion_Matrices_raw_dict[Patient_R]=\
                {Fold:Matrix_Counter_perFold}
                Confusion_Matrices_dict[Patient_R]=\
                {Fold:Conf_Mat_Patient_Fold}
        #-----------------------------------------------------------------------

        #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #==============================================================================

    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++

    #3. CREATE CONFUSION MATRICES
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    #3.1. Build confusion matrices per patient and compute precision,recall,accuracy
    #===============================================================================
    Confusion_Matrices_per_Patient={}
    for Patient_Matrix in Confusion_Matrices_dict:

        #Hacer la media de las matrices (medias) de cada fold
        #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        Confusion_Matrix_Patient_Avg=\
            sum(Confusion_Matrices_dict[Patient_Matrix].values())

        Confusion_Matrix_Patient_Avg=Normalize_Matrix(Confusion_Matrix_Patient_Avg)
        #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        #Hacer la suma de cada fold
        #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        Confusion_Matrix_Patient_Raw=\
            sum(Confusion_Matrices_raw_dict[Patient_Matrix].values())
        #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        #Normalizar el anterior
        #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        Confusion_Matrix_Patient_Agg=Normalize_Matrix(Confusion_Matrix_Patient_Raw)
        #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        try:        
            Confusion_Matrices_per_Patient[Patient_Matrix]['Raw']=\
                                        Confusion_Matrix_Patient_Raw

            Confusion_Matrices_per_Patient[Patient_Matrix]['Average']=\
                                      Confusion_Matrix_Patient_Avg

            Confusion_Matrices_per_Patient[Patient_Matrix]['Aggregated']=\
                                      Confusion_Matrix_Patient_Agg
        except KeyError:
                Confusion_Matrices_per_Patient[Patient_Matrix]=\
                {'Raw':Confusion_Matrix_Patient_Raw, \
                 'Average':Confusion_Matrix_Patient_Avg,\
                 'Aggregated':Confusion_Matrix_Patient_Agg}



        Patient_Classes[int(Patient_Matrix)].\
            raw_conf_matrix( Confusion_Matrices_per_Patient[Patient_Matrix]['Raw'] )
        Patient_Classes[int(Patient_Matrix)].Precision()
        Patient_Classes[int(Patient_Matrix)].Recall()
        Patient_Classes[int(Patient_Matrix)].Accuracy_Patients()
    #===============================================================================
    
    #Clases
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Patient_Confusions=[pat.raw_conf_matrix for pat in Patient_Classes]
    Patient_Precisions=[pat.Precision for pat in Patient_Classes]
    Patient_Recalls=   [pat.Recall    for pat in Patient_Classes]
    Patient_Accuracies=[pat.Accuracy_Patients  for pat in Patient_Classes]
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    #3.2. Build Precision, Recall, Accuracy for each microbe
    #===================================================================

    #3.2.2. Microbes
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    for Microbe_Tot in Confusion_Matrices_Microbes:

        Microbe_Classes[int(Microbe_Tot)].\
        raw_conf_matrix( Confusion_Matrices_Microbes[Microbe_Tot] )    
        Microbe_Classes[int(Microbe_Tot)].Precision()
        Microbe_Classes[int(Microbe_Tot)].Recall()
        Microbe_Classes[int(Microbe_Tot)].Accuracy_Microbes()
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #Clases
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Microbe_Confusions=[mic.raw_conf_matrix for mic in Microbe_Classes]
    Microbe_Precisions=[mic.Precision for mic in Microbe_Classes]
    Microbe_Recalls=   [mic.Recall    for mic in Microbe_Classes]
    Microbe_Accuracies=[mic.Accuracy_Microbes  for mic in Microbe_Classes]
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    #===================================================================


    #3.3. Global matrices over the whole sample
    #=================================================================
    Confusion_Matrix_Avg_Total=np.zeros((cols,rows))
    Confusion_Matrix_Agg_Total=np.zeros((cols,rows))
    Confusion_Matrix_Raw_Total=np.zeros((cols,rows))

    for p in Confusion_Matrices_per_Patient:
        Confusion_Matrix_Raw_Total+=Confusion_Matrices_per_Patient[p]['Raw']
        Confusion_Matrix_Agg_Total+=Confusion_Matrices_per_Patient[p]['Aggregated']
        Confusion_Matrix_Avg_Total+=Confusion_Matrices_per_Patient[p]['Average']

    Confusion_Matrix_Avg_Total=Normalize_Matrix( Confusion_Matrix_Avg_Total )
    Confusion_Matrix_Agg_Total=Normalize_Matrix( Confusion_Matrix_Agg_Total )

    print Confusion_Matrix_Raw_Total
    Raw_Accuracy_Total, Raw_Precision_Total, Raw_Recall_Total=Total_Validations(Confusion_Matrix_Raw_Total)
    print "Total Accuracy: "  ,Raw_Accuracy_Total
    print "Total Precision: " ,Raw_Precision_Total
    print "Total Recall: "    ,Raw_Recall_Total


    Confusion_Matrix_Raw_Total=Normalize_Matrix( Confusion_Matrix_Raw_Total )
    print Confusion_Matrix_Raw_Total
    Raw_Accuracy_Total, Raw_Precision_Total, Raw_Recall_Total=Total_Validations(Confusion_Matrix_Raw_Total)
    print Raw_Accuracy_Total
    print Raw_Precision_Total
    print Raw_Recall_Total
    #=================================================================

    #3.4. Held-Out-LogLikelihood
    #================================================================
    Held_Out_LogLikelihood_Patients_Average=[]
    for Patient_L in range(Patients):
         Held_Out_LogLikelihood_Patients_Average.append(\
         sum(Held_Out_LogLikelihood_Patients[str(Patient_L)].values())/Patients ) 

    Held_Out_LogLikelihood_Microbes_Average=[]
    for Microbe_L in range(Microbes):
        Held_Out_LogLikelihood_Microbes_Average.append(\
        sum(Held_Out_LogLikelihood_Microbes[Microbe_L].values())/Microbes ) 

    #================================================================


    #3.5. Correlations coefficients
    #======================================================================================================

    #3.5.1. Remove nans from microbes precision
    #------------------------------------------------------------------
    No_nan_Microbe_Entropies=[];No_nan_Microbe_Precisions=[]
    for element in zip(Microbe_Entropies,Microbe_Precisions):

        if math.isnan(element[1]) == True:
            continue
        else:
            No_nan_Microbe_Entropies.append( element[0] )
            No_nan_Microbe_Precisions.append( element[1] )
    #------------------------------------------------------------------

    #======================================================================================================


    #3.6. Save dataset-vectors to vectors
    #======================================================================================================

    #3.6.1. Patient Performance
    #--------------------------------------------------------------------
    Patient_Precisions_Datasets[Dataset]=Patient_Precisions
    Patient_Recalls_Datasets[Dataset]=   Patient_Recalls
    Patient_Accuracies_Datasets[Dataset]=Patient_Accuracies
    Patient_H_O_LogLikelihood_Datasets[Dataset]=Held_Out_LogLikelihood_Patients_Average
    # Patient_top_memberships_Datasets[Dataset]=Patient_top_memberships
    # Patient_Entropies_Datasets[Dataset]=Patient_Entropies
    #--------------------------------------------------------------------
    
    #3.6.3. Microbe Performance
    #--------------------------------------------------------------------
    Microbe_Precisions_Datasets[Dataset]=Microbe_Precisions
    Microbe_Recalls_Datasets[Dataset]=   Microbe_Recalls
    Microbe_Accuracies_Datasets[Dataset]=Microbe_Accuracies
    Microbe_H_O_LogLikelihood_Datasets[Dataset]=Held_Out_LogLikelihood_Microbes_Average
    # Microbe_top_memberships_Datasets[Dataset]=Microbe_top_memberships
    # Microbe_Entropies_Datasets[Dataset]=Microbe_Entropies
    #--------------------------------------------------------------------
    
    #======================================================================================================

    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

print Patient_Confusions[57]
print Patient_Precisions[57]
print Patient_Recalls[57]
