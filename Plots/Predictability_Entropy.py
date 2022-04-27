#20/11/2018. Este codigo calcula las correlaciones entre medidas de predictabilidad (precision, recall y accuracy)
#y la entropia de shannon de los vectores de membership de microbios y pacientes

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

NanSaver=1e-10


#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
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
#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""


#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
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
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#1.1. EXTERNAL PARAMETERS
#=====================
Dataset=sys.argv[1]#==
Number_of_Links=2  #==
#=====================

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
path_In= '../../Input_Data/' + 'Leave_One_Out_' + Dataset + '/'
#==============================================================

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#2. SHANNON ENTROPY
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#2.1. Choose best likelihood
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

File='Dataset_' + Dataset + '_Seed_' + str(Seed) + '_Parameters.txt'
Data=read_file(path_Out+File)
#======================================================================

#2.2. Compute thetas, etas and aggregated
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
Patient_memberships=[pat.membership for pat in Patient_Classes]
Patient_top_memberships=[pat.top_membership for pat in Patient_Classes]
Patient_Entropies=[pat.H_Patients for pat in Patient_Classes]
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
Microbe_memberships=[micro.membership for micro in Microbe_Classes]
Microbe_top_memberships=[micro.top_membership for micro in Microbe_Classes]
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

#1.3. Get Entropy Vectors
#======================================================
EntropyPatients,Patients_Dict=List_Entropy( theta, K );
EntropyMicrobes,Microbes_Dict=List_Entropy( eta  , L );
#======================================================

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#2. GENERATE RESULTS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#Read best results from best likelihoods
#--------------------------------------------
Results=read_file(path_In + 'Best_Seeds.txt')
#--------------------------------------------

#2.1. Internal parameters and structures
#==================================================================
Confusion_Matrices_raw_dict={}
Confusion_Matrices_dict={}

rows=Number_of_Links;cols=Number_of_Links;

Held_Out_LogLikelihood={}

#2.1.1. MICROBES
#--------------------------------------------------------------
Confusion_Matrices_Microbes={}
for Microbe in range(0,Microbes):
    Confusion_Matrices_Microbes[Microbe]=np.zeros((rows,cols))
#--------------------------------------------------------------

#==================================================================

#2.2. Generate+write results
#==============================================================================
Confusion_Matrix_Total=np.zeros((rows,cols))
for line_file in Results:
    
    #2.2.1. Read Input files
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    Patient=line_file[0];Fold=line_file[1];Seed=line_file[2]
    Parameters = str(Patient) + "_Fold_" + str(Fold) +  "_Seed_" + str(Seed)
    
    Scores = path_Out + Parameters +  '_scores.txt'
    raw_scores=read_file(Scores)
    Check=read_file(path_In + 'P_' + Patient + '_F_'+Fold+ '_TestCheck.txt')

    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #2.2.2. Read and compute MMSBM scores
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #2.2.2.1. Clean scores file
    #----------------------------
    scores=clean_data(raw_scores)
    #----------------------------

    #2.2.2.2. Collect scores in dictionary
    #---------------------------------------------------------------------------------    
    Held_Out_LogLikelihood_Patient=0
    scores_dict={}
    for line in scores:
        
        Microbe=line[0];Person=line[1];Bet=line[2:]
        Max_score=max(Bet);Pareja_Tupla=(Microbe, Person)

        if Number_of_Links==3:

            #.......................................................
	    Winner=Bet.index(max(Bet))
	    scores_dict[ Pareja_Tupla ]=Winner#Assign prediction
            #.......................................................
                
        elif Number_of_Links==2:
            #DOS LINKS
            #.......................................................
            New_line=[ Microbe, Person, Bet[0], Bet[1]+Bet[2] ]
            Winner=New_line[2:].index(max(New_line[2:]))
            scores_dict[ Pareja_Tupla ]=Winner
            #.......................................................

        #Compute Held_Out_LogLikelihood for individual patient (maximum score)
        #.....................................................................
        Held_Out_LogLikelihood_Patient+=math.log(Max_score + NanSaver )
        #.....................................................................

    #Save Held_Out_LogLikelihood in dictionary with patients and folds
    #.........................................................................
    try:
        Held_Out_LogLikelihood[Patient][Fold]=Held_Out_LogLikelihood_Patient
    except KeyError:
        Held_Out_LogLikelihood[Patient]={Fold:Held_Out_LogLikelihood_Patient}
    #.........................................................................

    #---------------------------------------------------------------------------------
    
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    
    #2.2.3. Build individual (per Fold) Confusion Matrices and Dictionaries
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #2.2.3.1. Confusion Matrices
    #----------------------------------------------------------------------
    Matrix_Counter_perFold=np.zeros((rows,cols))
    
    for line in Check:
        MicrobeCheck=float(line[0]); PersonCheck=float(line[1])
        Abundance=float(line[2])        
        TuplaCheck=(MicrobeCheck,PersonCheck)

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
        Confusion_Matrices_raw_dict[Patient][Fold]=Matrix_Counter_perFold
        Confusion_Matrices_dict[Patient][Fold]=Conf_Mat_Patient_Fold
    except KeyError:
        
        try:
            Confusion_Matrices_raw_dict[Patient].append(Fold)
            Confusion_Matrices_dict[Patient].append(Fold)
        except KeyError:
            Confusion_Matrices_raw_dict[Patient]=\
            {Fold:Matrix_Counter_perFold}
            Confusion_Matrices_dict[Patient]=\
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

            
    Patient_Classes[int(Patient_Matrix)].raw_conf_matrix( Confusion_Matrices_per_Patient[Patient_Matrix]['Raw'] )
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
#===============================================================================================
Held_Out_LogLikelihood_Average=[]
for Patient_L in range(Patients):
     Held_Out_LogLikelihood_Average.append( np.mean(Held_Out_LogLikelihood[str(Patient_L)].values()) )     
#===============================================================================================


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

#3.5.2. Compute correlation parameters
#-----------------------------------------------------------------------------------------------------
Spearman_Patients_Recall=stats.spearmanr(Patient_Entropies,Patient_Recalls)
Spearman_Patients_Precision=stats.spearmanr(Patient_Entropies,Patient_Precisions)
Spearman_Patients_Accuracy=stats.spearmanr(Patient_Entropies,Patient_Accuracies)
Spearman_Patients_H_O_LogLikelihood=stats.spearmanr(Patient_Entropies,Held_Out_LogLikelihood_Average)


Spearman_Microbes_Recall=stats.spearmanr(Microbe_Entropies,Microbe_Recalls)
Spearman_Microbes_Precision=stats.spearmanr(No_nan_Microbe_Entropies,No_nan_Microbe_Precisions)
Spearman_Microbes_Accuracy=stats.spearmanr(Microbe_Entropies,Microbe_Accuracies)
#-----------------------------------------------------------------------------------------------------
#======================================================================================================

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#4. PLOTS
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#4.0. Held Out LogLikelihood
#================================================================================================
plt.scatter(Held_Out_LogLikelihood_Average,Patient_Entropies)
plt.text(min(Held_Out_LogLikelihood_Average) , 0.9*max(Patient_Entropies), \
         'Spearman: %f' %Spearman_Patients_H_O_LogLikelihood[0] + '\n'+ 'p-value: %f' \
         %Spearman_Patients_H_O_LogLikelihood[1] , style='normal', bbox={'facecolor':'gray',\
                                                                         'alpha':0.5, 'pad':5})
plt.xlabel("Held-Out-LogLikelihood")
plt.ylabel("Shannon Entropy")
plt.title("Held_Out_LogLikelihood for Dataset %s" %Dataset)
plt.savefig(\
'../../Plots/First_Path/Held_Out_LogLikelihood_' + Dataset + '.png',dpi=300)
plt.show()
#================================================================================================

#4.1. 
#============================================================================

#4.1.1. Fontsizes and Parameters
#=============================
size_eje=12
size_ticks=13
size_title=15
size_letter=17
#=============================

fig=plt.figure(figsize=(15,10))
gs = gridspec.GridSpec(2, 3)
gs.update(left=0.1,right=0.95,bottom=0.1,top=0.94,wspace=0.25,hspace=0.3)

colors_patients = cm.Set3(np.linspace(0, 1, K))
colors_microbes = cm.Paired(np.linspace(0, 1, L/2))
                    
colors_shapes_microbes={}
shapes=["v", "^"]
counter=0
for group in range(L/2):
    for shape in range(len(shapes)):
        colors_shapes_microbes[counter]=( colors_microbes[group], shapes[shape] )
        counter+=1

for element in colors_shapes_microbes:
    print element, colors_shapes_microbes[element]
        
#(0,0) Plot
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Recall_Patients= plt.subplot(gs[0, 0])
for instance in Patient_Classes:
    plt.scatter(instance.H_Patients, instance.Recall, color=colors_patients[ instance.top_membership ], \
                alpha=1 - 0*instance.H_Patients,marker="s")

plt.title(r'$H_{patients}$ vs Recall in Dataset %s' %Dataset)
plt.xlabel(r'$H_{patients}$')
plt.ylabel('Recall')

#Legend
#------------------------------------------------
legend_x_1=min(Patient_Entropies)
legend_y_1=min(Patient_Recalls)

Recall_Patients.text(legend_x_1, legend_y_1, 'Spearman: %f' %Spearman_Patients_Recall[0] + '\n'+ 'p-value: %f' %Spearman_Patients_Recall[1] , style='normal', bbox={'facecolor':'gray', 'alpha':0.5, 'pad':5})
#------------------------------------------------

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#(0,1) Plot
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Precision_Patients= plt.subplot(gs[0, 1])
for instance in Patient_Classes:
    plt.scatter(instance.H_Patients, instance.Precision, color=colors_patients[ instance.top_membership ], \
                alpha=1 - 0*instance.H_Patients,marker="s")
plt.title(r'$H_{patients}$ vs Precision in Dataset %s' %Dataset)
plt.xlabel(r'$H_{patients}$')
plt.ylabel('Precision')
#Legend
#------------------------------------------------
legend_x_2=min(Patient_Entropies)
legend_y_2=min(Patient_Precisions)
Precision_Patients.text(legend_x_2, legend_y_2, 'Spearman: %f' %Spearman_Patients_Precision[0] + '\n'+ 'p-value: %f' %Spearman_Patients_Precision[1] , style='normal', bbox={'facecolor':'gray', 'alpha':0.5, 'pad':5})
#------------------------------------------------
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#(0,2) Plot
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Accuracy_Patients= plt.subplot(gs[0, 2])
for instance in Patient_Classes:
    plt.scatter(instance.H_Patients, instance.Accuracy_Patients, color=colors_patients[ instance.top_membership ], \
                alpha=1 - 0*instance.H_Patients, marker="s")
plt.title(r'$H_{patients}$ vs Accuracy in Dataset %s' %Dataset)
plt.xlabel(r'$H_{patients}$')
plt.ylabel('Accuracy')
#Legend
#------------------------------------------------
legend_x_3=min(Patient_Entropies)
legend_y_3=min(Patient_Accuracies)
Accuracy_Patients.text(legend_x_3, legend_y_3, 'Spearman: %f' %Spearman_Patients_Accuracy[0] + '\n'+ 'p-value: %f' %Spearman_Patients_Accuracy[1] , style='normal', bbox={'facecolor':'gray', 'alpha':0.5, 'pad':5})
#------------------------------------------------
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#(1,0) Plot
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Recall_Microbes= plt.subplot(gs[1, 0])
for instance in Microbe_Classes:
    plt.scatter(instance.H_Microbes, instance.Recall, color=colors_shapes_microbes[ instance.top_membership ][0],
                alpha=1 - 0*instance.H_Microbes,marker=colors_shapes_microbes[ instance.top_membership ][1] )
    
plt.title(r'$H_{microbes}$ vs Recall in Dataset %s' %Dataset)
plt.xlabel(r'$H_{microbes}$')
plt.ylabel('Recall')
#Legend
#------------------------------------------------
legend_x_3=min(Microbe_Entropies)
legend_y_3=min(Microbe_Recalls)
Recall_Microbes.text(legend_x_3, legend_y_3, \
'Spearman: %f' %Spearman_Microbes_Recall[0] + '\n'+ 'p-value: %f' %Spearman_Microbes_Recall[1] \
, style='normal', bbox={'facecolor':'gray', 'alpha':0.5, 'pad':5})
#------------------------------------------------
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#(1,1) Plot
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Precision_Microbes= plt.subplot(gs[1, 1])
for instance in Microbe_Classes:
    plt.scatter(instance.H_Microbes, instance.Precision, color=colors_shapes_microbes[instance.top_membership][0],
                alpha=1 - 0*instance.H_Microbes,marker=colors_shapes_microbes[ instance.top_membership ][1] )

plt.title(r'$H_{microbes}$ vs Precision in Dataset %s' %Dataset)
plt.xlabel(r'$H_{microbes}$')
plt.ylabel('Precision')
#Legend
#------------------------------------------------
legend_x_4=min(Microbe_Entropies)
legend_y_4=min(No_nan_Microbe_Precisions)

Precision_Microbes.text(legend_x_4, legend_y_4, 'Spearman: %f' %Spearman_Microbes_Precision[0] + '\n'+ 'p-value: %f' %Spearman_Microbes_Precision[1] , style='normal', bbox={'facecolor':'gray', 'alpha':0.5, 'pad':5})
#------------------------------------------------

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#(1,2) Plot
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Accuracy_Microbes= plt.subplot(gs[1, 2])
for instance in Microbe_Classes:
    plt.scatter(instance.H_Microbes, instance.Accuracy_Microbes, \
    color=colors_shapes_microbes[instance.top_membership][0],alpha=1 - 0*instance.H_Microbes,\
                marker=colors_shapes_microbes[ instance.top_membership ][1] )
    
plt.title(r'$H_{microbes}$ vs Accuracy in Dataset %s' %Dataset)
plt.xlabel(r'$H_{microbes}$')
plt.ylabel('Accuracy')
#Legend
#------------------------------------------------
legend_x_5=min(Microbe_Entropies)
legend_y_5=min(Microbe_Accuracies)
Accuracy_Microbes.text(legend_x_5, legend_y_5, 'Spearman: %f' %Spearman_Microbes_Accuracy[0] + '\n'+ 'p-value: %f' %Spearman_Microbes_Accuracy[1] , style='normal', bbox={'facecolor':'gray', 'alpha':0.5, 'pad':5})
#------------------------------------------------

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


plt.savefig(\
'../../Plots/First_Path/Entropy_vs_Precision_Recall_'+Dataset+'_Colors.pdf',dpi=300)
plt.show()
#============================================================================

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
