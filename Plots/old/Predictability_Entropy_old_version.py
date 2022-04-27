#20/11/2018. Con este script pretendemos ver si existe una correlacion entre la
#entropia de Shannon de los pacientes y su predictabilidad.

import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import scipy.cluster
import math

NanSaver=1e-10

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
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++


#1.GENERATE CONFUSION MATRICES
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#1.1. Parameters for reading data
#========================================================

#Seed=np.random.seed(1111)#--->Para escoger la muestra de 9 pacientes

#1.1.1. EXTERNAL PARAMETERS
#---------------------
Dataset=sys.argv[1]#--
Number_of_Links=3  #--
#---------------------

Patients_per_Dataset=\
{'S-8_Tot':107, 'V-10_Tot':92, 'V-22_Tot':467, 'V-23_24_Tot':222,\
 'V-25_Tot':883}

Microbes_per_Dataset=\
{'S-8_Tot':128, 'V-10_Tot':137, 'V-22_Tot':134, 'V-23_24_Tot':118,\
 'V-25_Tot':144}

Patients=Patients_per_Dataset[Dataset]
Microbes=Microbes_per_Dataset[Dataset]


#Paths for files
#------------------------------------------------------------
path_Out='../../Output_Data/'+ 'Leave_One_Out_' + Dataset + '/'
path_In= '../../Input_Data/' + 'Leave_One_Out_' + Dataset + '/'
#------------------------------------------------------------

#Read best results from best likelihoods
#--------------------------------------------
Results=read_file(path_In + 'Best_Seeds.txt')
#--------------------------------------------
#========================================================

#1.2. Internal parameters and structures
#========================================================
Confusion_Matrices_raw_dict={}
Confusion_Matrices_dict={}
rows=Number_of_Links;cols=Number_of_Links;
Link_Counter=[0]*Number_of_Links#->We count links present in the 30 pat.sample

#1.2.1. MICROBES
#--------------------------------------------
Confusion_Matrices_Microbes={}
for Microbe in range(0,Microbes):
    Confusion_Matrices_Microbes[Microbe]=np.zeros((rows,cols))
#--------------------------------------------

#========================================================

#1.3. Generate+write results
#==============================================================================
Confusion_Matrix_Total=np.zeros((rows,cols))
for line in Results:
    
    #1.3.1. Read Input files
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    Patient=line[0];Fold=line[1];Seed=line[2]
    Parameters = str(Patient) + "_Fold_" + str(Fold) +  "_Seed_" + str(Seed)
    
    Scores = path_Out + Parameters +  '_scores.txt'
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
    
    for line in scores:
        Microbe=line[0];Person=line[1];Bet=line[2:];Winner=Bet.index(max(Bet))

        Pareja_Tupla=(Microbe, Person)

        scores_dict[ Pareja_Tupla ]=Winner#Assign prediction
    #------------------------------------------
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
    #1.3.3. Build individual (per Fold) Confusion Matrices and Dictionaries
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #1.3.3.1. Confusion Matrices
    #----------------------------------------------------------------------
    Conf_Mat_Patient_Fold_raw=np.zeros((rows,cols))
    
    for line in Check:
        MicrobeCheck=float(line[0]); PersonCheck=float(line[1])
        Abundance=float(line[2])
        TuplaCheck=(MicrobeCheck,PersonCheck)
        
        Link_Counter[int(Abundance)]+=1#Abundance Counter

        #Individual Confusion Matrix
        #---------------------------------------------------------------
        Conf_Mat_Patient_Fold_raw[scores_dict[TuplaCheck],int(Abundance)]+=1
        #---------------------------------------------------------------

        Confusion_Matrices_Microbes[MicrobeCheck]\
            [scores_dict[TuplaCheck],int(Abundance)]+=1
    #----------------------------------------------------------------------

    #1.3.3.2. Normalize Confusion Matrix
    #----------------------------------------------------------------
    Conf_Mat_Patient_Fold=Normalize_Matrix(Conf_Mat_Patient_Fold_raw)
    #----------------------------------------------------------------


    
    #1.3.3.3. Fill dictionaries
    #-----------------------------------------------------------------------  
    try:
        Confusion_Matrices_raw_dict[Patient][Fold]=Conf_Mat_Patient_Fold_raw
        Confusion_Matrices_dict[Patient][Fold]=Conf_Mat_Patient_Fold
    except KeyError:
        try:
            Confusion_Matrices_raw_dict[Patient].append(Fold)
            Confusion_Matrices_dict[Patient].append(Fold)
        except KeyError:
            Confusion_Matrices_raw_dict[Patient]=\
            {Fold:Conf_Mat_Patient_Fold_raw}
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

#Raw, average, and aggregated confusion matrices
#===================================================================
Confusion_Matrices_Averaged={}
Confusion_Matrices_Aggregated={}
Confusion_Matrices_Raw={}

for Patient_Matrix in Confusion_Matrices_dict:
    Confusion_Matrix_Patient_Avg=np.zeros((cols,rows))
    Confusion_Matrix_Patient_Agg=np.zeros((cols,rows))
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    for Fold in Confusion_Matrices_dict[Patient_Matrix]:
        Confusion_Matrix_Patient_Avg+=\
            Confusion_Matrices_dict[Patient_Matrix][Fold]
        
        Confusion_Matrix_Patient_Agg+=\
            Confusion_Matrices_raw_dict[Patient_Matrix][Fold]
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    Confusion_Matrices_Raw[Patient_Matrix]=Confusion_Matrix_Patient_Agg
    
    Confusion_Matrix_Patient_Avg=Normalize_Matrix(Confusion_Matrix_Patient_Avg)
    Confusion_Matrix_Patient_Agg=Normalize_Matrix(Confusion_Matrix_Patient_Agg)
    
    Confusion_Matrices_Averaged[Patient_Matrix]=Confusion_Matrix_Patient_Avg
    Confusion_Matrices_Aggregated[Patient_Matrix]=Confusion_Matrix_Patient_Agg
    

Confusion_Matrix_Avg_Total=\
Normalize_Matrix( sum( Confusion_Matrices_Averaged.values() ) )

Confusion_Matrix_Agg_Total=\
Normalize_Matrix( sum( Confusion_Matrices_Aggregated.values() ) )    
#===================================================================

#Precision and recall
#===================================================================
Precision_Recall={}
for Patient_Prec_Rcll in Confusion_Matrices_Raw:
    Precision_Recall[Patient_Prec_Rcll]={'Precision':{}, 'Recall':{}}
    for Link in range(Number_of_Links):
        Actual=Confusion_Matrices_Raw[Patient_Prec_Rcll][Link]
        Predicted=np.transpose(Confusion_Matrices_Raw[Patient_Prec_Rcll])[Link]
        
        Precision_Recall[Patient_Prec_Rcll]['Precision'][Link]=\
        Confusion_Matrices_Raw[Patient_Prec_Rcll][Link][Link]/sum(Predicted)

        Precision_Recall[Patient_Prec_Rcll]['Recall'][Link]=\
        Confusion_Matrices_Raw[Patient_Prec_Rcll][Link][Link]/sum(Actual)

        
Precision_Recall_Microbes={}
for Microbe_Prec_Rcll in Confusion_Matrices_Microbes:
    Precision_Recall_Microbes[Microbe_Prec_Rcll]={'Precision':{}, 'Recall':{}}
    for Link in range(Number_of_Links):
        Actual=Confusion_Matrices_Microbes[Microbe_Prec_Rcll][Link]
        Predicted=\
        np.transpose(Confusion_Matrices_Microbes[Microbe_Prec_Rcll])[Link]
        
        Precision_Recall_Microbes[Microbe_Prec_Rcll]['Precision'][Link]=\
        Confusion_Matrices_Microbes[Microbe_Prec_Rcll][Link][Link]/\
        (sum(Predicted)+NanSaver)

        Precision_Recall_Microbes[Microbe_Prec_Rcll]['Recall'][Link]=\
        Confusion_Matrices_Microbes[Microbe_Prec_Rcll][Link][Link]/\
        (sum(Actual)+NanSaver)
#===================================================================

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#4. SHANNON ENTROPY
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
K=10;L=20

Likelihood_Seed={}
for Seed in range(1,10,2):
    Likelihood_File=\
    'Dataset_' + Dataset + '_Seed_' + str(Seed) + '_LogLikelihood.txt'
    LikelihoodFinal=read_file(path_Out + Likelihood_File)[-1]
    print LikelihoodFinal
    Likelihood_Seed[float(LikelihoodFinal[0])]=Seed

print Likelihood_Seed
print max(Likelihood_Seed)
print Likelihood_Seed[max(Likelihood_Seed)]
Seed=Likelihood_Seed[max(Likelihood_Seed)] 


File='Dataset_' + Dataset + '_Seed_' + str(Seed) + '_Parameters.txt'
Data=read_file(path_Out+File)

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

print contador_theta
print contador_eta
theta_aggregate=[0]*K
for lineT in theta:
#    print lineT
    for i in range(K):
#        print i
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

#1.4. Build final dictionaries containing all the information
#===========================================================
Final_Dictionary={}
for Patient_Final_Dict in Precision_Recall:
   
    Final_Dictionary[Patient_Final_Dict]=\
    {'H_pat':{},'Precision':{},'Recall':{}}
    
    Final_Dictionary[Patient_Final_Dict]['H_pat']=\
            Patients_Dict[int(Patient_Final_Dict)][-1]
    
    Final_Dictionary[Patient_Final_Dict]['Precision']=\
            Precision_Recall[Patient_Final_Dict]['Precision']
    
    Final_Dictionary[Patient_Final_Dict]['Recall']=\
            Precision_Recall[Patient_Final_Dict]['Recall']

Final_Dictionary_Microbes={}
for Microbe_Final_Dict in Precision_Recall_Microbes:
    Final_Dictionary_Microbes[Microbe_Final_Dict]=\
    {'H_mic':{},'Precision':{},'Recall':{}}
    
    Final_Dictionary_Microbes[Microbe_Final_Dict]['H_mic']=\
            Microbes_Dict[int(Microbe_Final_Dict)][-1]
    
    Final_Dictionary_Microbes[Microbe_Final_Dict]['Precision']=\
        Precision_Recall_Microbes[Microbe_Final_Dict]['Precision']
    
    Final_Dictionary_Microbes[Microbe_Final_Dict]['Recall']=\
        Precision_Recall_Microbes[Microbe_Final_Dict]['Recall']
#===========================================================

#1.5. Build vectors to plot
#========================================================================
Entropy_Vector_Patients=[];
Recall_Vector_Patients=[];Precision_Vector_Patients=[]

for Patient in Final_Dictionary:

    Recall_raw=sum(Final_Dictionary[Patient]['Recall'].values())
    
    Recall=\
    Recall_raw/len( Final_Dictionary[Patient]['Recall'].values() )

    Precision_raw=sum(Final_Dictionary[Patient]['Precision'].values())
    
    Precision=\
    Precision_raw/len(Final_Dictionary[Patient]['Precision'].values() )
    
    Entropy_Vector_Patients.append( Final_Dictionary[Patient]['H_pat'] )
    Recall_Vector_Patients.append(Recall)
    Precision_Vector_Patients.append(Precision)

Moscotronico=open('Moscovieyu.txt','w')
for data in zip(Entropy_Vector_Patients, Recall_Vector_Patients, Precision_Vector_Patients):
     print >>Moscotronico, data

Entropy_Vector_Microbes=[];Recall_Vector_Microbes=[];
Precision_Vector_Microbes=[]

for Microbe in Final_Dictionary_Microbes:

    Recall_raw=sum(Final_Dictionary_Microbes[Microbe]['Recall'].values())
    
    Recall=\
    Recall_raw/len( Final_Dictionary_Microbes[Microbe]['Recall'].values() )

    Precision_raw=sum(Final_Dictionary_Microbes[Microbe]['Precision'].values())
    
    Precision=\
    Precision_raw/len(Final_Dictionary_Microbes[Microbe]['Precision'].values())
    
    Entropy_Vector_Microbes.append(Final_Dictionary_Microbes[Microbe]['H_mic'])
    Recall_Vector_Microbes.append(Recall)
    Precision_Vector_Microbes.append(Precision)

#========================================================================
    
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#3. PLOTS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#3.1. 
#============================================================================

#3.1.1. Fontsizes and Parameters
#=============================
size_eje=12
size_ticks=13
size_title=15
size_letter=17
#=============================

fig=plt.figure(figsize=(10,10))
gs = gridspec.GridSpec(2, 2)
gs.update(\
          left=0.1,right=0.95,bottom=0.1,top=0.94,wspace=0.25,hspace=0.3)

Recall_Patients= plt.subplot(gs[0, 0])
plt.scatter(Entropy_Vector_Patients, Recall_Vector_Patients)
plt.title(r'$H_{patients}$ vs Recall in Dataset %s' %Dataset)
plt.xlabel(r'$H_{patients}$')
plt.ylabel('Recall')

Precision_Patients= plt.subplot(gs[0, 1])
plt.scatter(Entropy_Vector_Patients, Precision_Vector_Patients)
plt.title(r'$H_{patients}$ vs Precision in Dataset %s' %Dataset)
plt.xlabel(r'$H_{patients}$')
plt.ylabel('Precision')

Recall_Microbes= plt.subplot(gs[1, 0])
plt.scatter(Entropy_Vector_Microbes, Recall_Vector_Microbes)
plt.title(r'$H_{microbes}$ vs Recall in Dataset %s' %Dataset)
plt.xlabel(r'$H_{microbes}$')
plt.ylabel('Recall')

Precision_Microbes= plt.subplot(gs[1, 1])
plt.scatter(Entropy_Vector_Microbes, Precision_Vector_Microbes)
plt.title(r'$H_{microbes}$ vs Precision in Dataset %s' %Dataset)
plt.xlabel(r'$H_{microbes}$')
plt.ylabel('Precision')



plt.savefig(\
'../../Plots/First_Path/Entropy_vs_Precision_Recall_'+Dataset+'.pdf',dpi=300)
plt.show()
#============================================================================

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
