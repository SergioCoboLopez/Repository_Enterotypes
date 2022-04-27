#28/09/2018. Dibujamos una confusion matrix para ver en que links acierta y
#se equivoca el MMSBM. En principio, miramos lo que ocurre para K=10, L=20.

import sys
import numpy as np
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

#0. FUNCTIONS
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#0.1. Read File
#======================================================
def read_file(inFileName):
    lines = open(inFileName).readlines()
    d = [line.strip().split() for line in lines]
    return d
#======================================================

#0.1. Clean File
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

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#1. MAIN
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#1.1. Read and prepare files 
#====================================================================
K=str(10)
L=str(20)
Fold="0"
Seed= sys.argv[1]


path_raw_scores='../../Output_Data/' + 'LargeDataset/'  +'Categories/'
path_test_check='../../Input_Data/' + '5Fold_LargeData'+'/Categories/'

raw_scores_file= Fold+'_Seed_' + Seed + '_K_' + K +'L_' + L + '_scores.txt'

raw_scores=read_file(path_raw_scores + raw_scores_file)
Test_Check=read_file(path_test_check + Fold + 'TestCheck.txt')

#====================================================================


#1.2. Clean scores file
#============================
scores=clean_data(raw_scores)
#============================


#1.3. Assign scores in dictionary (link with highest score)
#==========================================
scores_assigned_dict={}
for line in scores:
    Microbe_Scores=line[0];Person_Scores=line[1]
    Bet=line[2:];Winner=Bet.index(max(Bet))
    Tupla_Scores=(Microbe_Scores, Person_Scores)
    scores_assigned_dict[ Tupla_Scores ]=Winner
#==========================================

#print scores_dict


True_Links={}
Guessed_Links={}

Confusion_Matrix_raw=np.zeros((4,4))

for line in Test_Check:
    Microbe=float(line[0]);Patient=float(line[1]);
    Tupla_Test=(Microbe, Patient)
    Abundance=float(line[2])

    #Build dictionary of True Links (in TestCheck file)
    #----------------------------
    try:
        True_Links[Abundance]+=1
    except KeyError:
        True_Links[Abundance]=1
    #----------------------------
    
    Confusion_Matrix_raw[ scores_assigned_dict[Tupla_Test] , int(Abundance)]+=1

    try:
        Guessed_Links[ scores_assigned_dict[Tupla_Test] ]+=1
    except KeyError:
            Guessed_Links[ scores_assigned_dict[Tupla_Test]  ]=1
    



print "Actual Links", True_Links, sum(True_Links.values())
print "Predicted Links", Guessed_Links, sum(Guessed_Links.values())

for type_of_link in True_Links:
    print float(True_Links[type_of_link])/sum(True_Links.values())
    print float(Guessed_Links[type_of_link])/sum(True_Links.values()), '\n'


Confusion_Matrix_norm_cols=np.zeros((4,4))
Confusion_Matrix_norm_rows=np.zeros((4,4))

for row in range(len(Confusion_Matrix_raw)):
    for column in range(len(Confusion_Matrix_raw)):
       Confusion_Matrix_norm_cols[row,column]\
           =Confusion_Matrix_raw[row,column]/True_Links[column]

       Confusion_Matrix_norm_rows[row,column]\
           =Confusion_Matrix_raw[row,column]/Guessed_Links[row]

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#PLOTS
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#2.1. Plotting Parameters
#====================================

#2.1.1. Ejes: limites y nombres
#:::::::::::::::::::::::
size_eje=15
size_ticks=9
#:::::::::::::::::::::::

#====================================

# fig=plt.figure(figsize=(10,8))
# sns.heatmap(Confusion_Matrix_raw,annot=True)
# plt.show()


#2.2. Gridspec parameters
#=============================
fig=plt.figure(figsize=(15,4))
gs = gridspec.GridSpec(1, 3)
gs.update(left=0.05, right=0.99, wspace=0.3,hspace=0.5, bottom=0.12, top=0.91)
#============================= 

#2.3. (0,0) PLOT
#==================================================
Confusion0=plt.subplot(gs[0,0])
sns.heatmap(Confusion_Matrix_norm_cols,annot=True,vmin=0,vmax=1)
Confusion0.set_title('Normalized over actual values', fontsize=size_eje)
Confusion0.set_ylabel('Predicted', fontsize=size_eje)
Confusion0.set_xlabel('Actual', fontsize=size_eje)
#==================================================

#2.3. (0,1) PLOT
#==================================================
Confusion1=plt.subplot(gs[0,1])
sns.heatmap(Confusion_Matrix_norm_rows,annot=True,vmin=0,vmax=1)
Confusion1.set_title('Normalized over predicted values', fontsize=size_eje)
Confusion1.set_ylabel('Predicted', fontsize=size_eje)
Confusion1.set_xlabel('Actual', fontsize=size_eje)
#==================================================
sns.set_style("white")
Calibration=plt.subplot(gs[0,2])
plt.setp(Calibration.get_yticklabels(), visible=False)
plt.setp(Calibration.get_xticklabels(), visible=False)

Calibration.annotate('Distribution of actual Links:',fontweight='bold', xy=(0.01,0.95), xytext=(0.01, 0.95) )

for type_of_link in True_Links:
    Calibration.annotate( str(type_of_link)+": "+str(True_Links[type_of_link]), xy=(0.01 + 0.25*type_of_link , 0.9), xytext=(0.01 + 0.25*type_of_link, 0.9) )

Calibration.annotate('CALIBRATION',fontweight='bold', xy=(0.01,0.80), xytext=(0.01, 0.80) )
    
Calibration.annotate('Actual Links:' , xy=(0.01,0.7), xytext=(0.01, 0.7) )

for type_of_link in True_Links:
        Calibration.annotate( str(type_of_link) +":  "+ str(float(True_Links[type_of_link])/sum(True_Links.values())), xy=(0.01, 0.6 - 0.1*type_of_link), xytext=(0.01, 0.6 - 0.1*type_of_link) )

Calibration.annotate('Predicted Links:' , xy=(0.51,0.8), xytext=(0.55, 0.7) )        
for type_of_link in True_Links:
        Calibration.annotate( str(float(Guessed_Links[type_of_link])/sum(Guessed_Links.values())), xy=(0.55, 0.6 - 0.1*type_of_link), xytext=(0.55, 0.6 - 0.1*type_of_link) )


plt.savefig('../../Plots/Confusion_Matrices_K_'+K+'_L_'+L+'.pdf')

plt.show()
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

