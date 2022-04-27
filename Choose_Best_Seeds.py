#11/10/2018. Este codigo selecciona la mejor likelihood para cada
#tirada de 5 seeds.
#Lo hace para bloques de un dataset entero


#import matplotlib.pyplot as plt
import sys
import os
import numpy as np
import pandas as pd

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

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++


#1.LIKELIHOOD AND ACCURACY
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++

#1.1. Read Data and store it
#========================================================

#1.1.1. EXTERNAL PARAMETERS
#----------------------
Dataset=sys.argv[1];#--
#----------------------

#1.1.2. Paths for input files
#::::::::::::::::::::::::::::::::::::::::::::::::
path='../Output_Data/' + 'Leave_One_Out_' + Dataset + '/'

Record_Of_Patients=read_file(\
'../Input_Data/Leave_One_Out_' +Dataset+'/'+'Patient_Record_'+Dataset+'.txt')
#::::::::::::::::::::::::::::::::::::::::::::::::

#1.1.3. Path for output file (DIRECTORIO DE MEJORES LIKELIHOODS)
#::::::::::::::::::::::::::::::::::::::::::::::::
Best_Likelihoods=open(path + 'Best_Likelihoods.txt', 'w')
#::::::::::::::::::::::::::::::::::::::::::::::::

#Read all likelihood files in loop
#::::::::::::::::::::::::::::::::::::::::::::::::
#Run over Patients
#------------------------------
for line in Record_Of_Patients:
#    print line
    Patient=line[0]
#------------------------------
    #Run over Folds
    #------------------------------ 
    for Fold in range(5):
        Batch_5_Likelihoods={}
#------------------------------ 
        #Run over Seeds
        #-------------------------
        for Seed in ['1357', '3571', '5713', '7135', '7777']:
#        for Seed in range(1,11,2):
        
            Parameters =str(Patient)+"_Fold_"+str(Fold)+ "_Seed_" + str(Seed)
            Likelihood_File=read_file(path + Parameters + '_LogLikelihood.txt')

            LogLikelihood=[];
            for value in Likelihood_File:
                LogLikelihood.append(float(value[0]))
                
            Final_Log_Likelihood=LogLikelihood[-1]

            Batch_5_Likelihoods[Final_Log_Likelihood]=Seed
#            print "Likelihood", Final_Log_Likelihood
        #-------------------------            
#        print Batch_5_Likelihoods
        Best_Likelihood=max(Batch_5_Likelihoods.keys())
        print Best_Likelihood
        print Patient, Fold, Batch_5_Likelihoods[Best_Likelihood]
        print >>Best_Likelihoods,\
            Patient, Fold, Batch_5_Likelihoods[Best_Likelihood]
#        print Batch_5_Likelihoods
#::::::::::::::::::::::::::::::::::::::::::::::::
            
#========================================================

#1.3. Plot LogLikelihood
#================================================
x=[index for index in range(len(LogLikelihood))]
plt.plot( LogLikelihood )

pathFig='../Plots/'+Dataset
plt.savefig(pathFig + 'Likelihood_' + 'Patient_'+ Patient + '_Fold_' + Fold + 'Seed_' + Seed + '.pdf', dpi=300)

#plt.show()
#================================================

#1.4. Clean scores file
#============================
#scores=clean_data(raw_scores)
#============================

#1.5. Collect scores in dictionary
#==========================================
# scores_dict={}
# for line in scores:
#     Microbe=line[0];Person=line[1]
#     Bet=line[2:];Winner=Bet.index(max(Bet))
#     Pareja_Tupla=(Microbe, Person)
#     scores_dict[ Pareja_Tupla ]=Winner
#==========================================

#1.4. Compare scores with test check file
#=========================================================
# Hits=0
# for line in Check:
#     MicrobeCheck=float(line[0]); PersonCheck=float(line[1])
#     TuplaCheck=(MicrobeCheck,PersonCheck)
    
#     if float(line[2])==scores_dict[TuplaCheck]:
#         Hits+=1
#=========================================================

# Accuracy=float(Hits)/len(scores)#--->Accuracy Rate
# print Accuracy

# print Hits, len(scores)



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++

#2. WRITE TO DATAFRAME
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++

#2.1. Create multiple indices
#=======================================================
# df=pd.read_csv('../Reports/Results_Sheet.csv',index_col=[0,1,2,3])

# print df.loc[Dataset,int(Patient),int(Fold),'Likelihood'][Seed]
# df.set_value((Dataset,int(Patient),int(Fold),'Likelihood'),Seed,\
#              LogLikelihood[-1] )

# df.set_value((Dataset,int(Patient),int(Fold),'Accuracy'),Seed,\
#              Accuracy )

# df.to_csv('../Reports/Results_Sheet.csv')



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
