#8/10/2018. Este script genera resultados de Likelihood
#y (en el futuro) de scores de las pruebas de Leave-One_Out.
#En el futuro, esperamos tambien que escriba los resultados
#directamente a un .csv

#Argumentos por linea de comandos: Seed K  L Id_Paciente

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

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++


#1.LIKELIHOOD AND ACCURACY
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++

#1.1. Read Data and store it
#========================================================
#1.1.1. EXTERNAL PARAMETERS
#---------------------------------------------------
# Patient=sys.argv[1];Fold=sys.argv[2]            #---
# Seed=sys.argv[3];Dataset=sys.argv[4];           #---
Dataset=sys.argv[1]
#---------------------------------------------------
    

# Parameters = str(Patient) + "_Fold_" + str(Fold) +  "_Seed_" + str(Seed)
# Output=read_file(path + Parameters + '_LogLikelihood.txt')

path_Out='../Output_Data/' + 'Leave_One_Out_' + Dataset + '/'
path_In= '../Input_Data/' + 'Leave_One_Out_' + Dataset + '/'

Results=read_file(path_In + 'Best_Seeds.txt')

for line in Results:
    print line
    Patient=line[0];Fold=line[1];Seed=line[2]

    Parameters = str(Patient) + "_Fold_" + str(Fold) +  "_Seed_" + str(Seed)
    
    LogLikelihood_File = read_file(path_Out+ Parameters + '_LogLikelihood.txt')

    Scores = path_Out + Parameters +  '_scores.txt'
    raw_scores=read_file(Scores)
    Check=read_file(path_In + 'P_' + Patient + '_F_'+Fold+ '_TestCheck.txt')
    


#========================================================

#1.2. Compute LogLikelihood
        #========================================================
    LogLikelihood=[];
    for value in LogLikelihood_File:
        LogLikelihood.append(float(value[0]))
        print "Likelihood", LogLikelihood[-1]    
        #========================================================

        #1.3. Plot LogLikelihood
        #================================================
    x=[index for index in range(len(LogLikelihood))]
#    plt.plot( LogLikelihood )

 #   pathFig='../Plots/'+Dataset
  #  plt.savefig(pathFig + 'Likelihood_' + 'Patient_'+ Patient + '_Fold_' + Fold + 'Seed_' + Seed + '.pdf', dpi=300)

    #plt.show()
        #================================================

        #1.4. Clean scores file
        #============================
    scores=clean_data(raw_scores)
        #============================

        #1.5. Collect scores in dictionary
        #==========================================
    scores_dict={}
    for line in scores:
        Microbe=line[0];Person=line[1]
        Bet=line[2:];Winner=Bet.index(max(Bet))
        Pareja_Tupla=(Microbe, Person)
        scores_dict[ Pareja_Tupla ]=Winner
        #==========================================

        #1.4. Compare scores with test check file
        #=========================================================
    Hits=0
    for line in Check:
        MicrobeCheck=float(line[0]); PersonCheck=float(line[1])
        TuplaCheck=(MicrobeCheck,PersonCheck)
    
        if float(line[2])==scores_dict[TuplaCheck]:
            Hits+=1
        #=========================================================

    Accuracy=float(Hits)/len(scores)#--->Accuracy Rate
    print Accuracy
    
    print Hits, len(scores)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++

#2. WRITE TO DATAFRAME
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++

#2.1. Create multiple indices
#=======================================================
    df=pd.read_csv('../Reports/Results_Sheet.csv',index_col=[0,1,2])

    print Seed
    print df.loc[Dataset,int(Patient),'Likelihood']['Fold'+str(Fold)]
    
    df.set_value((Dataset,int(Patient) ,'Likelihood'),'Fold'+Fold,\
             LogLikelihood[-1] )

    df.set_value((Dataset,int(Patient),'Seed'),'Fold'+Fold,\
             Seed )

    df.set_value((Dataset,int(Patient),'Accuracy'),'Fold'+Fold,\
             Accuracy )

    df.to_csv('../Reports/Results_Sheet.csv')



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
