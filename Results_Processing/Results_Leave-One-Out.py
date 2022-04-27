#8/10/2018. Este script genera resultados de Likelihood
#y (en el futuro) de scores de las pruebas de Leave-One_Out.
#En el futuro, esperamos tambien que escriba los resultados
#directamente a un .csv

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
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#1.1. Parameters for reading data
#========================================================

#1.1.1. EXTERNAL PARAMETERS
#---------------------
Dataset=sys.argv[1]#--
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

#========================================================


#1.2. Generate+write results
#========================================================
for line in Results:

    #1.2.1. Read Input files
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    Patient=line[0];Fold=line[1];Seed=line[2]
    Parameters = str(Patient) + "_Fold_" + str(Fold) +  "_Seed_" + str(Seed)
    
    LogLikelihood_File = read_file(path_Out+ Parameters + '_LogLikelihood.txt')

    Scores = path_Out + Parameters +  '_scores.txt'
    raw_scores=read_file(Scores)
    Check=read_file(path_In + 'P_' + Patient + '_F_'+Fold+ '_TestCheck.txt')
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


    #1.2. Compute LogLikelihood
    #::::::::::::::::::::::::::::::::::::::::::::
    LogLikelihood=[];
    for value in LogLikelihood_File:
        LogLikelihood.append(float(value[0]))
    #::::::::::::::::::::::::::::::::::::::::::::

    #1.3. Compute Accuracy
    #::::::::::::::::::::::::::::::::::::::::::::

    #1.3.1 Clean scores file
    #----------------------------
    scores=clean_data(raw_scores)
    #----------------------------

    #1.3.2. Collect scores in dictionary
    #------------------------------------------
    scores_dict={}
    for line in scores:
        Microbe=line[0];Person=line[1]
        Bet=line[2:];Winner=Bet.index(max(Bet))
        Pareja_Tupla=(Microbe, Person)
        #Asign prediction
        scores_dict[ Pareja_Tupla ]=Winner
    #------------------------------------------

    #1.3.3. Compare scores with test check file
    #------------------------------------------
    Hits=0
    Predictions_dict={}

    for line in Check:
        MicrobeCheck=float(line[0]); PersonCheck=float(line[1])
        Abundance=float(line[2])
        TuplaCheck=(MicrobeCheck,PersonCheck)

    
        LabelResult='P'+str(scores_dict[TuplaCheck])+'_A'+str( int(Abundance) )
        try:
	    Predictions_dict[ LabelResult ]+=1
        except KeyError:
            Predictions_dict[ LabelResult ]=1
            
        if float(line[2])==scores_dict[TuplaCheck]:
            Hits+=1            
    #------------------------------------------
    

    #1.3.4. Accuracy
    #------------------------------------------
    Accuracy=float(Hits)/len(scores)#--->Accuracy Rate
    #------------------------------------------
    #::::::::::::::::::::::::::::::::::::::::::::

    #1.4. Open Dataset
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    df=pd.read_csv('../../Reports/Results_Sheet_Tot.csv',index_col=[0,1,2])
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #1.5. Write results to dataset
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    df.set_value((Dataset,int(Patient) ,'Likelihood'),'Fold'+Fold,\
             LogLikelihood[-1] )

    df.set_value((Dataset,int(Patient),'Seed'),'Fold'+Fold,\
             Seed )

    #1.5.1. Write predicted vs. actual scores in a loop
    #--------------------------------------------------------
    for key in Predictions_dict:
        df.set_value((Dataset,int(Patient),key),'Fold'+Fold,\
                     Predictions_dict[key] )
    #--------------------------------------------------------
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #1.6. Save
    #::::::::::::::::::::::::::::::::::::::::
    df.to_csv('../../Reports/Results_Sheet_Tot.csv')
    #::::::::::::::::::::::::::::::::::::::::


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++

