#11/10/2018. Este codigo selecciona la mejor likelihood para cada tirada de
#5 seeds. Lo hace para bloques de un dataset entero

import matplotlib.pyplot as plt
import sys
import os
import numpy as np
import pandas as pd
import sys
sys.setrecursionlimit(5000)

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


#1. MAIN
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
path='../../Output_Data/' + 'Leave_One_Out_' + Dataset + '/'

Record_Of_Patients=read_file(\
'../../Input_Data/Leave_One_Out_'+Dataset+'/'+'Patient_Record_'+Dataset+'.txt')
#::::::::::::::::::::::::::::::::::::::::::::::::

#1.1.3. Path for output file (DIRECTORIO DE MEJORES LIKELIHOODS)
#::::::::::::::::::::::::::::::::::::::::::::::::
path_destin='../../Input_Data/' + 'Leave_One_Out_' + Dataset + '/'
Best_Likelihoods=open(path_destin + 'Best_Seeds.txt', 'w')

pathFig='../../Plots/First_Path/'+Dataset + '/'
#::::::::::::::::::::::::::::::::::::::::::::::::

#Read all likelihood files in loops
#::::::::::::::::::::::::::::::::::::::::::::::::

#Run over Patients
#------------------------------
for line in Record_Of_Patients:
    Patient=line[0]
#------------------------------
    #Run over Folds
    #------------------------------ 
    for Fold in range(5):
        Batch_5_Likelihoods={}
#------------------------------ 
        #Run over Seeds
        #-------------------------
        for Seed in range(1,11,2):
        
            Parameters =str(Patient)+"_Fold_"+str(Fold)+ "_Seed_" + str(Seed)
            Likelihood_File=read_file(path + Parameters + '_LogLikelihood.txt')

            LogLikelihood=[];
            for value in Likelihood_File:
                LogLikelihood.append(float(value[0]))
                
            Final_Log_Likelihood=LogLikelihood[-1]

            #Diccionario Likelihood-Semilla
            #.............................................
            Batch_5_Likelihoods[Final_Log_Likelihood]=Seed
            #.............................................

        #-------------------------

        #Choose best Likelihood and print data to file
        #--------------------------------------------------------
        Best_Likelihood=max(Batch_5_Likelihoods.keys())
#        print Batch_5_Likelihoods[Best_Likelihood], Patient

        print Patient

        print >>Best_Likelihoods,\
            Patient, Fold, Batch_5_Likelihoods[Best_Likelihood]
        #--------------------------------------------------------
        
        #Plot Best Likelihoods
        #--------------------------------------------------------
        Best_Likelihood_File=read_file(path + str(Patient) + "_Fold_" + \
        str(Fold) + "_Seed_" + str(Batch_5_Likelihoods[Best_Likelihood]) +\
        '_LogLikelihood.txt')

        BestLogLikelihood=[];
        for value in Best_Likelihood_File:
            BestLogLikelihood.append(float(value[0]))

        # fig=plt.figure()
        # plt.plot( BestLogLikelihood)
        # plt.savefig(pathFig+'Likelihood_'+'Patient_'+Patient+'_Fold_'+str(Fold)        +'Seed_'+str(Batch_5_Likelihoods[Best_Likelihood])+'.pdf', dpi=300)
        #--------------------------------------------------------
        
#::::::::::::::::::::::::::::::::::::::::::::::::
            
#========================================================

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
