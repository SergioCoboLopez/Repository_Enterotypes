#27/11/2018. Vamos a reducir los scores de 5 a 3 links para generar
#una comparativa "justa" entre los dos escenarios.

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
    
    Scores = path_Out + Parameters +  '_scores.txt'
    raw_scores=read_file(Scores)
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    print '\n', Patient, Fold

    NewScores=\
    open(path_Out + 'ThreeLinkReduced_'+ Parameters + '_scores.txt','w')

    for pat_microbe in raw_scores:
        
        if len(pat_microbe)>0:
        
            Microbe=pat_microbe[0];Patient=pat_microbe[1];
            Five_Ratings=[float(Rating) for Rating in pat_microbe[2:]]
#            print Five_Ratings

            Three_Ratings=[Microbe, Patient, Five_Ratings[0]+Five_Ratings[1], \
                           Five_Ratings[2]+Five_Ratings[3], Five_Ratings[4]]
            
            Three_Ratings_no_brackets=' '.join(map(str, Three_Ratings))
                        

            print >>NewScores, Three_Ratings_no_brackets
            print Three_Ratings

        else:
           print >>NewScores, " "
