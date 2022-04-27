import pandas as pd
import time
import random as rnd
from random import shuffle
import math
import numpy as np
from sklearn.model_selection import KFold
import pickle

import bisect

rnd.seed(int(3333))

#0. FUNCTIONS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#0.1. Identify how many links exist in the AdjacencyMatrix
#==========================================================
def Catalogue(AdjacencyMatrix):
    #0.1.1. Catalogo de abundancias
    #::::::::::::::::::::::::::::::::
    Links={}
    for row in AdjacencyMatrix:
        for column in row:
            try:
                Links[column]+=1
            except KeyError:
                Links[column]=1
    #::::::::::::::::::::::::::::::::
    
    return Links
#==========================================================

#0.2. Squash Adjacency Matrix (if types of links are missing)
#==========================================================
def Squasher(AdjacencyMatrix,Links):
    #0.2.1. Create a dictionary that maps old links to new ones
    #(sorting them from smaller to largest)
    #::::::::::::::::::::::::::::::::::::
    Indices={}
    SortedLinks=sorted(Links.keys())
    counterIndices=0
    for link in SortedLinks:
        Indices[link]=counterIndices
        counterIndices+=1
    #::::::::::::::::::::::::::::::::::::

    #0.2.2. Use the dictionary to rebuild Adjacency Matrix
    #::::::::::::::::::::::::::::::::::::
    for Original in Indices.keys():
        AdjacencyMatrix[AdjacencyMatrix==Original]=Indices[Original]
    #::::::::::::::::::::::::::::::::::::

    return AdjacencyMatrix
#==========================================================


#0.3. Classify values by categories
#==========================================================
def Categories(NumpyDataSet):
    
    #0.3.1. Declaramos matriz de adyacencia
    #::::::::::::::::::::::::::::::::::::::
    AdjacencyMatrix=np.zeros( (rows,columns) )
    #::::::::::::::::::::::::::::::::::::::
    
    #0.3.2. Rellenamos matriz de adyacencia
    #::::::::::::::::::::::::::::::::::::::
    for row in range(rows):
        for column in range(columns):
            #---------------------------------------------------------
            try:
                OrderOfMagnitude=\
                        math.floor(math.log10( NumpyDataSet[row][column]))
                
                #A bit
                #_______________________
                if OrderOfMagnitude<=-5:
                    TypeOfLink=1
                #_______________________
                
                #Medium
                #_______________________
                elif OrderOfMagnitude==-4 or OrderOfMagnitude==-3:
                    TypeOfLink=2
                #_______________________

                #Much
                #_______________________
                else:
                    TypeOfLink=3
                #_______________________
                
            except ValueError: #Prevenimos singularidades (log(0))
                TypeOfLink=0
            #---------------------------------------------------------
            AdjacencyMatrix[row][column]=TypeOfLink
    #::::::::::::::::::::::::::::::::::::::

    #0.3.3. Generamos catalogo de links
    #::::::::::::::::::::::::::::::::
    Links=Catalogue(AdjacencyMatrix)
    #::::::::::::::::::::::::::::::::

    #0.3.4. Checkeamos links no consecutivos o similares
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    NominalLinks=Links.keys()             #Links por su id
    ActualLinks=range(len(Links.keys()))  #Links por su numero
    if sum(NominalLinks)!=sum(ActualLinks):
        print "Aqui faltan links, corrijo Matriz de Adyacencia:"
        AdjacencyMatrix=Squasher(AdjacencyMatrix,Links)
        print AdjacencyMatrix
        Links=Catalogue(AdjacencyMatrix)
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
    return AdjacencyMatrix
#==========================================================


#0.4. Generate list of indices from a column of the Adj. Matrix
#==============================================================
def Column2List(Patient):
    List_of_Patient_Indices=[]

    #Generate list of info for given patient: row, column and abundances
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    for row in range(rows):
        List_of_Patient_Indices.append( [row, Patient, \
        AdjacencyMatrix[row,Patient], DataMatrix[row,Patient] ] )
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
    return List_of_Patient_Indices
#=============================================================

#0.5. Generar TestCheck con la informacion del paciente
#=============================================================
def FiveFold_Leave_One_Patient_Out(IndicesPatient):
    
    shuffle(IndicesPatient) #Randomize

    #0.5.1. Generate Folds
    #::::::::::::::::::::::::::::::::::::::::::::::::::::
    Block=len(IndicesPatient)/5
    Folds=[]
    for ind in range(5):
        Folds.append(IndicesPatient[ind*Block:(ind+1)*Block])
    #::::::::::::::::::::::::::::::::::::::::::::::::::::

    #0.5.2. Distribute remainders
    #:::::::::::::::::::::::::::::::::::::
    sobras=len(IndicesPatient)-(ind+1)*Block
    leftovers=IndicesPatient[-sobras:]

    for ind in range(sobras):
        Folds[ind].append(leftovers[ind])
    #:::::::::::::::::::::::::::::::::::::

    
    return Folds
#=============================================================

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



#1. MAIN
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#1.1. Read Dataset, Crop it and define global variables for matrices
#==================================================================
rnd.seed(1) #Generamos semilla

Dataset=['S-8', 'V-10', 'V-22', 'V-23_24', 'V-25']

#---------
Index=0#--
#---------

Data=(pd.read_csv('../../Input_Data/' + Dataset[Index] +'.txt',\
      sep="\t", header=None ))
DataMatrix=Data.values   #Convert to Matrix

global rows;   rows=DataMatrix.shape[0]
global columns;columns=DataMatrix.shape[1]
print columns, rows
#==================================================================

#1.2. Get Adjacency Matrices
#===============================================
AdjacencyMatrix=Categories(DataMatrix)
print '\n',AdjacencyMatrix

#1.3.2. Write Adjacency Matrix to file
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
path='../../Input_Data/Leave_One_Out_' + Dataset[Index] + '/'

with open(path+'Adyacencia_'+ Dataset[Index] + '.txt','w') as file:
    for line in AdjacencyMatrix:
        line=' '.join(map(str, line))
        file.write(line)
        file.write('\n')
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
#==================================================================


#1.4. Pseudo-leave-one-out
#====================================================================

#1.4.1. Randomly choose a sample of N patients
#::::::::::::::::::::::::::::::::::::::::::::::::::
Sample_of_Patients=30
Patients_Chosen=rnd.sample(range(columns - 1), Sample_of_Patients)
HiddenInfo=0.20 #--->Esconder el 20%
TestSize=int( HiddenInfo*rows )
#::::::::::::::::::::::::::::::::::::::::::::::::::

#1.4.2. Generate Cross-Validation Datasets
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
F_PatientRecord=open(path+'Patient_Record_'+str(Dataset[Index]) + '.txt','w')
for Patient in Patients_Chosen:
  #  PatientRecord.append(Patient)
    print >>F_PatientRecord, Patient
    #1.4.2.1. Generate Folds
    #-----------------------------------------------
    #Row(microbe),column(patient),link type and abundance (matrix elements)
    IndicesPatient=Column2List(Patient)
    Folds=FiveFold_Leave_One_Patient_Out(IndicesPatient)
    #-----------------------------------------------

    #1.4.2. Write Train, Test, TestCheck to files
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    for Fold in range(len(Folds)):
#        print Fold

        #1.4.2.1. Create files and copy Adj Matrix
        #----------------------------------------------
        f0=open(path+'P_'+str(Patient)+'_F_'+str(Fold)+'_TestCheck.txt','w')
        f1=open(path+'P_'+str(Patient)+'_F_'+str(Fold)+ '_Test.txt','w')
        f2=open(path+'P_'+str(Patient)+'_F_'+str(Fold)+ '_Train.txt','w')
        Train=np.copy(AdjacencyMatrix)
        #----------------------------------------------

        #1.4.2.2. Generate Test, TestCheck and put nans in Adj Matrix
        #--------------------------------------------------
        for line in Folds[Fold]:
#            print line
            #Remove brackets from array
            line_no_brackets=' '.join(map(str, line))
            test_no_brackets=' '.join(map(str, line[0:2]))

            print >>f0, line_no_brackets
            print >>f1, test_no_brackets

            row=line[0];column=line[1]
            Train[row][column]=np.nan
        #--------------------------------------------------

        #1.4.2.3. Generate Train File
        #--------------------------------------------------
        for row in Train:
            row_no_brackets=' '.join(map(str, row))#Remove brackets from array
            print >>f2, row_no_brackets
        #--------------------------------------------------

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#====================================================================

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



