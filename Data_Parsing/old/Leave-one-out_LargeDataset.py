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
            #------------------------------------------------------------------------
            try:
                OrderOfMagnitude=math.floor(math.log10( NumpyDataSet[row][column]  ))
                #A bit
                #_______________________
                if OrderOfMagnitude<=-4:
                    TypeOfLink=1
                #_______________________
                
                #Medium
                #_______________________
                elif OrderOfMagnitude==-3:
                    TypeOfLink=2
                #_______________________

                #Much
                #_______________________
                else:
                    TypeOfLink=3
                #_______________________
                
            except ValueError: #Prevenimos singularidades (log(0))
                TypeOfLink=0
            #------------------------------------------------------------------------
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


#0.4. Build Adjacency Matrix by Percentiles WITHOUT ZEROS
#========================================================
def Percentiles_no_zeros(NumpyDataSet,NumberOfLinks):

    #0.3.0. Get only non-zero elements
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    No_Zeros_Indices=np.nonzero(d1)
    No_Zeros=[]

    for element in zip( No_Zeros_Indices[0],No_Zeros_Indices[1] ):
        No_Zeros.append(d1[element])
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
    #0.3.1. Get percentiles of non-zero elements
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    index=float( 100/float(NumberOfLinks) ) 
    
    x=np.arange(index,100+index, index)

    Percentiles=[np.percentile(No_Zeros,paso) for paso in x ]
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
    #0.3.2. Classify new Adjacency Matrix by percentiles
    #:::::::::::::::::::::::::::::::::::::::::::::::::::
    AdjacencyMatrix=np.zeros( (rows,columns) )
    for row in range(rows):
        for column in range(columns):
            
            for index in range(len(Percentiles)):
                if NumpyDataSet[row][column]<=Percentiles[index] and NumpyDataSet[row][column]!=0:
                    AdjacencyMatrix[row][column]=index+1
                    break
    #:::::::::::::::::::::::::::::::::::::::::::::::::::

    #0.3.3. Generamos catalogo de links
    #::::::::::::::::::::::::::::::::
    Links=Catalogue(AdjacencyMatrix)
    print Links
    #::::::::::::::::::::::::::::::::

    #0.3.4. Checkeamos links no consecutivos o similares
    #::::::::::::::::::::::::::::::::
    NominalLinks=Links.keys()           #Links por su id
    ActualLinks=range(len(Links.keys()))#Links por su numero
    if sum(NominalLinks)!=sum(ActualLinks):
        print "Aqui faltan links, corrijo Matriz de Adyacencia:"
        AdjacencyMatrix=Squasher(AdjacencyMatrix,Links)
        print AdjacencyMatrix
        Links=Catalogue(AdjacencyMatrix)
    #::::::::::::::::::::::::::::::::
    
    return AdjacencyMatrix,Percentiles
#========================================================

#0.7. Generate list of indices from a column of the Adj. Matrix
#==============================================================
def Column2List(Patient):
    List_of_Patient_Indices=[]

    #Generate a list of information of the given patient: row of microbe, column and abundances
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    for row in range(rows):
        List_of_Patient_Indices.append( [row, Patient, AdjacencyMatrix[row,Patient], d1[row,Patient] ] )
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
    return List_of_Patient_Indices
#=============================================================

#0.8. Generar TestCheck con la informacion del paciente
#=============================================================
def Pseudo_Leave_One_Out(List,TestSize):
    
    shuffle(IndicesPatient) #Randomize

    TestCheck=[]
    for ind in range(TestSize):
        TestCheck.append( IndicesPatient[ind] )
        
    print len(TestCheck), len(IndicesPatient)
    
    return TestCheck
#=============================================================

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#1. MAIN
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#1.1. Read Dataset, Crop it and define global variables for matrices
#==================================================================
Data=(pd.read_csv('/export/home/shared/Projects/Microbiome/Input_Data/Healthy_relative_abundance.csv'))
Data = Data.drop('#', 1) #Remove the names of bacteria
DataMatrix=Data.values   #Convert to Matrix
d1=0.01*DataMatrix       #Convert to percentages

print d1
global rows;   rows=d1.shape[0]
global columns;columns=d1.shape[1]
print columns, rows
#==================================================================

#1.2. Get Adjacency Matrices
#===============================================
#By orders
#::::::::::::::::::::::::::::::::::::
AdjacencyMatrix_Orders=Categories(d1)
#::::::::::::::::::::::::::::::::::::

#By quantiles
#:::::::::::::::::::::::::::::::::::::::::::::
NumberOfLinks=4 #Number of links without zeros
AdjacencyMatrix_Quartiles,Percentiles=Percentiles_no_zeros(d1,NumberOfLinks)
print "Percentiles", Percentiles
#:::::::::::::::::::::::::::::::::::::::::::::
#===============================================

#1.3. Choose Adjacency Matrix to generate datasets
#==================================================================
Types_of_Data=['Categories', 'Quartiles']

#1.3.1. Choose type of Data
#:::::::::::::::::::::::::::::::::::::::::::::
Choice_Data=Types_of_Data[0]
print Choice_Data

if Choice_Data=='Categories':
    AdjacencyMatrix=AdjacencyMatrix_Orders
elif Choice_Data=='Quartiles':
    AdjacencyMatrix=AdjacencyMatrix_Quartiles
#:::::::::::::::::::::::::::::::::::::::::::::

print '\n',AdjacencyMatrix

#1.3.2. Write Adjacency Matrix to file
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
path='/export/home/shared/Projects/Microbiome/Input_Data/Leave-one-out_LargeData/' + Choice_Data + '/'

with open(path+'Adyacencia'+ Choice_Data + '.txt','w') as file:
    for line in AdjacencyMatrix:
        line=' '.join(map(str, line))
        file.write(line)
        file.write('\n')
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
#==================================================================


#1.4. Pseudo-leave-one-out
#============================================

#1.4.1. Randomly choose a sample of N patients
#::::::::::::::::::::::::::::::::::::::::::::::::::
Sample_of_Patients=100
Patients_Chosen=rnd.sample(range(columns - 1), Sample_of_Patients)
print Patients_Chosen

HiddenInfo=0.25 #--->Esconder el 25%
TestSize=int( HiddenInfo*columns )
#::::::::::::::::::::::::::::::::::::::::::::::::::

#1.4.2. Generate Cross-Validation Datasets
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
for Patient in Patients_Chosen:

    #1.4.2.1. Generate TestCheck
    #-----------------------------------------------
    IndicesPatient=Column2List(Patient)   #List of indices for the patient
    # HiddenInfo=0.25 #--->Esconder el 25%
    # TestSize=int( HiddenInfo*len(IndicesPatient) )
    print TestSize
    print Patient
    
    TestCheck=Pseudo_Leave_One_Out(IndicesPatient, TestSize)
    print TestCheck
    #-----------------------------------------------
        
    #1.4.2. Write Train, Test, TestCheck to files
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #1.4.2.1. Create files and copy Adj Matrix
    #----------------------------------------------
    f0=open(path+ str(Patient) + 'TestCheck.txt', 'w')
    f1=open(path+ str(Patient) + 'Test.txt', 'w')
    f2=open(path+ str(Patient) + 'Train.txt', 'w')
    Train=np.copy(AdjacencyMatrix)
    #----------------------------------------------

    #1.4.2.2. Generate Test, TestCheck and put nans in Adj Matrix
    #--------------------------------------------------
    for line in TestCheck:
        line_no_brackets=' '.join(map(str, line))     #Remove brackets from array
        test_no_brackets=' '.join(map(str, line[0:2]))#Remove brackets from array

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
    

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#===============================================================
print columns, rows
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

