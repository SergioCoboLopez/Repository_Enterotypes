import pandas as pd
import time
import random as rnd
from random import shuffle
import math
import numpy as np
from sklearn.model_selection import KFold
import pickle

import bisect

rnd.seed(int(1111))

#0. FUNCTIONS
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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

    # Abundances=Catalogue(AdjacencyMatrix)
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
    index=float(100/float(NumberOfLinks) ) 
    
    x=np.arange(index,100+index, index)
#    print x

    Percentiles=[np.percentile(No_Zeros,paso) for paso in x ]
    print Percentiles
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
    #0.3.2. Classify new Adjacency Matrix by percentiles
    #:::::::::::::::::::::::::::::::::::::::::::::::::::
    AdjacencyMatrix=np.zeros( (rows,columns) )
    for row in range(rows):
        for column in range(columns):
            
            for index in range(len(Percentiles)):
#                if NumpySubset[row][column]<=Percentiles[index]:
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

            
#0.3. Classify values by order of magnitude
#==========================================================
def Order_of_Magnitude(NumpyDataSet):
    
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
                TypeOfLink=-int(OrderOfMagnitude)
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
    #::::::::::::::::::::::::::::::::
    NominalLinks=Links.keys()           #Links por su id
    ActualLinks=range(len(Links.keys()))#Links por su numero
    if sum(NominalLinks)!=sum(ActualLinks):
        print "Aqui faltan links, corrijo Matriz de Adyacencia:"
        AdjacencyMatrix=Squasher(AdjacencyMatrix,Links)
        print AdjacencyMatrix
        Links=Catalogue(AdjacencyMatrix)
    #::::::::::::::::::::::::::::::::
    
    return AdjacencyMatrix
#==========================================================

#0.4. Build Adjacency Matrix by Percentiles
#==========================================
def Percentiles(NumpyDataSet,NumberOfLinks):

    AdjacencyMatrix=np.zeros( (rows,columns) )
    
    index=float(100/float(NumberOfLinks) ) 
    print index
    
    x=np.arange(index,100+index, index)


    Percentiles=[np.percentile(NumpyDataSet,paso) for paso in x ]
    print Percentiles
    #Classify by percentiles
    #::::::::::::::::::::::::::::::::::::::::::::
    for row in range(rows):
        for column in range(columns):
            
            for index in range(len(Percentiles)):
                if NumpySubset[row][column]<=Percentiles[index]:
                    AdjacencyMatrix[row][column]=index
                    break
    #::::::::::::::::::::::::::::::::::::::::::::

    #0.3.3. Generamos catalogo de links
    #::::::::::::::::::::::::::::::::
    Links=Catalogue(AdjacencyMatrix)
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
#==========================================


#0.5. Generate list of indices from the Adj.Matrix
#===========================================
def AdjMatrix2List(AdjacencyMatrix):
    List=[]
    for row in range(rows):
        for column in range(columns):
#            List.append( [row, column, AdjacencyMatrix[row,column], NumpySubset[row,column] ] )
            List.append( [row, column, AdjacencyMatrix[row,column], d1[row,column] ] )
            
    return List
#===========================================

#0.6. Generate K-Fold
#===========================================
def GetFolds(ListIndices,K):
    #0.6.1.Shuffle List
    #:::::::::::::::::::
    shuffle(ListIndices)
    #::::::::::::::::::

    #0.6.2. Generate Folds
    #::::::::::::::::::::::::::::::::::::::::::::::::::::
    Block=len(ListIndices)/K
    Folds=[]
    for ind in range(K):
        Folds.append(ListIndices[ind*Block:(ind+1)*Block])
    #::::::::::::::::::::::::::::::::::::::::::::::::::::


    #0.2.3. Distribute remainders
    #:::::::::::::::::::::::::::::::::::::
    sobras=len(ListIndices)-(ind+1)*Block
    leftovers=ListIndices[-sobras:]
    
    for ind in range(sobras):
        Folds[ind].append(leftovers[ind])
    #:::::::::::::::::::::::::::::::::::::
    
    return Folds
#===========================================

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#1. MAIN
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#1.1. Read Dataset, Crop it and define global variables for matrices
#==================================================================
#d1=0.01*(pd.read_csv('/export/home/shared/Projects/Microbiome/Input_Data/toy_data.txt', sep='\t', header=None))
#Subset=d1[ [column for column in range(10)] ][0:10]


Data=(pd.read_csv('../../Input_Data/Healthy_relative_abundance.csv'))
Data = Data.drop('#', 1) #Remove the names of bacteria
DataMatrix=Data.values   #Convert to Matrix
d1=0.01*DataMatrix #Convert to percentages

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
AdjacencyMatrix_Quartiles,Percentiles=\
Percentiles_no_zeros(d1,NumberOfLinks)
print "Percentiles", Percentiles
#:::::::::::::::::::::::::::::::::::::::::::::
#===============================================

#1.3. Choose Adjacency Matrix to generate datasets
#============================================
Types_of_Data=['Categories', 'Quartiles']

#Choose type of Data
#:::::::::::::::::::::::::::::::::::::::::::::
Choice_Data=Types_of_Data[1]
print Choice_Data

if Choice_Data=='Categories':
    AdjacencyMatrix=AdjacencyMatrix_Orders
elif Choice_Data=='Quartiles':
    AdjacencyMatrix=AdjacencyMatrix_Quartiles
#:::::::::::::::::::::::::::::::::::::::::::::


print '\n',AdjacencyMatrix
#============================================

#1.4. Get List of indices of Adjacency Matrix
#==========================================
ListIndices=AdjMatrix2List(AdjacencyMatrix)
#==========================================

'''
#1.5. Cargar datos a pickle
#==========================================================
pickle.dump([AdjacencyMatrix_Orders,AdjacencyMatrix_Quartiles,d1,ListIndices,d1,Percentiles],\
            open('../../Input_Data/cargoDataAnalysis_LargeDataset.p','wb'))
#==========================================================        


#1.6. Generate K Folds
#=============================
K=5
KFolds=GetFolds(ListIndices,K)
for fold in KFolds:
    print len(fold)
#=============================


#1.7. Generate Train+Test
#===============================================================
path='/export/home/shared/Projects/Microbiome/Input_Data/5Fold_LargeData/' + Choice_Data + '/'

with open(path+'Adyacencia'+ Choice_Data + '.txt','w') as file:
    for line in AdjacencyMatrix:
        line=' '.join(map(str, line))
        file.write(line)
        file.write('\n')
        
#1.7.1. Run over folds
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
for fold in range(K):

    #1.7.1.1. Create files and copy Adj Matrix
    #----------------------------------------------
    f0=open(path+ str(fold) + 'TestCheck.txt', 'w')
    f1=open(path+ str(fold) + 'Test.txt', 'w')
    f2=open(path+ str(fold) + 'Train.txt', 'w')
    Train=np.copy(AdjacencyMatrix)
    #----------------------------------------------

    #1.7.1.2. Generate Test, TestCheck and put nans in Adj Matrix
    #--------------------------------------------------
    for line in KFolds[fold]:
        line_no_brackets=' '.join(map(str, line))     #Remove brackets from array
        test_no_brackets=' '.join(map(str, line[0:2]))#Remove brackets from array

        print >>f0, line_no_brackets
        print >>f1, test_no_brackets

        row=line[0];column=line[1]
        Train[row][column]=np.nan
    #--------------------------------------------------
    
    #1.7.1.3. Generate Train Matrices
    #--------------------------------------------------
    for row in Train:
        row_no_brackets=' '.join(map(str, row))#Remove brackets from array
        print >>f2, row_no_brackets
    #--------------------------------------------------
    
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#===============================================================

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

'''
