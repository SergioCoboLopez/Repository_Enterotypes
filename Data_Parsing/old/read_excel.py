#4/06/2018. Por ahora, este codigo lee el archivo 'Healthy_relative_abundance.csv' y lo convierte en algo
#tratable para C++. La segunda parte es construir el 5-fold, pero es muy pronto para eso.
import pandas as pd
import time
import random as rnd
from random import shuffle
import math
import numpy as np
from sklearn.model_selection import KFold
import pickle
np.set_printoptions(threshold=np.nan)

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

    x=np.arange(index,100+index, index)


    Percentiles=[np.percentile(NumpyDataSet,paso) for paso in x ]
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
            List.append( [row, column, AdjacencyMatrix[row,column], NumpySubset[row,column] ] )
            
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
Data=(pd.read_csv('../Input_Data/Healthy_relative_abundance.csv'))
Data = Data.drop('#', 1) #Remove the names of bacteria
DataMatrix=Data.values   #Convert to Matrix

print DataMatrix.shape[0], DataMatrix.shape[1]

    
global rows;   rows=DataMatrix.shape[0]
global columns;columns=DataMatrix.shape[1]

d1=0.01*DataMatrix
print d1
print np.amax(d1)



#1.2. Get Adjacency Matrix
#=====================================
AdjacencyMatrix_Orders=Order_of_Magnitude(NumpySubset)
NumberOfLinks=5
AdjacencyMatrix,Percentiles=Percentiles(NumpySubset,NumberOfLinks)
print '\n',AdjacencyMatrix, '\n', '\n',Percentiles
#=====================================

#1.2. Get List of indices of Adjacency Matrix
#======================================
ListIndices=AdjMatrix2List(AdjacencyMatrix)
#======================================

#1.3. Cargar datos a pickle
#==========================================================
pickle.dump([AdjacencyMatrix,AdjacencyMatrix_Orders,NumpySubset,ListIndices,NumpySubset,Percentiles],\
            open('../Input_Data/cargoDataAnalysis.p','wb'))
#==========================================================
'''
#1.4. Generate K Folds
#=============================
K=5
KFolds=GetFolds(ListIndices,K)
#=============================


#1.4. Generate Train+Test
#===============================================================
path='/export/home/shared/Projects/Microbiome/Input_Data/5Fold_ToyData/'

with open(path+'AdyacenciaPercentiles.txt','w') as file:
    for line in AdjacencyMatrix:
        line=' '.join(map(str, line))
        file.write(line)
        file.write('\n')
        
#1.4.1. Run over folds
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
for fold in range(K):

    #1.4.1.1. Create files and copy Adj Matrix
    #----------------------------------------------
    f0=open(path+ str(fold) + 'TestCheck.txt', 'w')
    f1=open(path+ str(fold) + 'Test.txt', 'w')
    f2=open(path+ str(fold) + 'Train.txt', 'w')
    Train=np.copy(AdjacencyMatrix)
    #----------------------------------------------

    #1.4.1.2. Generate Test, TestCheck and put nans in Adj Matrix
    #--------------------------------------------------
    for line in KFolds[fold]:
        line_no_brackets=' '.join(map(str, line))     #Remove brackets from array
        test_no_brackets=' '.join(map(str, line[0:2]))#Remove brackets from array

        print >>f0, line_no_brackets
        print >>f1, test_no_brackets

        row=line[0];column=line[1]
        Train[row][column]=np.nan
    #--------------------------------------------------
    
    #1.4.1.3. Generate Train Matrices
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
