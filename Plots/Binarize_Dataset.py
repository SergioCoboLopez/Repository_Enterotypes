#1/02/2019. Convertimos las abundancias en 0 o 1 y analizamos luego la nestedness del sistema

import pandas as pd
import time
import random as rnd
from random import shuffle
import math
import numpy as np
from sklearn.model_selection import KFold
import pickle
import sys
import operator
import seaborn as sns
import matplotlib.pyplot as plt
import bisect
import matplotlib.gridspec as gridspec

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

#0.3. Classify values by categories: NOVEMBER. 18, SuperMed Link. 1 + 2 
#======================================================================
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
                OrderOfMagnitude=math.floor(math.log10( NumpyDataSet[row][column]))

                print(OrderOfMagnitude)
                #Zero(negligible)
                #_______________________
                if NumpyDataSet[row][column]< 5e-6:
#                if OrderOfMagnitude<-5:
#                if OrderOfMagnitude<-4:
                    TypeOfLink=0
                #_______________________
                
                #Non-zero
                #_______________________
                else:
                    TypeOfLink=1
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
#======================================================================

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#1. MAIN
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#1.1. Read Dataset, Crop it and define global variables for matrices
#==================================================================

#1.1.1. EXTERNAL PARAMETERS
#---------------------
Dataset=sys.argv[1]#--
#---------------------

Data=(pd.read_csv('../../Input_Data/' + Dataset + '.txt',\
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
path='../../Input_Data/Binarized_Datasets/'

with open(path+'Binarized_'+ Dataset + '_Threshold_5e-6.txt','w') as file:
    for line in AdjacencyMatrix:
#        print line
        line=' '.join(map(str, line))
        file.write(line)
        file.write('\n')
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::    
#==================================================================

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
