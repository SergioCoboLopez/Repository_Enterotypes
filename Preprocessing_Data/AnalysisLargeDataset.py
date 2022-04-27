#1/06/2018. En este codigo vamos a ver como estan distribuidas las abundancias del dataset con nombre:
#'Healthy_relative_abundance.csv'.
import pandas as pd
import time
import random as rnd
from random import shuffle
import math
import numpy as np
from sklearn.model_selection import KFold
import pickle
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

#np.set_printoptions(threshold=np.nan)


#0. FUNCTIONS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def abundancias(vector):
    contador=[0,0,0,0,0,0,0,0]
    for abundancy in vector:
        #Contador de 0
        #======================
        if abundancy<10**(-7):
            contador[0]+=1
        #======================

        #Contador de nx10^-7
        #=============================================
        if abundancy<10**(-6) and abundancy>=10**(-7):
            contador[1]+=1
        #=============================================
    
        #Contador de nx10^-6
        #=============================================
        if abundancy<10**(-5) and abundancy>=10**(-6):
            contador[2]+=1
        #=============================================

        #Contador de nx10^-5
        #=============================================
        if abundancy<10**(-4) and abundancy>=10**(-5):
            contador[3]+=1
        #=============================================
        
        #Contador de nx10^-4
        #=============================================
        if abundancy<10**(-3) and abundancy>=10**(-4):
            contador[4]+=1
        #=============================================

        #Contador de nx10^-3
        #=============================================
        if abundancy<10**(-2) and abundancy>=10**(-3):
            contador[5]+=1
        #=============================================

        #Contador de nx10^-2
        #=============================================
        if abundancy<10**(-1) and abundancy>=10**(-2):
            contador[6]+=1
        #=============================================
        
        #Contador de nx10^-1
        #=============================================
        if abundancy<10**(0) and abundancy>=10**(-1):
            contador[7]+=1
        #=============================================
    return contador
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#1. MAIN
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#1.1. Read Dataset, Crop it and define global variables for matrices
#==================================================================
Data=(pd.read_csv('../../Input_Data/Healthy_relative_abundance.csv'))
Data = Data.drop('#', 1) #Remove the names of bacteria
DataMatrix=Data.values   #Convert to Matrix

d1=0.01*DataMatrix #Convert to percentages

print "rows", d1.shape[0], "columns ", d1.shape[1]

global rows;   rows=d1.shape[0]
global columns;columns=d1.shape[1]

length=rows*columns

vector_ratios=np.reshape(d1,length) #Convert matrix into vector
#==================================================================

#1. Clasificamos las abundancias por ordenes de magnitud
#==================================================================
contador1=abundancias(vector_ratios)
print contador1, sum(contador1)

Ratios=[];
for category in contador1:
    Ratios.append( float(category)/float(length) )

print Ratios, sum(Ratios)

x=np.arange(10,110,10)
Percentiles1=[np.percentile(vector_ratios,paso) for paso in x ]
print Percentiles1
#==================================================================

fig=plt.figure(figsize=(15,7))
gs = gridspec.GridSpec(1, 2)
gs.update(left=0.1, right=0.9, wspace=0.3,hspace=0.4)

#Pintar Histograma
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Fijamos los bins
#=========================================================================
bins=(10**(-7),10**(-6),10**(-5),10**(-4),10**(-3),10**(-2),10**(-1),10**0)
#=========================================================================
#Fijamos xticks y posiciones
#================================================
xticks=[tick for tick in bins[:-1] ]
xticks_position=[2.5*tick for tick in bins[:-1] ]
#================================================

#[0,0] Histograma Healthy Subject Data
#================================================
ax00=plt.subplot(gs[0,0])
plt.xscale('log')  #Escala logaritmica
#plt.yscale('log')  #Escala logaritmica
plt.hist(vector_ratios,bins=bins)
plt.xticks(xticks_position,xticks,rotation='-45')
plt.title('Distribution of orders of magnitude of relative abundances (Healthy_relative_abundance)')
plt.ylabel('occurrences')
plt.xlabel('abundance (log scale)')
plt.axhline(y=contador1[0], xmin=0, xmax=1, hold=None, color='k') #Linea para marcar ceros
#================================================

#[0,1] Pintar grafica de percentiles
#================================================
ax10=plt.subplot(gs[0,1])
plt.scatter(x, Percentiles1)
plt.title('Percentile Plot (10%) (Healthy_relative_abundance)')
ax10.set_ylim((1e-7, 1))
plt.ylabel('Percentile value (log scale)')
plt.xlabel('Percentage')
plt.yscale('log')
#================================================ 
plt.show()

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
