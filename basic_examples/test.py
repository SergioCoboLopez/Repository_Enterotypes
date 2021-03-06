import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import time
import math
import numpy as np
np.set_printoptions(threshold=np.nan)

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


#+++++++++++++++++++++++++++++++++++++++++++++++++++++
d1=pd.read_csv('../../Input_Data/toy_data.txt', sep='\t', header=None)
d2=pd.read_csv('../../Input_Data/toy_data.txt', sep='\t',header=None)
array_raw1=d1.values             #Take values of the dataframe
array_raw2=d2.values             #Take values of the dataframe

array_ratios1=array_raw1/100     #Take ratios
array_ratios2=array_raw2/100     #Take ratios

array_log1=np.log(array_ratios1) #Take logs of ratios
array_log2=np.log(array_ratios2) #Take logs of ratios

length1=451400
length2=length1

vector_ratios1=np.reshape(array_ratios1,length1) #Convert matrix into vector
vector_ratios2=np.reshape(array_ratios2,length2) #Convert matrix into vector
#+++++++++++++++++++++++++++++++++++++++++++++++++++++


#1. Clasificamos las abundancias por ordenes de magnitud
#+++++++++++++++++++++++++++++++++++++++++++++++++++++
contador1=abundancias(vector_ratios1)
print contador1, sum(contador1)

Ratios=[];
for category in contador1:
        Ratios.append( float(category)/float(length1) )

print Ratios, sum(Ratios)

contador2=abundancias(vector_ratios2)

x=np.arange(10,110,10)
Percentiles1=[np.percentile(vector_ratios1,paso) for paso in x ]
Percentiles2=[np.percentile(vector_ratios2,paso) for paso in x ]
#+++++++++++++++++++++++++++++++++++++++++++++++++++++

fig=plt.figure(figsize=(15,12))
gs = gridspec.GridSpec(2, 2)
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

#[0,0] Histograma Toy-Data
#================================================
ax00=plt.subplot(gs[0,0])
plt.xscale('log')  #Escala logaritmica
#plt.yscale('log')  #Escala logaritmica
plt.hist(vector_ratios1,bins=bins)
plt.xticks(xticks_position,xticks,rotation='-45')
plt.title('Distribution of orders of magnitude of relative abundances (Toy Data)')
plt.ylabel('occurrences')
plt.xlabel('abundance (log scale)')
plt.axhline(y=contador1[0], xmin=0, xmax=1, hold=None, color='k') #Linea para marcar ceros
#================================================

#[0,1] Histograma Taxa
#================================================
ax01=plt.subplot(gs[0,1])
plt.xscale('log')  #Escala logaritmica
plt.yscale('log')  #Escala logaritmica
plt.hist(vector_ratios2,bins=bins)
plt.xticks(xticks_position,xticks,rotation='-45')
plt.title('Distribution of orders of magnitude of relative abundances (Taxa)')
plt.ylabel('occurrences (log scale)')
plt.xlabel('abundance (log scale)')
plt.axhline(y=contador2[0], xmin=0, xmax=1, hold=None, color='k') #Linea para marcar ceros
#================================================

#[1,0] Pintar grafica de percentiles
#================================================
ax10=plt.subplot(gs[1,0])
plt.scatter(x, Percentiles1)
plt.title('Percentile Plot (10%) (Toy Data)')
ax10.set_ylim((1e-7, 1))
plt.ylabel('Percentile value (log scale)')
plt.xlabel('Percentage')
plt.yscale('log')
#================================================

#[1,1] Pintar grafica de percentiles
#================================================
ax11=plt.subplot(gs[1,1])
plt.scatter(x, Percentiles2)
plt.title('Percentile Plot (10%) (Taxa)')
plt.ylabel('Percentile value')
plt.xlabel('Percentage')
ax11.set_ylim((1e-7, 1))
#plt.yscale('log')
#================================================
plt.show()

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
