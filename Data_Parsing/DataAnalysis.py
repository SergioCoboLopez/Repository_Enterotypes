import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pickle

#0.FUNCTIONS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#0.1. Read Data
#===============================================
def read_file(inFileName):
    lines = open(inFileName).readlines()
    d = [line.strip().split() for line in lines]
    return d
#===============================================


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#1.MAIN
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    

#1.0.Cargar Datos
#===============================================================
Data=pickle.load(open('../../Input_Data/cargoDataAnalysis.p','rb'))
AdjacencyMatrix=Data[0];AdjacencyMatrix_Orders=Data[1];AdjacencyMatrix_raw=Data[2];
ListIndices=Data[3];Abundances=Data[4];Percentiles=Data[5];
#===============================================================


#1.1. Calcular frecuencia de links
#==============================
NumLinks=len(Percentiles)
Frecuencia=NumLinks*[0]
for index in ListIndices:
        Frecuencia[int(index[2])]+=1
#===============================


#1.2. Figuras
#==================================================================
fig=plt.figure(figsize=(12,12) )
gs = gridspec.GridSpec(2,2)
gs.update(left=0.1, right=0.9, wspace=0.75,hspace=0.5)


#[0,0] Clustermap Percentiles
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::
ax00=plt.subplot(gs[0,0])
sns.heatmap(AdjacencyMatrix,vmin=0,vmax=2,linewidths=.05)
ax00.set_xlabel('People');ax00.set_ylabel('Microbes')
ax00.set_title('Percentile Heatmap')
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#[0,1] Clustermap Ordenes de Magnitud
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::
ax01=plt.subplot(gs[0,1])
sns.heatmap(AdjacencyMatrix_Orders,linewidths=.05)
ax01.set_xlabel('People');ax01.set_ylabel('Microbes')
ax01.set_title('Heatmap of Magnitude Orders')
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#[1,0] Pintar grafica de percentiles
#::::::::::::::::::::::::::::::::::::::::::::::

#Definir indices y vector de percentiles
#-----------------------------------
index=float(100/float(NumLinks) )
x=np.arange(index,100+index, index)
#-----------------------------------

#Dibujamos la grafica
#-----------------------------------
ax10=plt.subplot(gs[1,0])
plt.scatter(x, Percentiles)
plt.title('Percentile Plot (33.33%) (10x10 Subset)')
ax10.set_ylim((1e-7, 1))
plt.ylabel('Percentile value (log scale)')
plt.xlabel('Percentage')
plt.yscale('log')
#-----------------------------------
#::::::::::::::::::::::::::::::::::::::::::::::

#[1,1] Histograma Ordenes de Magnitud
#::::::::::::::::::::::::::::::::::::::::::::::

#Fijamos los bins
#--------------------------------------------------------------------------
bins=(10**(-7),10**(-6),10**(-5),10**(-4),10**(-3),10**(-2),10**(-1),10**0)
VectorAbundances=np.reshape(AdjacencyMatrix_raw,len(AdjacencyMatrix_raw)**2)
zerocounter=0;
for abundance in VectorAbundances:
    if abundance==0:
        zerocounter+=1
#--------------------------------------------------------------------------

#Fijamos xticks y posiciones
#------------------------------------------------
xticks=[tick for tick in bins[:-1] ]
xticks_position=[2.5*tick for tick in bins[:-1] ]
#------------------------------------------------

#Dibujamos la grafica
#------------------------------------------------
ax11=plt.subplot(gs[1,1])
plt.xscale('log')  #Escala logaritmica
#plt.yscale('log')  #Escala logaritmica
plt.hist(VectorAbundances,bins=bins)
plt.xticks(xticks_position,xticks,rotation='-45')
plt.title('Distribution of orders of magnitude')
plt.ylabel('occurrences')
plt.xlabel('abundance (log scale)')
plt.axhline(y=zerocounter, xmin=0, xmax=1, hold=None, color='k') #Linea para marcar ceros
plt.axvline(x=Percentiles[0], color='r')
plt.axvline(x=Percentiles[1], color='r')
plt.axvline(x=Percentiles[2], color='r')
#------------------------------------------------

#::::::::::::::::::::::::::::::::::::::::::::::

print Percentiles
plt.show()



#==================================================================

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
