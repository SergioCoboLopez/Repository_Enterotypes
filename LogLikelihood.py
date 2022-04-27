#27/02/2018. Este script toma como entrada los datos de
#la LogLikelihood generados por el codigo en c++ del
#MixMembership y hace una grafica de la misma en funcion
#de las iteraciones.

#Argumentos por linea de comandos: Seed K  L

import matplotlib.pyplot as plt
import sys

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


#1.MAIN
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++

#1.1. Read Data and store it
#========================================================
Datasets=['LargeDataset/Categories/', 'LargeDataset/Quartiles/', 'ToyData/']


#1.1.1. EXTERNAL PARAMETERS
#---------------------------------------------------
Fold=0                                          #---
Seed=sys.argv[1];K=sys.argv[2];L=sys.argv[3]    #---
Dataset=sys.argv[4];#Categories, Quartiles, #Toy#---
#---------------------------------------------------


if Dataset=='Categories':
    Dataset='LargeDataset/Categories/'
elif Dataset=='Quartiles':
    Dataset='LargeDataset/Quartiles/'
elif Dataset=='Toy':
    Dataset='ToyData/'
else:
    print "Aqui hay algo mal, muchacho"
    
path='../Output_Data/'+Dataset
Parameters = str(Fold) + "_Seed_" + str(Seed) + "_K_" + str(K) + "_L_" + str(L)
Output=read_file(path + Parameters + '_LogLikelihood.txt')


LogLikelihood=[];
for value in Output:
    LogLikelihood.append(float(value[0]))
#========================================================

print "Likelihood", LogLikelihood[-1]


#1.2. Plot LogLikelihood
#================================================
x=[index for index in range(len(LogLikelihood))]
plt.plot( LogLikelihood )

pathFig='../Plots/'+Dataset
plt.savefig(pathFig + 'Likelihood_' + 'Seed_' + Seed + '_K_' + K + '_L_' + L + '.pdf', dpi=300)

plt.show()
#================================================


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++

