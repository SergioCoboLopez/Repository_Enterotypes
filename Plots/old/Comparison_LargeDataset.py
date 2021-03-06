import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib
from pylab import *
from matplotlib import colors
import matplotlib.patches as patches
import seaborn as sns
import pickle
import time
import copy
import operator


#1. Data (results)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#1.1. CATEGORIES
#==========================================================================================================

#1.1.1. Data
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
CategoriesBestLikelihood=\
[0.8284457877, 0.833074661, 0.8346403093, 0.8285547024, 0.8339868213,\
 0.8396503839, 0.832829603, 0.8389015956,0.8394189403,
 0.8398192017, 0.8375837281, 0]#, 0, 0, 0, 0, 0, 0]
        
CategoriesAverage=\
[ 0.8281843925, 0.8306948756, 0.8341828677, 0.8280700321, 0.8329711921,\
  0.8371943582, 0.8330038665, 0.8367641453,0.8381609759, 0.8363039808,
  0.837341393, 0.8398192017, 0.8375837281, 0.8392365082]#, 0, 0, 0, 0, 0, 0]

CategoriesError=\
[0.0001039871, 0.0009443226, 0.0003082256, 0.0003414349, 0.0003444076, 0.0010010248, 0.0001617415,\
 0.0009491153, 0.0008034602, 0.0005162759, 0.0010638066, 0.0004425098, 0.0007442995, 0.0004960664]#, 0, 0, 0, 0, 0, 0]
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#1.1.2. Parameters for plot
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
ymin=0.775
ymax=0.845
yticks_Labels_Categories=[0.78, 0.80, 0.82, 0.84]

xticks_Labels_Categories=\
['K=3,L=5','K=5,L=5', 'K=5,L=10', 'K=10,L=5', 'K=10,L=10', 'K=10,L=20', 'K=20,L=10', 'K=20, L=20',\
 'K=20, L=30', 'K=30, L=20', 'K=30, L=30', 'K=30, L=40', 'K=40, L=30', 'K=40, L=40', 'K=40, L=50',\
 'K=50, L=40', 'K=50, L=50']
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#==========================================================================================================

#1.2. QUARTILES
#==========================================================================================================

#1.2.1. Data
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
QuartilesAverage=\
[0.8195283995, 0.8196509285, 0.8222485433, 0.8185154931, 0.8221287371, 0.8230218374, 0.821355443, 0.8230517889,\
 0.8243288134, 0.8251021075, 0.825929859, 0.8271687633]  
                       #, 0, 0, 0, 0, 0, 0, 0, 0, 0]

QuartilesError=\
[0.0001051922, 0.0002191443, 0.0002292874, 0.0002249292, 0.0002374911, 0.0005673862, 0.0003156764, 0.0006316141,\
 0.0003324407, 0.0005292212, 0.0005511098, 0.0006077912]
          #, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


#1.2.2. Parameters for plot
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
ymin=0.775
ymax=0.845
yticks_Labels_Quartiles=[0.78, 0.80, 0.82, 0.84]

xticks_Labels_Quartiles=\
['K=3,L=5','K=5,L=5', 'K=5,L=10', 'K=10,L=5', 'K=10,L=10', 'K=10,L=15', 'K=15,L=10', 'K=15, L=15',\
 'K=20, L=20',  'K=25, L=25', 'K=30, L=30', 'K=40, L=40' ]
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#==========================================================================================================

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#2. PLOTS
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#2.1. Plotting Parameters
#==============================================

#2.1.1. Ejes: limites y nombres
#:::::::::::::::::::::::
size_eje=15
size_ticks=9
#:::::::::::::::::::::::

#2.1.2.Xticks positions Categories
#::::::::::::::::::::::::::::::::::::::::
xticks_Position_Categories=[];contador=0
for i in range(len(CategoriesAverage)):
    xticks_Position_Categories.append(contador)
    contador+=0.75
#::::::::::::::::::::::::::::::::::::::::

#2.1.3.Xticks positions Quartiles
#::::::::::::::::::::::::::::::::::::::::
xticks_Position_Quartiles=[];contador=0
for i in range(len(QuartilesAverage)):
    xticks_Position_Quartiles.append(contador)
    contador+=0.75
#::::::::::::::::::::::::::::::::::::::::

#==============================================

#2.2. Gridspec parameters
#=============================
fig=plt.figure(figsize=(17,7))
gs = gridspec.GridSpec(1, 2)
#=============================

#2.3. (0,0) PLOT
#=========================================================================================
ax0=plt.subplot( gs[0,0] )

plt.errorbar(xticks_Position_Categories, CategoriesAverage, CategoriesError, marker='o',markersize=8, ls='solid')
plt.hlines(0.7876418947, -0.5, 100, colors='k', linestyles='solid')

ax0.set_title('Categories')
ax0.set_ylabel('Accuracy',fontsize=size_eje)
ax0.set_xlabel('Number of Groups',fontsize=size_eje)

#xticks and yticks
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
xTickMarks_Cat = xticks_Labels_Categories
ax0.set_xticks([pos for pos in xticks_Position_Categories])#Fija la posicion de los xticks
xtickNames = ax0.set_xticklabels(xTickMarks_Cat)
plt.setp(xtickNames, rotation=45, fontsize=size_ticks)

yTicksMarks=yticks_Labels_Categories
ax0.set_yticks(yTicksMarks)
yticknames=ax0.set_yticklabels(yTicksMarks)
plt.yticks(fontsize=size_ticks)
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#Limite ejes
#::::::::::::::::::::::::::::::::::::::
ax0.set_xlim( xticks_Position_Categories[0]-0.5, xticks_Position_Categories[-1]+0.5 )
ax0.set_ylim(ymin,ymax)
#::::::::::::::::::::::::::::::::::::::

#=========================================================================================


#2.4. (0,1) PLOT
#=========================================================================================
ax1=plt.subplot( gs[0,1] )

plt.errorbar(xticks_Position_Quartiles,QuartilesAverage, yerr=QuartilesError, marker='o',markersize=8, ls='solid')
plt.hlines(0.7876418947, -0.5, 100, colors='k', linestyles='solid')

ax1.set_title('Quartiles')

#xticks and yticks
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
xTickMarks_Qrt = xticks_Labels_Quartiles
ax1.set_xticks([pos for pos in xticks_Position_Quartiles] )#Fija la posicion de los xticks
xtickNames = ax1.set_xticklabels(xTickMarks_Qrt)
plt.setp(xtickNames, rotation=45, fontsize=size_ticks)

yTicksMarks=yticks_Labels_Quartiles
ax1.set_yticks(yTicksMarks)
yticknames=ax1.set_yticklabels(yTicksMarks)
plt.yticks(fontsize=size_ticks)
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#Limite ejes
#::::::::::::::::::::::::::::::::::::::
ax1.set_xlim( xticks_Position_Quartiles[0]-0.5, xticks_Position_Quartiles[-1]+0.5 )
ax1.set_ylim(ymin,ymax)
#::::::::::::::::::::::::::::::::::::::

#=========================================================================================

#2.5. Save and show
#================================================
plt.savefig('../../Plots/Comparison.pdf')
plt.show()
#================================================

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
