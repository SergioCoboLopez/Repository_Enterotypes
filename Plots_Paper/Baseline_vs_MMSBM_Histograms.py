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

sns.set_style('ticks') #cambiar el estilo de las graficas


#PLOTS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#fig=plt.figure(figsize=(35,7))
fig=plt.figure(figsize=(17.5,3.5))

gs = gridspec.GridSpec(1, 5)
gs.update(left=0.06, right=0.96, wspace=0.4,hspace=0.2, bottom=0.15, top=0.88)

#Fontsizes
#:::::::::::::
size_eje=16
size_ticks=16
size_title=18
size_letter=23
size_legend=10.1
#:::::::::::::

#2. Localizacion y ancho de las barras
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::
Bar_buffer=0.05 #Distancia al eje y
ind = [0+Bar_buffer, 0.5 + Bar_buffer , 1 + Bar_buffer] #Posicion de las barras
width = 0.13          #Ancho de las barras
indBaseline=[posB - 0.5*width for posB in ind]
indMMSBM =  [posM + 0.5*width for posM in ind]
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# ticks_posx0=[0.5,2.5,4.5,6.5,8.5,10.5]
# xticks0=[5,7,9,11,13,15]
ticks_posy0=[10.5,8.5,6.5,4.5,2.5,0.5]
yticks0=[0,2,4,6,8,10]
xTickMarks = ['Accuracy    ','  Precision', '     Recall']
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#Data
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#S-8
#--------------------------------------------------------------------
Baseline_S8=[0.792, 0.758, 0.697]
MMSBM_S8=[0.817, 0.798, 0.727]
#Error_Baseline_S8=[0.00243602712299 ,0.0077012991 ,0.00265328900475]
#Error_MMSBM_S8=[0.00243602712299 ,0.0077012991 ,0.00265328900475]
#--------------------------------------------------------------------

#V-10
#--------------------------------------------------------------------
Baseline_V10=[0.7919100467289719, 0.7578187650360866, 0.6971597196606418]
MMSBM_V10=[0.7936369406537607, 0.7776645608314653 , 0.7164851671047691]
#Error_Baseline_V10=[0.00243602712299 ,0.0077012991 ,0.00265328900475]
#Error_MMSBM_V10=[0.00243602712299 ,0.0077012991 ,0.00265328900475]
#--------------------------------------------------------------------

#V-22
#--------------------------------------------------------------------
Baseline_V22=[0.803141679184378, 0.7811390837804375, 0.7295609605673977]
MMSBM_V22=[0.8125858928057784, 0.7948069185781244, 0.7386192807308329]
#Error_Baseline_V22=[0.00243602712299 ,0.0077012991 ,0.00265328900475]
#Error_MMSBM_V22=[0.00243602712299 ,0.0077012991 ,0.00265328900475]
#--------------------------------------------------------------------

#V-23-24
#--------------------------------------------------------------------
Baseline_V23_24=[0.7564131928538709, 0.6822890971879625, 0.5924815251151333]
MMSBM_V23_24=[0.7944342647732479, 0.7477432296890673, 0.6387490628681589]
#Error_Baseline_V23_24=[0.00243602712299 ,0.0077012991 ,0.00265328900475]
#Error_MMSBM_V23_24=[0.00243602712299 ,0.0077012991 ,0.00265328900475]
#--------------------------------------------------------------------

#V-25
#--------------------------------------------------------------------
Baseline_V25=[0.5, 0.5, 0.5]
MMSBM_V25=[0.808253114382786, 0.7954888712241653, 0.7396478131524973]
#Error_Baseline_V25=[0.00243602712299 ,0.0077012991 ,0.00265328900475]
#Error_MMSBM_V25=[0.00243602712299 ,0.0077012991 ,0.00265328900475]
#--------------------------------------------------------------------

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


# PLOT (0,0)
#=============================================================================================
ax = plt.subplot(gs[0,0])
sns.despine(offset=0, ax=ax)#Poner el fondo de la grafica en blanco

#1. Definir colores
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::
cmaptest=matplotlib.cm.get_cmap("BrBG") #Selecciona el heatmap que quieres trocear                           
color1=cmaptest(0.95);color2=cmaptest(0.05)
cmap1=matplotlib.cm.get_cmap("RdBu_r") #Selecciona el heatmap que quieres trocear
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::

#Plot
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Baseline_Bars = ax.bar(indBaseline, Baseline_S8, width, color=[cmap1(0.9),cmap1(0.9),cmap1(0.9)],
                            yerr=0, error_kw=dict(elinewidth=2,ecolor='black'), label="Baseline")

MMSBM_Bars = ax.bar(indMMSBM, MMSBM_S8, width, color=[cmap1(0.1),cmap1(0.1),cmap1(0.1)],
                            yerr=0, error_kw=dict(elinewidth=2,ecolor='black'), label="MMSBM")
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#Ejes: limites y nombres
#::::::::::::::::::::::::::::::::::::::::::::::
max_y_0=0.84
ax.set_xlim(-width,ind[-1]+width+0.5*width)
ax.set_ylim(0.65,max_y_0)
ax.set_ylabel('Performance',fontsize=size_eje)
#::::::::::::::::::::::::::::::::::::::::::::::

#Ticks
#::::::::::::::::::::::::::::::::::::::::::::::::
ax.set_xticks( [pos for pos in ind] )#Fija la posicion de los xticks con una list comprehension
xtickNames = ax.set_xticklabels(xTickMarks)
plt.setp(xtickNames, rotation=0, fontsize=size_ticks)

#yTicksMarks=[ 0.67, 0.69, 0.71, 0.73, 0.75, 0.77,0.79, 0.81, 0.83]
yTicksMarks=[ 0.67, 0.71, 0.75, 0.79, 0.83]
ax.set_yticks(yTicksMarks)
yticknames=ax.set_yticklabels(yTicksMarks)
plt.yticks(fontsize=size_ticks)
ax.legend(prop={'size': size_legend})
#plot.legend(loc=2, prop={'size': 6})
#::::::::::::::::::::::::::::::::::::::::::::::::                                                                 

#Letra para la caption
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
ax.text(-0.5, max_y_0 + 0.01, 'a', fontsize=size_letter,fontweight='bold')
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#=============================================================================================


# PLOT (1,0)
#=============================================================================================
ax = plt.subplot(gs[0,1])
sns.despine(offset=0, ax=ax)#Poner el fondo de la grafica en blanco


#Plot
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Baseline_Bars = ax.bar(indBaseline, Baseline_V10, width, color=[cmap1(0.9),cmap1(0.9),cmap1(0.9)],
                            yerr=0, error_kw=dict(elinewidth=2,ecolor='black'), label="Baseline")

MMSBM_Bars = ax.bar(indMMSBM, MMSBM_V10, width, color=[cmap1(0.1),cmap1(0.1),cmap1(0.1)],
                            yerr=0, error_kw=dict(elinewidth=2,ecolor='black'), label="MMSBM")
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#Ejes: limites y nombres
#::::::::::::::::::::::::::::::::::::::::::::::
max_y_1=0.84
ax.set_xlim(-width,ind[-1]+width+0.5*width)
ax.set_ylim(0.65, max_y_1)
#ax.set_ylabel('Performance',fontsize=size_eje)
#::::::::::::::::::::::::::::::::::::::::::::::

#Ticks
#::::::::::::::::::::::::::::::::::::::::::::::::
#xTickMarks = ['Accuracy  ','  Precision', '  Recall']
#ax.set_xticks([pos + 0.5*width for pos in ind] )#Fija la posicion de los xticks con una list comprehension
ax.set_xticks( [pos for pos in ind] )#Fija la posicion de los xticks con una list comprehension
xtickNames = ax.set_xticklabels(xTickMarks)
plt.setp(xtickNames, rotation=0, fontsize=size_ticks)

#yTicksMarks=[ 0.67, 0.69, 0.71, 0.73, 0.75, 0.77,0.79, 0.81, 0.83]
yTicksMarks=[ 0.67, 0.71, 0.75, 0.79, 0.83]
ax.set_yticks(yTicksMarks)
yticknames=ax.set_yticklabels(yTicksMarks)
plt.yticks(fontsize=size_ticks)
ax.legend(prop={'size': size_legend})
#::::::::::::::::::::::::::::::::::::::::::::::::                                                                 

#Letra para la caption
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
ax.text(-0.5, max_y_1 + 0.01, 'b', fontsize=size_letter,fontweight='bold')
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#=============================================================================================


# PLOT V-22 (0,2)
#=============================================================================================
ax = plt.subplot(gs[0,2])
sns.despine(offset=0, ax=ax)#Poner el fondo de la grafica en blanco

#Plot
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Baseline_Bars = ax.bar(indBaseline, Baseline_V22, width, color=[cmap1(0.9),cmap1(0.9),cmap1(0.9)],
                            yerr=0, error_kw=dict(elinewidth=2,ecolor='black'), label="Baseline")

MMSBM_Bars = ax.bar(indMMSBM, MMSBM_V22, width, color=[cmap1(0.1),cmap1(0.1),cmap1(0.1)],
                            yerr=0, error_kw=dict(elinewidth=2,ecolor='black'), label="MMSBM")
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#Ejes: limites y nombres
#::::::::::::::::::::::::::::::::::::::::::::::
max_y_2=0.84
ax.set_xlim(-width,ind[-1]+width+0.5*width)
ax.set_ylim(0.69,max_y_2)
#ax.set_ylabel('Performance',fontsize=size_eje)
#::::::::::::::::::::::::::::::::::::::::::::::

#Ticks
#::::::::::::::::::::::::::::::::::::::::::::::::
#xTickMarks = ['Accuracy  ','  Precision', '  Recall']
#ax.set_xticks([pos + 0.5*width for pos in ind] )#Fija la posicion de los xticks con una list comprehension
ax.set_xticks( [pos for pos in ind] )#Fija la posicion de los xticks con una list comprehension
xtickNames = ax.set_xticklabels(xTickMarks)
plt.setp(xtickNames, rotation=0, fontsize=size_ticks)

#yTicksMarks=[ 0.67, 0.69, 0.71, 0.73, 0.75, 0.77,0.79, 0.81, 0.83]
yTicksMarks=[ 0.71, 0.75, 0.79, 0.83]
ax.set_yticks(yTicksMarks)
yticknames=ax.set_yticklabels(yTicksMarks)
plt.yticks(fontsize=size_ticks)
ax.legend(prop={'size': size_legend})

#::::::::::::::::::::::::::::::::::::::::::::::::                                                                 

#Letra para la caption
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
ax.text(-0.5, max_y_2 + 0.01, 'c', fontsize=size_letter,fontweight='bold')
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#=============================================================================================


# PLOT V-23_24 (0,3)
#=============================================================================================
ax = plt.subplot(gs[0,3])
sns.despine(offset=0, ax=ax)#Poner el fondo de la grafica en blanco

#Plot
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Baseline_Bars = ax.bar(indBaseline, Baseline_V23_24, width, color=[cmap1(0.9),cmap1(0.9),cmap1(0.9)],
                            yerr=0, error_kw=dict(elinewidth=2,ecolor='black'), label="Baseline")

MMSBM_Bars = ax.bar(indMMSBM, MMSBM_V23_24, width, color=[cmap1(0.1),cmap1(0.1),cmap1(0.1)],
                            yerr=0, error_kw=dict(elinewidth=2,ecolor='black'), label="MMSBM")
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#Ejes: limites y nombres
#::::::::::::::::::::::::::::::::::::::::::::::
max_y_3=0.81
ax.set_xlim(-width,ind[-1]+width+0.5*width)
ax.set_ylim(0.53,max_y_3)
#ax.set_ylabel('Performance',fontsize=size_eje)
#::::::::::::::::::::::::::::::::::::::::::::::

#Ticks
#::::::::::::::::::::::::::::::::::::::::::::::::
#xTickMarks = ['Accuracy  ','  Precision', '  Recall']
#ax.set_xticks([pos + 0.5*width for pos in ind] )#Fija la posicion de los xticks con una list comprehension
ax.set_xticks( [pos for pos in ind] )#Fija la posicion de los xticks con una list comprehension
xtickNames = ax.set_xticklabels(xTickMarks)
plt.setp(xtickNames, rotation=0, fontsize=size_ticks)

#yTicksMarks=[ 0.67, 0.69, 0.71, 0.73, 0.75, 0.77,0.79, 0.81, 0.83]
yTicksMarks=[0.55, 0.63,  0.71, 0.79]
ax.set_yticks(yTicksMarks)
yticknames=ax.set_yticklabels(yTicksMarks)
plt.yticks(fontsize=size_ticks)
ax.legend(prop={'size': size_legend})
#::::::::::::::::::::::::::::::::::::::::::::::::                                                                 

#Letra para la caption
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
ax.text(-0.5, max_y_3 + 0.01, 'd', fontsize=size_letter,fontweight='bold')
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#=============================================================================================


# PLOT (4,0)
#=============================================================================================
ax = plt.subplot(gs[0,4])
sns.despine(offset=0, ax=ax)#Poner el fondo de la grafica en blanco

#Plot
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Baseline_Bars = ax.bar(indBaseline, Baseline_V25, width, color=[cmap1(0.9),cmap1(0.9),cmap1(0.9)],
                            yerr=0, error_kw=dict(elinewidth=2,ecolor='black'), label="Baseline")

MMSBM_Bars = ax.bar(indMMSBM, MMSBM_V25, width, color=[cmap1(0.1),cmap1(0.1),cmap1(0.1)],
                            yerr=0, error_kw=dict(elinewidth=2,ecolor='black'), label="MMSBM")
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#Ejes: limites y nombres
#::::::::::::::::::::::::::::::::::::::::::::::
max_y_4=0.84
ax.set_xlim(-width,ind[-1]+width+0.5*width)
ax.set_ylim(0.65,max_y_4)
#ax.set_ylabel('Performance',fontsize=size_eje)
#::::::::::::::::::::::::::::::::::::::::::::::

#Ticks
#::::::::::::::::::::::::::::::::::::::::::::::::
#xTickMarks = ['Accuracy  ','  Precision', '  Recall']
#ax.set_xticks([pos + 0.5*width for pos in ind] )#Fija la posicion de los xticks con una list comprehension
ax.set_xticks( [pos for pos in ind] )#Fija la posicion de los xticks con una list comprehension
xtickNames = ax.set_xticklabels(xTickMarks)
plt.setp(xtickNames, rotation=0, fontsize=size_ticks)

#yTicksMarks=[ 0.67, 0.69, 0.71, 0.73, 0.75, 0.77,0.79, 0.81, 0.83]
yTicksMarks=[ 0.67, 0.71, 0.75, 0.79, 0.83]
ax.set_yticks(yTicksMarks)
yticknames=ax.set_yticklabels(yTicksMarks)
plt.yticks(fontsize=size_ticks)
ax.legend(prop={'size': size_legend})
#::::::::::::::::::::::::::::::::::::::::::::::::                                                                 

#Letra para la caption
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
ax.text(-0.5, max_y_4+0.01, 'e', fontsize=size_letter,fontweight='bold')
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#=============================================================================================

#plt.savefig('/home/scobo/Dropbox/Tesis/Figures/Comparison_Microbes.pdf', dpi=300)
plt.show()

