#20/06/2018. Este codigo tiene como objetivo adaptar el fichero de metadata de los microbios a un formato que
#podamos usar para nuestros scripts.

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

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


#1. READ DATA
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#1.1. Leer datos y hacer diccionarios reciprocos: Microbios:{genus} y Genus:{Microbios}
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
filename= '/export/home/shared/Projects/Microbiome/Input_Data/taxonomies_for_toy_dataset.txt'
Data=read_file(filename)

Microbes_Genus={}
Genus_Microbes={}
contador_microbes=0
for line in Data:
    #Get genus and remove quotes
    #=========================================
    genus_with_quote=line[0].split('_')[0] #==
    genus=genus_with_quote[1:]             #==
    #=========================================

    Microbes_Genus[contador_microbes]=genus
    try:
        Genus_Microbes[genus].append(contador_microbes)
    except KeyError:
        Genus_Microbes[genus]=[contador_microbes]

    contador_microbes+=1
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


#1.2. Diccionario con keys=numero de bacterias en un genus. values=numero de genuses con ese numero de bacterias
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Genus_Frequency={}
for line in Genus_Microbes:
    try:
        Genus_Frequency[ len( Genus_Microbes[line] ) ]+=1
    except KeyError:
        try:
            Genus_Frequency[ len( Genus_Microbes[line] )]=1
        except KeyError:
            Genus_Frequency[ len( Genus_Microbes[line] )]={}
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

Genus_Abundance=[]
for genus in Genus_Microbes.values():   
    Genus_Abundance.append( len(genus) )



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
color1=(float(141)/float(255),float(160)/float(255),float(153)/float(203),1)


#PLOTS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#PLOT_1: Distribution of genuses by number of bacteria in them
#====================================================================================

#Plot
#----------------------------------------------
ax0=plt.figure(figsize=(12,12) )
x0 = Genus_Frequency.keys()
y1=Genus_Frequency.values()
sns.barplot(x0, y1, color=color1)
#----------------------------------------------

#Axes
#----------------------------------------------
plt.title("Distribution of genuses with given number of bacteria")
plt.xlabel("Number of bacteria in a genus")
plt.ylabel("Number of genuses")
#----------------------------------------------
plt.show()
#====================================================================================


#Plot_0: Histogram of all genuses with corresponding number of bacteria 
#====================================================================================

#Plot
#-----------------------------------------------------
ax0=plt.figure(figsize=(12,12) )
x0 = [index0 for index0 in range(len(Genus_Microbes))]
y1=sorted(Genus_Abundance, reverse=True)
sns.barplot(x0, y1, color=color1)
#-----------------------------------------------------

#Axes
#----------------------------------------------
plt.xlabel("genus")
plt.ylabel("number of bacteria per genus")
plt.xticks(range(len(Genus_Microbes)), [], rotation='90')
#----------------------------------------------
plt.show()
#====================================================================================

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
