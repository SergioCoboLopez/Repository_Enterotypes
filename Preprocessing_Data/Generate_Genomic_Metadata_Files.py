#31/3/2020. Cruzamos la metadata de los microbios con la genomica para generar
#un fichero que nos diga especie microbial + informacion genetica.


import numpy as np
import random as rnd
import copy as cp
import sys
import matplotlib.gridspec as gridspec
import matplotlib
from pylab import *
from matplotlib import colors
import matplotlib.patches as patches
import pandas as pd
import scipy.cluster
from scipy import stats
from sklearn import datasets, linear_model
from sklearn.metrics import mean_squared_error, r2_score
from scipy.special import expit


import seaborn as sns
import matplotlib.pyplot as plt
import operator
import time
import collections
import jellyfish



Datasets=['S-8_Tot', 'V-10_Tot', 'V-22_Tot', 'V-23_24_Tot', 'V-25_Tot']
Dataset=Datasets[0]

#1.1. Taxonomic and Genomic Information

#1.1.1. Taxonomic Information
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
File_Taxon='Names_Microbes_' + Dataset + '.csv'
Taxon_Data_Raw=pd.read_csv('../../Input_Data/' + File_Taxon,header=None)
Species_Names=Taxon_Data_Raw[6]
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#1.1.2. Genomic Information
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
File_Genomics='Genomes_Full.csv'
Genomic_Data_Raw=pd.read_csv('../../Input_Data/'+File_Genomics)#,header=None)
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

Name_Similarity={};counter=0
for Name in Species_Names:

    for index,row in Genomic_Data_Raw.iterrows():
        # print row[1]
        # print row[4]
        Name_Genomics=row[2]
        Similarity= jellyfish.levenshtein_distance(unicode(Name_Genomics)\
                                                   ,unicode(Name))
        Name_Similarity[Name_Genomics]=Similarity
        

    sorted_Name_Similarity = \
    sorted(Name_Similarity.items(),key=operator.itemgetter(1))

    if sorted_Name_Similarity[0][1]>0:
        counter+=1
        print counter
        with open("Different_Microbes.txt",'a') as file:
            file.write(Name + '\n')
            file.write( str(sorted_Name_Similarity[0])+ '\n')
            file.write( str(sorted_Name_Similarity[1])+ '\n')
            file.write( str(sorted_Name_Similarity[2])+ '\n')
            file.write( str(sorted_Name_Similarity[3])+ '\n')
            file.write( str(sorted_Name_Similarity[4])+ '\n')
            file.write("______________________________" + '\n')

print counter

    
