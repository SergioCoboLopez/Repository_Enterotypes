#17/04/2022. Adaptacion y correcciones de"Save_Intra-Inter_Distances_Hosts.py"
#Generamos tres plots alternativos con unas comparaciones entre distancias intra e inter de acuerdo a
#determinados criterios. En principio, nos quedamos con el logaritmo neperariano del cociente entre la media de 
#las distancias intra y las distancias inter. A tal efecto, tenemos que calcular el error de esa funcion
#logaritmica. Lo hacemos con lo que se llama la propagacion de errores.

import pickle
import time
from itertools import product
import math
import numpy as np
import scipy.cluster
from scipy import stats
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


start_time = time.time()

nan_saver=1e-10

#1. MAIN
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#1.1. Load Data
#==================================================================
Distance_Data = pickle.load( open( "Pairs_of_Distances_Hosts.p", "rb" ) )
#==================================================================

#1.2. Compute averages
#===================================================================
Datasets=['V-10_Tot','V-23_24_Tot']
Type_of_Metadata_Hosts=['Sex']

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Metrics=['log_of_means', 'mean_of_logs', 'ratio_of_means']
Metric='ratio_of_means'
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


#CALCULATIONS
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Dist={}
Distances={}
Error_Dict={}
for Dataset in Distance_Data:
    print Distance_Data[Dataset]

    Intra=Distance_Data[Dataset]['Intra']
    Inter=Distance_Data[Dataset]['Inter']

    #ERROR
    #-------------------------------------------------------------------------
    #We do a propagation of errors to compute the actual error:
    #sqrt[ SEM(Intra)^2/Intra^2 + SEM(Inter)^2/Inter^2 ]
    Inter_Mean=np.mean(Inter); SEM_Inter=scipy.stats.sem(Inter,ddof=0)
    Intra_Mean=np.mean(Intra); SEM_Intra=scipy.stats.sem(Intra,ddof=0)

    # print Dataset
    # print SEM_Intra
    # print Intra_Mean
    # print SEM_Inter
    # print Inter_Mean

    #Error
    #........................................................
    Error= math.sqrt( (SEM_Intra/Intra_Mean)**2 + (SEM_Inter/Inter_Mean)**2 )
    #.......................................................
    print "Error", Error
    #-------------------------------------------------------

    Distances[Dataset]={"Intra_Mean":Intra_Mean,"Inter_Mean":Inter_Mean,"Error":Error }

    Error_Dict[Dataset]={"Err_Intra":SEM_Intra, "Err_Inter":SEM_Inter }
    Dist[Dataset]={"Intra_Mean":Intra_Mean,"Inter_Mean":Inter_Mean}

    print(Dist)
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#====================================================================

Test=pd.DataFrame(Dist)
print(Dist)

Error_df=pd.DataFrame(Error_Dict)
# print(Error)


# df = pd.DataFrame([\
# ['V-10','Intra',Dist['V-10_Tot']['Intra_Mean'],Dist['V-10_Tot']['Err_Intra']],\
# ['V-10','Inter',Dist['V-10_Tot']['Inter_Mean'],Dist['V-10_Tot']['Err_Inter']],\
# ['V-23_24','Intra',Dist['V-23_24_Tot']['Intra_Mean'],Dist['V-23_24_Tot']['Err_Intra']],\
# ['V-23_24','Inter',Distances['V-23_24_Tot']['Inter_Mean'],Dist['V-23_24_Tot']['Err_Inter']]]\
# , columns=['Dataset', 'Distance', 'val', 'Error'])

#print df

#yerr=df['Error'].tolist()
yerr=Error_df
Test.plot(kind='bar',yerr=yerr, rot=0)
plt.show()
#yerr=
                      
# sns.barplot(data=df, x='Dataset', y='val', hue='Distance')
# plt.show()
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


elapsed_time = time.time() - start_time
print elapsed_time


