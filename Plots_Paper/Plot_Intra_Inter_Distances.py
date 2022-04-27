#29/07/2019. En este script cargamos los datos que generamos en el script anterior "Save_Intra-Inter_Distances.py"
#y generamos los plots que nos interesan.
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


start_time = time.time()

nan_saver=1e-10

#1. MAIN
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#1.1. Load Data
#==================================================================
Distance_Data = pickle.load( open( "Pairs_of_Distances.p", "rb" ) )
#==================================================================

#1.2. Compute averages
#===================================================================
Datasets=['S-8_Tot','V-10_Tot','V-22_Tot','V-23_24_Tot','V-25_Tot' ]
Taxonomic_Orders=['Phylum','Class','Order','Family','Genus']

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Metrics=['log_of_means', 'mean_of_logs', 'ratio_of_means']
Metric='ratio_of_means'
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if Metric=='mean_of_logs':

    #CALCULATIONS
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    Distance_logRatios={}
    for Dataset in Distance_Data:
        for taxonomy in Distance_Data[Dataset]:

            A=Distance_Data[Dataset][taxonomy]['Intra']
            B=Distance_Data[Dataset][taxonomy]['Inter']
            All_Pairs=list(product( A,B ))

            LogRatios=[]            
            for pair in All_Pairs:
                Ratio=pair[0]/(pair[1]+nan_saver)
                LogRatio=math.log( Ratio, math.e)
                LogRatios.append( LogRatio )

            Mean_LogRatios=np.mean(LogRatios)
            SEM_LogRatios=scipy.stats.sem(LogRatios,ddof=0)

            try:
                Distance_logRatios[Dataset][taxonomy]={"Mean": Mean_LogRatios, "SEM": SEM_LogRatios }
            except KeyError:
                Distance_logRatios[Dataset]={taxonomy:{"Mean": Mean_LogRatios, "SEM": SEM_LogRatios }}
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #PLOTS
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    xticks_pos=[tick for tick in Taxonomic_Orders  ]  

    #Plot
    #---------------------------------------------------------------
    for Data_plot in Datasets:
        Means=[];SEMs=[]
        for clade in Taxonomic_Orders:
            print Data_plot, clade

            Means.append( Distance_logRatios[Data_plot][clade]["Mean"] )
            SEMs.append( Distance_logRatios[Data_plot][clade]["SEM"] )

        print Means
        plt.errorbar(xticks_pos, Means, yerr=0, marker='o') 
    #---------------------------------------------------------------
    Datasets_Names=['S-8','V-10','V-22','V-23_24','V-25' ]

    #Cosmetics
    #---------------------------------------------------------------
    plt.legend(Datasets)
    plt.ylabel(r'$\langle$log(Intra/Inter)$\rangle$')      
    plt.xlabel('Clades')           
    #---------------------------------------------------------------

    #Save and show                                                   
    #---------------------------------------------------------------
    path_out='/export/home/shared/Projects/Microbiome/Plots_Paper/'
    plt.savefig(path_out + 'Mean_of_Log_Intra-Inter_Averages.pdf',dpi=300)
    plt.show()                                                     
    #---------------------------------------------------------------
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if Metric=='log_of_means':
    

    #CALCULATIONS
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    Distance_Ratios={}

    for Dataset in Distance_Data:
        for taxonomy in Distance_Data[Dataset]:

            A=Distance_Data[Dataset][taxonomy]['Intra']
            B=Distance_Data[Dataset][taxonomy]['Inter']
            All_Pairs=list(product( A,B ))

            Ratios=[]

            for pair in All_Pairs:
                Ratio=pair[0]/(pair[1]+nan_saver)
                Ratios.append( Ratio )


            Mean_Ratios=np.mean(Ratios)
            SEM_Ratios=scipy.stats.sem(Ratios,ddof=0)

            try:
                Distance_Ratios[Dataset][taxonomy]={"Mean": Mean_Ratios, "SEM": SEM_Ratios }
            except KeyError:
                Distance_Ratios[Dataset]={taxonomy:{"Mean": Mean_Ratios, "SEM": SEM_Ratios }}
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


    #PLOTS
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    xticks_pos=[tick for tick in Taxonomic_Orders  ]  

    #Plot
    #---------------------------------------------------------------
    for Data_plot in Datasets:
        Means=[];SEMs=[]
        for clade in Taxonomic_Orders:
            print Data_plot, clade

            Means.append( math.log(Distance_Ratios[Data_plot][clade]["Mean"], math.e) )
            SEMs.append(  math.log(Distance_Ratios[Data_plot][clade]["SEM"] , math.e) )

        print Means
        plt.errorbar(xticks_pos, Means, yerr=0, marker='o') 
    #---------------------------------------------------------------

    #Cosmetics
    #---------------------------------------------------------------
    plt.legend(Datasets)
    plt.ylabel(r'log$\langle$Intra/Inter$\rangle$')      
    plt.xlabel('Clades')                                            
    #---------------------------------------------------------------

    #Save and show
    #----------------------------------------------------------------
    path_out='../../Plots_Paper/'
    plt.savefig(path_out + 'Log_of_Mean_Intra-Inter_Averages.pdf',dpi=300)
    plt.show()                                                     
    #---------------------------------------------------------------
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if Metric=='ratio_of_means':

    #CALCULATIONS
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    Distances={}

    for Dataset in Distance_Data:
        for taxonomy in Distance_Data[Dataset]:

            Intra=Distance_Data[Dataset][taxonomy]['Intra']
            print "Intra", len(Intra)
            Inter=Distance_Data[Dataset][taxonomy]['Inter']
            print "Inter", len(Inter)
            #ERROR
            #-------------------------------------------------------------------------------
            #We do a propagation of errors to compute the actual error:
            #sqrt[ SEM(Intra)^2/Intra^2 + SEM(Inter)^2/Inter^2 ]
            Inter_Mean=np.mean(Inter); SEM_Inter=scipy.stats.sem(Inter,ddof=0)
            Intra_Mean=np.mean(Intra); SEM_Intra=scipy.stats.sem(Intra,ddof=0)

            print Dataset, taxonomy
            print SEM_Intra
            print Intra_Mean
            print SEM_Inter
            print Inter_Mean
            
            #Error
            #........................................................
            Error= math.sqrt( (SEM_Intra/Intra_Mean)**2 + (SEM_Inter/Inter_Mean)**2 )
            #.......................................................
            print "Error", Error

            #-------------------------------------------------------

            try:
                Distances[Dataset][taxonomy]={"Intra_Mean":Intra_Mean,"Inter_Mean":Inter_Mean,"Error":Error }
            except KeyError:
                Distances[Dataset]={taxonomy:{"Intra_Mean":Intra_Mean,"Inter_Mean":Inter_Mean,"Error":Error } }
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


    #PLOTS
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    xticks_pos=[tick for tick in Taxonomic_Orders  ]  

    #Plot
    #---------------------------------------------------------------
    for Data_plot in Datasets:
        Means=[];Error=[]
        for clade in Taxonomic_Orders:
            print Data_plot, clade

            Means.append(\
            math.log(Distances[Data_plot][clade]["Intra_Mean"]/Distances[Data_plot][clade]["Inter_Mean"],math.e) )

            Error.append( (Distances[Data_plot][clade]["Error"]) )

        print Means
        plt.errorbar(xticks_pos, Means, yerr=Error, marker='o') 
    #---------------------------------------------------------------
    Datasets_Names=['S-8','V-10','V-22','V-23_24','V-25' ]
    #Cosmetics
    #---------------------------------------------------------------
    plt.legend(Datasets_Names)
    plt.ylabel(r'log( $\langle$ $d_{\it{Intra}}$ $\rangle$/$\langle d_{\it{Inter}} \rangle$ )')      
    plt.xlabel('Clades')
    #---------------------------------------------------------------

    #Save and show                                    
    #---------------------------------------------------------------
    path_out='../../Plots_Paper/'
    plt.savefig(path_out + 'Mean_Intra-Inter_Averages.pdf',dpi=300)
    plt.show()                                                     
    #---------------------------------------------------------------
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#====================================================================


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


elapsed_time = time.time() - start_time
print elapsed_time

