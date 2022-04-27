#19/04/2022. Measure Intra/Inter distances for V-10 with BMI

import numpy as np
import pandas as pd
import scipy.cluster
from scipy import stats
from scipy.spatial import distance
import math
import seaborn as sns
import matplotlib.pyplot as plt

#0.FUNCTIONS

#0.1. Read Data                                                                  
def read_file(inFileName):
    lines = open(inFileName).readlines()
    d = [line.strip().split() for line in lines]
    return d

#Intra-class distance
def Measure_Intra_Distances(Dict1): 
    Dict2={}
    #Run over phylum
    for phylum in Dict1:                                                       
        Microbes_per_phylum=len(Dict1[phylum])             
        #Run over pairs of microbes                                
        for microbe_a in range(Microbes_per_phylum):  
            for microbe_b in range(microbe_a+1,Microbes_per_phylum):
               a=np.array(Dict1[phylum][microbe_a][1])               
               b=np.array(Dict1[phylum][microbe_b][1])
               dist = np.linalg.norm(a - b)

               Ids_distance_tuple=(\
            Dict1[phylum][microbe_a][0],Dict1[phylum][microbe_b][0], dist)     

               try:                                                  
                   Dict2[phylum].append(Ids_distance_tuple)   

               except KeyError:                      
                   Dict2[phylum]=[Ids_distance_tuple] 

    return Dict2

#interclass distances
def Measure_Inter_Distances(Dict1):
    Dict_Inter={}
    Orders=Dict1.keys()
    for order_0 in Dict1:
        Orders.remove(order_0)
        for order_1 in Orders:
            for microbe_0 in Dict1[order_0]:
                for microbe_1 in Dict1[order_1]:
                   Dist_0=np.array(microbe_0[1])
                   Dist_1=np.array(microbe_1[1])
                   Dist = np.linalg.norm(Dist_0 - Dist_1)
                   Ids_distance_tuple=(microbe_0[0], microbe_1[0], Dist)
                   try:
                       Dict_Inter[order_0+'-'+order_1].append(Ids_distance_tuple)
                   except KeyError:
                       Dict_Inter[order_0+'-'+order_1]=[Ids_distance_tuple]
    return Dict_Inter 

def Extract_Distances(Dict_Distances,Hosts=0):
    Distances_fun=[]
    Means=[];SEMs=[];stds=[]                      
    xticks=[]; Sizes_Plot=[]
    for phylum in Dict_Distances:
        if Hosts==0:
            #Extract names of phyla and sizes of groups for plotting purposes
            #----------------------------------------------      
            phylum_ticks=phylum.split('__')                        
            xticks.append(phylum_ticks[1])                                        
            Sizes_Plot.append( len(Dict_Distances[phylum]))
            #-----------------------------------------------
        else:
            xticks.append(phylum)
        #---------------------------------------------- 
        Distances_per_Phylum_fun=[]                           
        for distances in Dict_Distances[phylum]:                                 
            Distances_per_Phylum_fun.append(distances[2])                        
        Distances_fun.append(Distances_per_Phylum_fun)                          
        Mean=sum(Distances_per_Phylum_fun)/float(len(Distances_per_Phylum_fun))
        std=np.std(Distances_per_Phylum_fun)    
        SEM=stats.sem(Distances_per_Phylum_fun,ddof=0)
        #----------------------------------------------
        
        #-------------------                                                      
        SEMs.append(SEM)                                                          
        stds.append(std)                                                          
        Means.append(Mean)                                                        
        #-------------------     
    return Distances_fun, Means, SEMs, stds, xticks, Sizes_Plot



#Read Data
K=10;L=20
Hosts=92;Microbes=137

path_in='../../Input_Data/'
path_Out='../../Output_Data/'+ 'Leave_One_Out_V-10_Tot/'              

#Read model data
Likelihood_Seed={}                                                             
for Seed in range(1,10,2):                                                     
    Likelihood_File='Dataset_V-10_Tot_Seed_'+str(Seed)+'_LogLikelihood.txt'
    LikelihoodFinal=read_file(path_Out + Likelihood_File)[-1]                  
    Likelihood_Seed[float(LikelihoodFinal[0])]=Seed                            
    Seed=Likelihood_Seed[max(Likelihood_Seed)]                                 
File_Data='Dataset_V-10_Tot_Seed_' + str(Seed) + '_Parameters.txt'
Data=read_file(path_Out+File_Data)

#Membership vectors
theta=[];eta=[]                                                                
theta=Data[:Hosts]                                                          
eta=Data[Hosts+1:Hosts+1+Microbes]                                       

for microbe in range(len(eta)):                                                
    eta[microbe] = [ float(membership) for membership in eta[microbe] ]        

for host in range(len(theta)):                                                 
    theta[host] = [ float(membership) for membership in theta[host] ]

#Read Metadata and discretize
Metadata_Hosts='Metadata_Hosts_V-10_Tot.csv'
Metadata_Hosts_Raw=pd.read_csv(path_in + Metadata_Hosts)
BMI_Data_Raw=Metadata_Hosts_Raw['BMI']
print(BMI_Data_Raw)

BMI_Data=pd.cut(BMI_Data_Raw,
       bins=[18.5, 20.66, 22.83, 25], 
       labels=["Low", "Medium", "High"])

print(BMI_Data)


#Group by BMI
BMI_Host_Dict={}
for host1 in range(len(BMI_Data)):
    Id_theta_tuple=( host1, theta[host1] )                                   
    try:                                                                       
        BMI_Host_Dict[ BMI_Data[host1] ].append( Id_theta_tuple )            
    except KeyError:                                                           
        BMI_Host_Dict[ BMI_Data[host1] ]=[ Id_theta_tuple ]
print(BMI_Host_Dict)


#2.5.2.3. COMPUTE EUCLIDEAN INTRA-DISTANCES FOR SAME SEX                       
#-------------------------------------------------------------------------     
Intra_Host_Dict={}                                    
Intra_Host_Dict=Measure_Intra_Distances(BMI_Host_Dict) 
Distances,Mean_Distances,SEM_Distances,std_Distances,phylums_xticks,sizes=\
Extract_Distances(Intra_Host_Dict,1)                                           

Intra_Flat_List = [item for sublist in Distances for item in sublist]          
Intra_Flat_List = np.asarray(Intra_Flat_List)
#-------------------------------------------------------------------------     


#3.3.2.3. COMPUTE EUCLIDEAN INTER-DISTANCES FOR PAIRS OF DIFFERENT SEX         
#-------------------------------------------------------------------------     
Inter_Host_Dict={}                                                             
Inter_Host_Dict=Measure_Inter_Distances(BMI_Host_Dict)                        

Distances_Inter,Mean_Distances_Inter,SEM_Distances_Inter,std_Distances_Inter,phylums_xticks_Inter,sizes_Inter=Extract_Distances(Inter_Host_Dict,1) 
Inter_Flat_List = [item for sublist in Distances_Inter for item in sublist]    
Inter_Flat_List = np.asarray(Inter_Flat_List)
#-------------------------------------------------------------------------     

Statistics_Dist_Hosts={'Same BMI':Intra_Flat_List,'Different BMI':Inter_Flat_List}

#Convert dictionary to dataframe for plots 
#--------------------------------------------------------------------------
List_df=[]
for key1 in Statistics_Dist_Hosts:                                            
    print(key1)
    for element in Statistics_Dist_Hosts[key1]:
    	List_df.append([key1,element])

Dataframe_Plot=pd.DataFrame(List_df,columns=\
["Type of Distance","Euclidian Distance"])
print(Dataframe_Plot)

    
#--------------------------------------------------------------------------

#Plot 
#--------------------------------------------------------------------------
sns.set(font_scale=1.2)
sns.set_style("ticks")
sns.barplot(x="Type of Distance",y="Euclidian Distance",data=Dataframe_Plot)
plt.legend(loc='upper left',frameon=False)
sns.despine(top=True, right=True, left=False, bottom=False) 
plt.show()
#--------------------------------------------------------------------------
