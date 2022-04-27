#10/01/2019. Hacemos un plot de las matrices p. Deformamos esa matriz para tener en cuenta la prevalencia de cada membership.

import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib
from pylab import *
from matplotlib import colors
import matplotlib.patches as patches

import seaborn as sns
import networkx as nx
from networkx.algorithms import bipartite
import matplotlib.pyplot as plt
import pickle
import time
import copy
import operator
import sys


#0. FUNCTIONS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#0.1. Read Data
#===============================================
def read_file(inFileName):
    lines = open(inFileName).readlines()
    d = [line.strip().split() for line in lines]
    return d
#===============================================


#0.2.Heatmap Fenotipos
#=====================================================
def Twist_Pmatrix(MatrizBase,MembershipGames,MembershipPlayers,Patients,Microbes):
    rows=Patients;cols=Microbes

    #0.2.1. Convertir las membresias en enteros (y corregir fluctuaciones decimales)
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
    #0.2.1.1. Juegos
    #--------------------------------------------------------------------------------------------
    
    MembershipGames_int=[round(membresia) for membresia in MembershipGames]

    #Corregir fluctuaciones de los decimales
    #....................................................................................
    
    #Mayor
    if sum(MembershipGames_int)>cols:
        #Diccionario con las partes decimales de la suma de las memberships
        #...................................................................
        Decimales={}
        for numero in MembershipGames:
            #Diccionario: {Indice del vector: (numero, parte decimal)}
            Decimales[ MembershipGames.index(numero) ]=math.modf(numero)[0]
        #Diccionario ordenado por partes decimales
        sortedDec = sorted(Decimales.items(), key=operator.itemgetter(1))
        #...................................................................

        Decimal_Dict_counter=0
        while sum(MembershipGames_int)>cols:
            MembershipGames_int[sortedDec[Decimal_Dict_counter][0]]=\
            int(MembershipGames[sortedDec[Decimal_Dict_counter][0]])
            Decimal_Dict_counter+=1
            
    #Menor
    if sum(MembershipGames_int)<cols:
        #Diccionario con las partes decimales de la suma de las memberships
        #..................................................................
        Decimales={}
        for numero in MembershipGames:
            #Diccionario: {Indice del vector: (numero, parte decimal)}
            Decimales[ MembershipGames.index(numero) ]=math.modf(numero)[0]
            #Diccionario ordenado por partes decimales
        sortedDec = sorted(Decimales.items(), key=operator.itemgetter(1),reverse=True)
        #..................................................................

        Decimal_Dict_counter=0
        while sum(MembershipGames_int)<cols:
            MembershipGames_int[sortedDec[Decimal_Dict_counter][0]]=\
            ceil(MembershipGames[sortedDec[Decimal_Dict_counter][0]])
            Decimal_Dict_counter+=1
    #....................................................................................

    #--------------------------------------------------------------------------------------------

    #0.2.2.1. Jugadores
    #--------------------------------------------------------------------------------------------
    MembershipPlayers_int=[round(membresia) for membresia in MembershipPlayers]

    #Corregir fluctuaciones de los decimales
    #....................................................................................

    #Mayor
    if sum(MembershipPlayers_int)>rows:
        #Diccionario con las partes decimales de la suma de las memberships
        #...................................................................
        Decimales={}
        for numero in MembershipPlayers:
            #Diccionario: {Indice del vector: (numero, parte decimal)}
            Decimales[ MembershipPlayers.index(numero) ]=math.modf(numero)[0]
        #Diccionario ordenado por partes decimales
        sortedDec = sorted(Decimales.items(), key=operator.itemgetter(1))
        #...................................................................

        Decimal_Dict_counter=0
        while sum(MembershipPlayers_int)>rows:
            MembershipPlayers_int[sortedDec[Decimal_Dict_counter][0]]=\
            int(MembershipPlayers[sortedDec[Decimal_Dict_counter][0]])
            Decimal_Dict_counter+=1

    #Menor
    if sum(MembershipPlayers_int)<rows:
        #Diccionario con las partes decimales de la suma de las memberships
        #..................................................................
        Decimales={}
        for numero in MembershipPlayers:
            #Diccionario: {Indice del vector: (numero, parte decimal)}
            Decimales[ MembershipPlayers.index(numero) ]=math.modf(numero)[0]
            #Diccionario ordenado por partes decimales
        sortedDec = sorted(Decimales.items(), key=operator.itemgetter(1),reverse=True)
        #..................................................................

        Decimal_Dict_counter=0
        while sum(MembershipPlayers_int)<rows:
            MembershipGames_int[sortedDec[Decimal_Dict_counter][0]]=\
            ceil(MembershipGames[sortedDec[Decimal_Dict_counter][0]])
            Decimal_Dict_counter+=1
    #....................................................................................

    #--------------------------------------------------------------------------------------------
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #NO TENGO NADA CLARA LA UTILIDAD DE ESTE TROZO DE CODIGO
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #0.5.2.Ordenamos el vector de memberships de jugadores de mayor a menor (conservando indices del original)
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #DICCIONARIO DE INDICES: {valor elemento del vector:indice} (ES UNA MEMORIA)
    #---------------------------------------------------------------
#    IndiceMembershipPlayers=dict( (valor,indice) for (valor,indice) in zip(MembershipPlayers_int,range(K)) )
    # #---------------------------------------------------------------

    #ORDENAR DE MAYOR A MENOR
    #---------------------------------------------------------------
    MembershipPlayers_int=sorted(MembershipPlayers_int,reverse=True)
    MembershipGames_int=sorted(MembershipGames_int,reverse=True)
    #---------------------------------------------------------------
        
    # #Construimos diccionario de indices del vector: {valor elemento del vector:indice}
    # IndiceMembershipGames=dict( (valor,indice) for (valor,indice) in zip(MembershipGames_int,range(L)) )
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    #0.2.3. Construimos la pre-matriz en forma de un vector
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::
    fila_matriz=list(); #Inicializar vector y contador para correr sobre el fichero

    contador_i=0    
    for fila in MembershipPlayers_int:
        contador_j=0
        for columna in MembershipGames_int:
            columna=int(columna);#Formatear a entero
            fila_matriz.extend( \
            [float(MatrizBase[ contador_i ][ contador_j ]) ]*columna)
            contador_j+=1
        fila_matriz.extend(fila_matriz[-cols:]*int(round(fila-1)))#Extender vector en filas
        contador_i+=1
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #0.2.4. Reshape a matriz de rowsxcols
    #:::::::::::::::::::::::::::::::::::::::::::::::::::
    MatrizFenotipos=np.reshape(fila_matriz, (rows,cols) )
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::

    return MatrizFenotipos
#=====================================================


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#1. PARAMETERS FOR INPUT/OUTPUT DATA
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#1.1. EXTERNAL PARAMETERS
#=====================
Dataset=sys.argv[1]#==
number_of_links=3  #==
#=====================

#1.2. INTERNAL PARAMETERS
#==================================================================
K=10;L=20

Patients_per_Dataset=\
{'S-8_Tot':107, 'V-10_Tot':92, 'V-22_Tot':467, 'V-23_24_Tot':222,\
 'V-25_Tot':883}

Microbes_per_Dataset=\
{'S-8_Tot':128, 'V-10_Tot':137, 'V-22_Tot':134, 'V-23_24_Tot':118,\
 'V-25_Tot':144}

Patients=Patients_per_Dataset[Dataset]
Microbes=Microbes_per_Dataset[Dataset]
#================================================================== 

#1.3. PATHS FOR FILES
#==============================================================
Path_out='../../Output_Data/'+ 'Leave_One_Out_' + Dataset + '/'
Path_in= '../../Input_Data/' + 'Leave_One_Out_' + Dataset + '/'
#==============================================================

#1.4. CHOOSE BEST LIKELIHOOD
#======================================================================
Likelihood_Seed={}
for Seed in range(1,10,2):
    likelihood_file=\
    'Dataset_' + Dataset + '_Seed_' + str(Seed) + '_LogLikelihood.txt'
    likelihoodfinal=read_file(Path_out + likelihood_file)[-1]
    Likelihood_Seed[float(likelihoodfinal[0])]=Seed

print Likelihood_Seed
print max(Likelihood_Seed)
print Likelihood_Seed[max(Likelihood_Seed)]
Seed=Likelihood_Seed[max(Likelihood_Seed)]

file='Dataset_' + Dataset + '_Seed_' + str(Seed) + '_Parameters.txt'
Data=read_file(Path_out+file)
#====================================================================== 


#1.5. READ PMATRICES AND FORMAT INTO MATRICES
#======================================================================
counter=0
pMatrices=[]
for line in Data:
    if 'slice' in line:
        pmatrix_raw=np.zeros((K,L))
        for row in range(K):
            for col in range(L):
                pmatrix_raw[row][col]=float(Data[counter+row+1][col])

        pMatrices.append(pmatrix_raw)

    counter+=1
#======================================================================

#1.6. COMPUTE THETAS, ETAS AND AGGREGATED
#=================================================
theta=[];contador_theta=0
eta=[];contador_eta=0

for line in Data[:Patients]:
    line.insert(0, contador_theta )
    contador_theta+=1

theta=Data[:Patients]

for line in Data[Patients+1:Patients+1+Microbes]:
    line.insert(0, contador_eta )
    contador_eta+=1

eta=Data[Patients+1:Patients+1+Microbes]

theta_aggregate=[0]*K
for lineT in theta:
    for i in range(K):
        theta_aggregate[i]+=float(lineT[i+1])

eta_aggregate=[0]*L
for lineE in eta:
    for j in range(L):
        eta_aggregate[j]+=float(lineE[j+1])
#=================================================

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#2. Compute reshaped p-matrices
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
MatrizFenotipos_0=Twist_Pmatrix(pMatrices[0],eta_aggregate,theta_aggregate,Patients,Microbes)
MatrizFenotipos_1=Twist_Pmatrix(pMatrices[1],eta_aggregate,theta_aggregate,Patients,Microbes)
MatrizFenotipos_2=Twist_Pmatrix(pMatrices[2],eta_aggregate,theta_aggregate,Patients,Microbes)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#3 Reduce to a single matrix. Only abundance vs. non-abundance
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Contrast=np.ones((K, L))
AbundanceMatrix=Contrast - pMatrices[0]

#Dictionary of rows: sum of columns
#::::::::::::::::::::::::::::::::::::::::
sum_cols={}
for row in range(K):
    sum_cols[row]=sum(AbundanceMatrix[row])
#::::::::::::::::::::::::::::::::::::::::

#Order the dictionary by ascending sum of columns
#::::::::::::::::::::::::::::::::::::::::::::::::
sorted_sum_cols = sorted(sum_cols.items(), key=operator.itemgetter(1),reverse=True)
#::::::::::::::::::::::::::::::::::::::::::::::::

#Reorder matrix rows by the sum of the columns
#::::::::::::::::::::::::::::::::::::::::::::::::::::::
print sorted_sum_cols
Dict_row_membership={}
AbundanceMatrix_row_ordered=np.zeros(( K,L ))
for new_row in range(K):
    AbundanceMatrix_row_ordered[new_row] = AbundanceMatrix[ sorted_sum_cols[new_row][0] ]
    Dict_row_membership[new_row]=sorted_sum_cols[new_row][0]
print Dict_row_membership
#::::::::::::::::::::::::::::::::::::::::::::::::::::::

#Dictionary of cols: sum of rows
#::::::::::::::::::::::::::::::::::::::::
sum_rows={}
for col in range(L):
    sum_rows[col]=sum( np.transpose(AbundanceMatrix_row_ordered)[col] )
#::::::::::::::::::::::::::::::::::::::::

#Order the dictionary by ascending sum of rows
#::::::::::::::::::::::::::::::::::::::::::::::::
sorted_sum_rows = sorted(sum_rows.items(), key=operator.itemgetter(1),reverse=True)
#::::::::::::::::::::::::::::::::::::::::::::::::

print sorted_sum_rows

#Reorder matrix rows by the sum of the columns
#::::::::::::::::::::::::::::::::::::::::::::::::::::::
Dict_column_membership={}
AbundanceMatrix_ordered=np.zeros(( K,L ))
for new_col in range(L):
    np.transpose(AbundanceMatrix_ordered)[new_col] = \
    np.transpose(AbundanceMatrix_row_ordered)[ sorted_sum_rows[new_col][0] ]
    
    Dict_column_membership[new_col]=sorted_sum_rows[new_col][0]

print Dict_column_membership
#::::::::::::::::::::::::::::::::::::::::::::::::::::::

#Binarize Abundance Matrix
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
for row_binary in range(K):
    for col_binary in range(L):
        if AbundanceMatrix_ordered[row_binary][col_binary] < 1.0 and \
           AbundanceMatrix_ordered[row_binary][col_binary] >= 0.5:
            AbundanceMatrix_ordered[row_binary][col_binary]=1

        elif AbundanceMatrix_ordered[row_binary][col_binary] > 0.0 and \
           AbundanceMatrix_ordered[row_binary][col_binary] < 0.5:
            AbundanceMatrix_ordered[row_binary][col_binary]=0
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#Shapes and colors for groups of microbes/patients
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
colors_shapes_microbes={}
shapes=["^", "v"]
counter=0
        
for shape in range(len(shapes)):
    for group in range(L/2):
        colors_shapes_microbes[counter]=( group, shapes[shape] )
        counter+=1
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


#Compute Nestedness from ordered matrix
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
d_ij_microbes=0;d_ij_patients=0;
Min_Interactions_microbes=0;Min_Interactions_patients=0;

for row_i in range(L-1):
    for row_j in range(row_i+1,L):

        #Numerator
        #---------------------------------------------------------------------------------------
        #Microbes
        #.......................................................................
        if row_j<K:
            product_microbes=np.multiply(AbundanceMatrix_ordered[row_i], AbundanceMatrix_ordered[row_j])
            shared_interactions_microbes=np.sum( product_microbes)
            d_ij_microbes+=shared_interactions_microbes
        #.......................................................................

        #Patients
        #.......................................................................
        product_patients=np.multiply(np.transpose(AbundanceMatrix_ordered)[row_i],\
                                     np.transpose(AbundanceMatrix_ordered)[row_j])
        shared_interactions_patients=np.sum( product_patients)
        d_ij_patients+=shared_interactions_patients
        #.......................................................................

        #---------------------------------------------------------------------------------------
        
        #Denominator
        #---------------------------------------------------------------------------------------
        if row_j<K:
            Interactions_microbes_i=np.sum(AbundanceMatrix_ordered[row_i])
            Interactions_microbes_j=np.sum(AbundanceMatrix_ordered[row_j])
            Min_Interactions_microbes+=min(Interactions_microbes_i, Interactions_microbes_j)

        Interactions_patients_i=np.sum(np.transpose(AbundanceMatrix_ordered)[row_i])
        Interactions_patients_j=np.sum(np.transpose(AbundanceMatrix_ordered)[row_j])
        Min_Interactions_patients+=min(Interactions_patients_i, Interactions_patients_j)
        #---------------------------------------------------------------------------------------


Nestedness=float(d_ij_microbes+d_ij_patients)/float(Min_Interactions_microbes+Min_Interactions_patients)        
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#3. PLOTS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Information=["Single_Pmatrix_Complenet_Presentation", "Single_Pmatrix_w_membershipts" , "Pmatrices", "Network"]
Plot=Information[1]

#3.1. Plot abundance vs non-abundance matrix with nestedness and membership labels
#===============================================================================================================
if Plot=="Single_Pmatrix_Complenet_Presentation":
    ax=sns.heatmap(AbundanceMatrix_ordered,cmap="RdBu_r",linewidths=1,vmin=0,vmax=1,\
                   xticklabels=False,yticklabels=False)

    ax.text(0.75, 0.75, 'nestedness= %f' %Nestedness,color='white',fontsize=12,
                    bbox={'facecolor':'red', 'alpha':0.5, 'pad':2.6})

    #---------------------------------------------
    ax.set_xlabel("Communities of Microbes",fontsize=15)
    ax.set_ylabel("Communities of Patients",fontsize=15)
    #---------------------------------------------
    plt.savefig(\
    '/home/scobo/Dropbox/beamer_template/COMPLENET/Plots/Mutualistic_P_Matrix_'+Dataset+'.png',dpi=300)
    plt.show()
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#===============================================================================================================


#3.2. Plot abundance vs non-abundance matrix with nestedness and membership labels
#===============================================================================================================
if Plot=="Single_Pmatrix_w_membershipts":

    #3.2.1. Gridspec Parameters
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    fig=plt.figure(figsize=(16,7))
    gs = gridspec.GridSpec(2, 2, height_ratios=[0.15,1],width_ratios=[0.15,1])
    gs.update(left=0.05, right=0.999,top=0.95,bottom=0.1, wspace=0.05,hspace=0.1)
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
    #3.2.2. Plot Triangles
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    triangles_up = plt.subplot(gs[0, 1])
    sns.despine(offset=0, ax=triangles_up, top=True, left=True,bottom=True)#Poner fondo de la grafica en blanco

    #3.2.2.1. Triangle parameters and colormap
    #-----------------------------------------
    Base=0.03; Height=0.45
    Buffer=0.005
    Right_Vertex=[-Buffer, 0]

    colors_microbes = cm.Paired(np.linspace(0, 1, K))
    #-----------------------------------------

    #3.2.2.2. Loop to draw triangles
    #---------------------------------------------------------------------------------------------------
    for triangle in range(L):

        #Get membership corresponding to column
        #....................................................
        Membership_of_column=Dict_column_membership[triangle]
        #....................................................

        #Set up,down triangles according to membership: shape, color dictionary
        #.....................................................................................
        if colors_shapes_microbes[Membership_of_column][1]=="^":
            Origin=[Right_Vertex[0]+2*Buffer,0.0*Height];Right_Vertex=[Origin[0]+Base, Origin[1]]
            Top_Vertex=[Right_Vertex[0]-0.5*Base,Origin[1]+ Height]

        else:
            Origin=[Right_Vertex[0]+2*Buffer,Height];Right_Vertex=[Origin[0]+Base, Origin[1]]
            Top_Vertex=[Right_Vertex[0]-0.5*Base,Origin[1]- Height]
        #.....................................................................................

        #Draw and color triangle
        #.............................................................................
        Vertices=[ Origin, Right_Vertex, Top_Vertex, Origin]
        triangles_up.add_patch(patches.Polygon( Vertices,\
        facecolor=colors_microbes[ colors_shapes_microbes[Membership_of_column][0]] ))
        #.............................................................................
        
    #---------------------------------------------------------------------------------------------------

    #3.2.2.3. Put title and remove ticks and axis
    #-----------------------------------------------------------------
    plt.title(Dataset+"            ",fontsize=18)

    #Quitar xticks,yticks y ejes
    triangles_up.tick_params(
        axis='both',         # changes apply to the x-axis
        which='both',        # both major and minor ticks are affected
        bottom='off',        # ticks along the bottom edge are off
        top='off',           # ticks along the top edge are off
        labelbottom='off',
        left='off',
        labelleft='off'     )
    #-----------------------------------------------------------------
    
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #3.2.3. PLOT SQUARES
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    squares_left=plt.subplot(gs[1, 0])
    sns.despine(offset=0, ax=squares_left, top=True, left=True,bottom=True)#Poner fondo de la grafica en blanco
    
    #3.2.3.1. Square parameters and colormap
    #----------------------------------------------
    Origin_x_left=0.55;
    Width_left=0.18;Height_left=0.06;
    Separation=0.02
    Origin_y_left=1-Separation
    colors_patients = cm.Set3(np.linspace(0, 1, K))
    #----------------------------------------------

    #3.1.3.2. Loop to draw squares
    #-----------------------------------------------------------------
    for square in range(K):
        print square, Dict_row_membership[square]
        
        squares_left.add_patch(patches.Rectangle(\
        (Origin_x_left,Origin_y_left),Width_left,-Height_left,\
        facecolor=colors_patients[ Dict_row_membership[square] ]))

        Origin_y_left=Origin_y_left - Height_left - 2*Separation
    #-----------------------------------------------------------------
    
    #3.2.3.3. Put title and remove ticks and axis
    #-----------------------------------------------------------------
    #Quitar xticks,yticks y ejes
    squares_left.tick_params(
        axis='both',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom='off',      # ticks along the bottom edge are off
        top='off',         # ticks along the top edge are off
        labelbottom='off',
        left='off',
        labelleft='off'     )
    #-----------------------------------------------------------------
    
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #3.2.3. Plot Mutualistic Network
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ax = plt.subplot(gs[1, 1])
    ax=sns.heatmap(AbundanceMatrix_ordered,cmap="RdBu",linewidths=1,vmin=0,vmax=1)

    ax.text(0.75, 0.75, 'n= %f' %Nestedness,color='white',fontsize=15,
                    bbox={'facecolor':'red', 'alpha':0.5, 'pad':2.6})

    #yticks, xticks, ylabels, xlabels
    #---------------------------------------------
    yticks_pos=[row + 0.5 for row in range(K)]
    yticks_labels=Dict_row_membership.values()
    ax.set_yticks(yticks_pos)
    ax.set_yticklabels(yticks_labels,fontsize=15)

    xticks_pos=[col + 0.5 for col in range(L)]
    xticks_labels=Dict_column_membership.values()
    ax.set_xticks(xticks_pos)
    ax.set_xticklabels(xticks_labels,fontsize=15)
    
    ax.set_xlabel("Groups of Microbes",fontsize=18)
    ax.set_ylabel("Groups of Patients",fontsize=18)
    #---------------------------------------------

    plt.savefig(\
    '../../Plots/First_Path/Mutualistic_P_Matrix_w_Memberships'+Dataset+'.pdf',dpi=300)
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    plt.show()
#===============================================================================================================

#3.3. Plot all three p-matrices
#==============================================================================    
elif Plot=="Pmatrices":
    
    fig=plt.figure(figsize=(17,5))
    gs = gridspec.GridSpec(1, 3)
    gs.update(left=0.05,right=0.97,bottom=0.1,top=0.94,wspace=0.25,hspace=0.3)

    ticks_x=range(0,Microbes)
    ticks_y=range(0,Patients)

    #3.3.0. (0,0) Plot
    #:::::::::::::::::::::::::::::::::::::::::::::::::
    pmatrix0=plt.subplot(gs[0,0])
    sns.heatmap(MatrizFenotipos_0,cmap="RdBu",vmin=0,vmax=1)
    plt.title("0 Abundance")
    #:::::::::::::::::::::::::::::::::::::::::::::::::

    #3.3.1. (0,1) Plot
    #:::::::::::::::::::::::::::::::::::::::::::::::::
    pmatrix1=plt.subplot(gs[0,1])
    sns.heatmap(MatrizFenotipos_1,cmap="RdBu",vmin=0,vmax=1)
    plt.title("Low Abundance")
    #:::::::::::::::::::::::::::::::::::::::::::::::::

    #3.3.2. (0,2) Plot
    #:::::::::::::::::::::::::::::::::::::::::::::::::
    pmatrix2=plt.subplot(gs[0,2])
    sns.heatmap(MatrizFenotipos_2,cmap="RdBu",vmin=0,vmax=1)
    plt.title("High Abundance")
    #:::::::::::::::::::::::::::::::::::::::::::::::::

    plt.savefig(\
    '../../Plots/First_Path/Adjusted_P_Matrices_'+Dataset+'.pdf',dpi=300)

    plt.show()
#==============================================================================


#3.3. Plot bipartite network of groups of patients/microbes
#==============================================================================
elif Plot=="Network":

    Bipartite_Network=nx.Graph()

    #3.3.1. Add Nodes
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ScalerPatients=3000
    Dict_node_patients={}
    for node_patients in range(0,K):
        Label_P='Patients_'+str(node_patients)
        Size_P=ScalerPatients*theta_aggregate[node_patients]/sum(theta_aggregate)
        Dict_node_patients[Label_P]={'bipartite':0, 'size':Size_P,'shape':'d', 'color':node_patients,\
                                     'cmap':'Reds'}
        #Add nodes one-by-one
        Bipartite_Network.add_node(Label_P,size=Size_P,bipartite=0, color=node_patients, shape='d',cmap='Reds')


    ScalerMicrobes=2*ScalerPatients
    Dict_node_microbes={}
    for node_microbes in range(0,L):
        Label_M='Microbes_'+str(node_microbes)
        Size_M=ScalerMicrobes*eta_aggregate[node_microbes]/sum(eta_aggregate)
        Dict_node_microbes[Label_M]={'bipartite':1, 'size':Size_M,'shape':'o', 'color':node_microbes,\
                                     'cmap':'Greys'}
        #Add nodes one-by-one
        Bipartite_Network.add_node(Label_M,size=Size_M,bipartite=1,color=node_microbes, shape='o',cmap='Greys')

    nodeShapes = set((aShape[1]["shape"] for aShape in Bipartite_Network.nodes(data = True)))
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    
    #3.3.3. Get node sizes
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::
    size=[];color=[]
    for node in Bipartite_Network.nodes:
        size.append(Bipartite_Network.nodes[node]['size'])
        color.append(Bipartite_Network.nodes[node]['color'])
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::

    #3.3.3. Add edge weights
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::
    List_Links=[]
    for node_patient in Dict_node_patients.keys():
        for node_microbe in Dict_node_microbes.keys():

            Number_Patient=int(node_patient.split('_')[1])
            Number_Microbe=int(node_microbe.split('_')[1])

            Real_Link=1 - pMatrices[0][Number_Patient][Number_Microbe]

            if Real_Link==0:
                continue

            Node_Link_Node=(node_patient, node_microbe, Real_Link)
            List_Links.append(Node_Link_Node)

    Bipartite_Network.add_weighted_edges_from(List_Links)
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::

    #3.3.4. Draw Network
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::

    #3.3.4.1. Assign positions for plot
    #---------------------------------------------------------------------------------
    nodePos = nx.layout.spring_layout(Bipartite_Network,k=2.0,iterations=100,seed=123)
    pos_labels={}
    for position in nodePos:
        pos_labels[position]=[ nodePos[position][0], nodePos[position][1]+0.05 ]
    #---------------------------------------------------------------------------------

    plt.figure(figsize=(20,12))
    plt.title(str(Dataset))

    #3.3.4.2 Draw nodes
    #---------------------------------------------------------------------------------
    for aShape in nodeShapes:
        #...filter and draw the subset of nodes with the same symbol in the
        #positions that are now known through the use of the layout.
        nx.draw_networkx_nodes(Bipartite_Network,nodePos,node_size=size,node_shape = aShape,\
        nodelist = [sNode[0] for sNode in \
                    filter(lambda x: x[1]["shape"]==aShape,Bipartite_Network.nodes(data = True))],\
        node_color=[sNode[1]['color'] for sNode in \
                    filter(lambda x: x[1]["shape"]==aShape,Bipartite_Network.nodes(data = True))],
        cmap=sNode[1]['cmap'])
    #---------------------------------------------------------------------------------

    #3.3.4.3. Draw edges
    #---------------------------------------------------------------------------------
    Full_Link=[(u,v) for (u,v,d) in Bipartite_Network.edges(data=True) if d['weight'] ==1.0]
    elarge=[(u,v) for (u,v,d) in Bipartite_Network.edges(data=True) if d['weight'] >0.5 and d['weight']<1.0]
    esmall=[(u,v) for (u,v,d) in Bipartite_Network.edges(data=True) if d['weight'] <=0.5]

    nx.draw_networkx_edges(Bipartite_Network,nodePos,edgelist=esmall,width=0.5,edge_color='k',style='-.')
    nx.draw_networkx_edges(Bipartite_Network,nodePos,edgelist=elarge,width=1.0, edge_color='k', style='--')
    nx.draw_networkx_edges(Bipartite_Network,nodePos,edgelist=Full_Link,width=1.0, edge_color='k')
    #---------------------------------------------------------------------------------

    #3.3.4.4. Draw labels
    #----------------------------------------------------------------
    labels=nx.draw_networkx_labels(Bipartite_Network,pos=pos_labels)
    #----------------------------------------------------------------

    #3.3.4.5. Remove frame and axes
    #---------------
    plt.axis('off')
    #---------------

    plt.savefig(\
    '../../Plots/First_Path/Metanetworks_'+Dataset+'.pdf',dpi=300)
    plt.show()
#==============================================================================

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

