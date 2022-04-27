import pandas as pd
import numpy as np
import os

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
path='../../Input_Data/Leave_One_Out_'
Apellido='_Tot'
Datasets=['S-8'+Apellido,'V-10'+Apellido,'V-22'+Apellido,'V-23_24'+Apellido,\
          'V-25'+Apellido]

#1.1. Generar NOMBRES de filas para Excel
#==============================================================
measure_test=["Likelihood", "Seed"]
#--------------------
Links_Considered=5#--
#--------------------

for Link_Predicted in range(0,Links_Considered):
    for Link_Actual in range(0,Links_Considered):
        measure_test.append('P'+str(Link_Predicted) + '_A'+str(Link_Actual) )
print measure_test
        

tuples=[]
for Dataset in Datasets:
    folder=path+Dataset+'/'
    print folder
    Patients=read_file(folder+'Patient_Record_'+Dataset+'.txt')
    #----------------------------------------------------
    for line in Patients:
        Patient=line[0]
    #----------------------------------------------------
        for measure in measure_test:
            tuples.append((str(Dataset), str(Patient), measure))
#==============================================================

#1.2. Generar Multiindex (niveles de indexado)
#==================================================
index = pd.MultiIndex.from_tuples(tuples, names=\
                ['Dataset', 'Patient','Magnitude'])
#==================================================

#1.3. Crear Dataframe y nombrar a las columnas
#==================================================
df = pd.DataFrame(index=index, \
columns=['Fold0', 'Fold1','Fold2','Fold3','Fold4'] )
#==================================================

#1.4. Guardar fichero
#=========================================
df.to_csv('../../Results_Sheet_' + Apellido + '.csv')
#=========================================

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
