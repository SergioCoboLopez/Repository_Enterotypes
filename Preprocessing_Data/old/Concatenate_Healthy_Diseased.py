import numpy as np

#0.1. Read Data
#===============================================
def read_file(inFileName):
    lines = open(inFileName).readlines()
    d = [line.strip().split() for line in lines]
    return d
#===============================================




Input_Path='/export/home/shared/Projects/Microbiome/Input_Data/Datasets_Diseases/'
Datasets=['He_CD','Nielsen_IBD','Vogtmann_CRC','Zeller_CRC','Zhang_RA']

for Dataset in Datasets:
    Name=Dataset.split("_")[0]
    
    Healthy=read_file(Input_Path  + 'Adyacencia' + Name + '_Healthy.txt')
    Diseased=read_file(Input_Path + 'Adyacencia' + Dataset+ '.txt')
    print np.shape(Healthy)
    print np.shape(Diseased)

    
    Whole_Dataset=np.concatenate((Healthy,Diseased), axis=1)

    with open('Adjacency_All_' + Name + '.txt', 'w') as file:
        for line in Whole_Dataset:
            line=' '.join(map(str, line))
            file.write(line)
            file.write('\n')





