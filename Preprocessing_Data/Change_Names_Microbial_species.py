import pandas as pd

path='../../Input_Data/'

#Read file
#===========================================
File=pd.read_csv(path + 'Genomes_Full.csv')
#===========================================

#Change species name to match original format
#==============================================================
#File["#Organism Name"]=File["#Organism Name"].apply(lambda x: \
#                                "s__" + x.replace(" ","_"))

File["#Organism Name"]=File["#Organism Name"].apply(lambda x: \
                            x.replace("[",""))

File["#Organism Name"]=File["#Organism Name"].apply(lambda x: \
                            x.replace("]",""))
#==============================================================

#Save modified file
#===========================================
File.to_csv(path + "Genomes_Full_pangolin.csv")
#===========================================    
