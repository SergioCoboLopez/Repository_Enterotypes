import pandas as pd


path="../../Input_Data/"
Prokaryotes=pd.read_csv(path + 'prokaryotes_clean.csv')
Genomes=pd.read_csv(path + 'genomes.csv')

Full = Prokaryotes.merge(Genomes, on='#Organism Name')
Full.to_csv(path + "Genomes_Full.csv", index=False)
