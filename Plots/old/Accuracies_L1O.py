#17/10/2018. Generamos una distribucion de accuracies de los pacientes
#seleccionados para nuestro leave-one-out.

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns



#Read Dataset
df=pd.read_csv('../../Reports/Results_Sheet.csv',index_col=[0,1,2])

#Reduce it to a given Dataset
SubDataset=df.loc['S-8']

#Reduce Dataset to only Accuracies
#==================================
#Sort Patients in Dataset
SubDataset.sort_index(inplace=True)
#Take Only Accuracies by slicing
Accuracies_DF=SubDataset.loc[( slice(None), 'Accuracy') ,]

Means_DF= Accuracies_DF.mean(axis=1)
SEMs_DF = Accuracies_DF.sem(axis=1)

Means=Means_DF.to_dict()
SEMs =SEMs_DF.to_dict()

Results={}
for Patient in zip(Means, SEMs):
    Results[ Patient[0][0]  ]= (Means[Patient[0]],  SEMs[Patient[1]] )

Values=sorted(Results.values())

#for line in Values:
Accuracies=[ line[0] for line in Values]
print Accuracies, min(Accuracies), max(Accuracies)

ax = plt.hist(Accuracies,bins=30)
plt.axvline(x=0.802275238124, linewidth=4, color='k'  )
plt.show()
#print Accuracies


#==================================





    
