#23/07/2018. Esto es una adaptacion del script scores para jugar
#con los scores del baseline

import sys

#0. FUNCTIONS
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#0.1. Read File
#======================================================
def read_file(inFileName):
   lines = open(inFileName).readlines()
   d = [line.strip().split() for line in lines]
   return d
#======================================================


#0.1. Clean File
#======================================================
def clean_data(data):
   scores=[]
   for line in data:
      if len(line)==0:
         continue
      line=map(float,line)
      scores.append(line) 
      
   return scores
#======================================================

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#1. MAIN
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Fold="4"

#1.1. Name and read files
#====================================================================
Outputs =['LargeDataset/Categories/','LargeDataset/Quartiles/']
Datasets=['5Fold_LargeData/Categories/','5Fold_LargeData/Quartiles/' ]
Index=1
Output =Outputs[Index]
Dataset=Datasets[Index]

path_scores='/export/home/shared/Projects/Microbiome/Output_Data/'+Output

Name='Baseline_' + Fold + '_scores.txt'

path_check='/export/home/shared/Projects/Microbiome/Input_Data/' + Dataset

raw_scores=read_file(path_scores+Name)
Check=read_file(path_check + Fold + 'TestCheck.txt')
#====================================================================

#1.2. Clean scores file
#============================
scores=clean_data(raw_scores)
#============================


#1.3. Collect scores in dictionary
#==========================================
scores_dict={}
for line in scores:
   Microbe=line[0];Person=line[1]
   Bet=line[2:];Winner=Bet.index(max(Bet))
   Pareja_Tupla=(Microbe, Person)
   scores_dict[ Pareja_Tupla ]=Winner
#==========================================

#1.4. Compare scores with test check file
#=========================================================
Hits=0
#log=open(Seed+'log.txt','w')
for line in Check:
   MicrobeCheck=float(line[0]); PersonCheck=float(line[1])
   TuplaCheck=(MicrobeCheck,PersonCheck)
 #  print >>log, line[0:3], scores_dict[TuplaCheck]
   if float(line[2])==scores_dict[TuplaCheck]:
      Hits+=1
#=========================================================

print float(Hits)/len(scores)#--->Accuracy Rate

print Hits, len(scores)
   
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
