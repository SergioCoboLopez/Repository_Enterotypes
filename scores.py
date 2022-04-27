#1/03/2018. Este script toma como input las scores generadas por el
#MixMembership, las compara con el test set y da el porcentaje de
#aciertos

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

#1.1. Name and read files
#====================================================================

#1.1.1. EXTERNAL PARAMETERS
#------------------------------------------------------------
Fold="0"                                                 #---
Seed=sys.argv[1]; K=sys.argv[2]; L=sys.argv[3]           #---
Dataset=sys.argv[4]#Categories, Quartiles, #Toy, #Toy_Old#---
#------------------------------------------------------------

Index=0

if Dataset=='Categories':
    Dataset='5Fold_LargeData/Categories/'
    Output ='LargeDataset/Categories/'
elif Dataset=='Quartiles':
    Dataset='5Fold_LargeData/Quartiles/'
    Output ='LargeDataset/Quartiles/'
elif Dataset=='Toy':
    Dataset='5Fold_ToyData/'
    Output= 'ToyDataset/'
elif Dataset=='Toy_Old':
    Dataset='5Fold_ToyData/OldData/'
    Output ='ToyDataset/OldData'
else:
    print "Aqui hay algo mal, muchacho"

path_scores='../Output_Data/'+Output

Name=Fold+'_Seed_' + Seed + '_K_' + K +'L_' + L + '_scores.txt'

path_check='../Input_Data/' + Dataset

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
   
#print scores_dict


#1.4. Compare scores with test check file
#=========================================================
Hits=0
#log=open(Seed+'log.txt','w')
for line in Check:
   MicrobeCheck=float(line[0]); PersonCheck=float(line[1])
   TuplaCheck=(MicrobeCheck,PersonCheck)
#   print >>log, line[0:3], scores_dict[TuplaCheck]
   if float(line[2])==scores_dict[TuplaCheck]:
      Hits+=1
#=========================================================

print float(Hits)/len(scores)#--->Accuracy Rate

print Hits, len(scores)
   
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
