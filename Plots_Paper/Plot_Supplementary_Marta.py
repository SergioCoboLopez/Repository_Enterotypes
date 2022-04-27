import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import math 
import numpy as np
from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition,
                                                  mark_inset)
from math import *
import scipy.stats
import pandas as pd


df = pd.read_csv ('Comparison_MMSBM_RClustering.csv')
#print(df.head())

Datasets=['S-8','V-10','V-22','V-23_24','V-25']
Datasetst=['Liu','Qin','Schirmer','Huttenhower/Lloyd-Price','Zeevi' ]
size_letter=15                                                      
Letters=['a','b','c','d','e']

dfdr={}
for d in Datasets:
    dfd=df[df.Dataset==d]
    dfdsbm=dfd[dfd.Method=='MMSBM']
    dfdc=dfd[dfd.Method=='Clustering']
    dfdcr=dfdc[dfdc.Metric=='Recall']
    dfdsbmr=dfdsbm[dfdsbm.Metric=='Recall']
    dfdr[d]=pd.merge(dfdcr,dfdsbmr,on='Patient')

dfdp={}
for d in Datasets:
    dfd=df[df.Dataset==d]
    dfdsbm=dfd[dfd.Method=='MMSBM']
    dfdc=dfd[dfd.Method=='Clustering']
    dfdcp=dfdc[dfdc.Metric=='Precision']
    dfdsbmp=dfdsbm[dfdsbm.Metric=='Precision']
    dfdp[d]=pd.merge(dfdcp,dfdsbmp,on='Patient')

dfda={}
for d in Datasets:
    dfd=df[df.Dataset==d]
    dfdsbm=dfd[dfd.Method=='MMSBM']
    dfdc=dfd[dfd.Method=='Clustering']
    dfdca=dfdc[dfdc.Metric=='Accuracy']
    dfdsbma=dfdsbm[dfdsbm.Metric=='Accuracy']
    dfda[d]=pd.merge(dfdca,dfdsbma,on='Patient')
   
 
sns.set(font_scale=1.75,style='ticks')
sns.set(style='ticks',font_scale=1.25)
fig = plt.figure(constrained_layout=True,figsize=(13,7.5))
spec = gridspec.GridSpec(ncols=5, nrows=3, figure=fig)
cmap2=sns.color_palette("magma", 3)

Counter_Plt=0
for d in Datasets:
    f_ax = fig.add_subplot(spec[0, Counter_Plt])

    plt.ylim(0.2, 0.9)
    plt.xlim(0.2, 0.9)

    Lim_x_up=plt.gca().get_xlim()                                    
    Lim_y_up=plt.gca().get_ylim()
    f_ax.plot (Lim_x_up,Lim_y_up,color='grey')

    # x= dfda[d]['Performance_x']
    # y= dfda[d]['Performance_y']

    x= dfdp[d]['Performance_x']
    y= dfdp[d]['Performance_y']
    z1=x/y # do opposite because of colors red
    z2=[1 if z>1 else 0 for z in z1 ]
    z3=[-math.log(z) for z in z1]
    print(d,np.mean(z3),scipy.stats.sem(z3))
    nn=sum(z2)/len(z2) #x= clustering y = MMSBM
    npos=1-nn
    print(npos,nn)
    ax1=sns.scatterplot(x,y, s=6,alpha=1,hue=z2,legend=False)
    plt.title(Datasetst[Counter_Plt])
    plt.ylabel('')#, fontsize=15)
    plt.xlabel('')#, fontsize=15)
    if Counter_Plt==0:                                              
        plt.ylabel('MMSBM')#, fontsize=15)
    if Counter_Plt==4:
        plt.ylabel('ACCURACY',labelpad=10,fontsize=20)
        ax1.yaxis.set_label_position("right")
   
    ##Inset
   
    axins = inset_axes(f_ax, width="20%", height="40%", loc=4,borderpad=1)
    axi=sns.barplot([0,0.5],[npos,nn],hue=[0,1],dodge=False)
    sns.despine(trim=True,ax=axi)
    #plt.xlim(-0.5,1)
    for bar in axi.patches:
        width=0.75
        bar.set_width(width)
    axins.set_xticklabels([])
    axins.legend_.remove()
    #axins.tick_params(labelleft=False, labelbottom=False)
   # sns.despine
   

   
    f_ax = fig.add_subplot(spec[1, Counter_Plt])

    plt.ylim(0.2, 0.9)
    plt.xlim(0.2, 0.9)

    Lim_x_up=plt.gca().get_xlim()                                    
    Lim_y_up=plt.gca().get_ylim()
    f_ax.plot (Lim_x_up,Lim_y_up,color='grey')

    x= dfdp[d]['Performance_x']
    y= dfdp[d]['Performance_y']
    z1=x/y # do opposite because of colors red
    z2=[1 if z>1 else 0 for z in z1 ]
    nn=sum(z2)/len(z2) #x= clustering y = MMSBM
    npos=1-nn
    #print(np,nn)

    ax1=sns.scatterplot(x,y, s=6,alpha=1,hue=z2,legend=False)
    plt.ylabel('')#, fontsize=15)
    plt.xlabel('')#, fontsize=15)
 

    if Counter_Plt==0:                                              
        plt.ylabel('MMSBM')#, fontsize=15)
    if Counter_Plt==4:
        plt.ylabel('PRECISION',labelpad=10,fontsize=20)
        ax1.yaxis.set_label_position("right")
               
    ##Inset
    axins = inset_axes(f_ax, width="20%", height="40%", loc=4,borderpad=1)
    axi=sns.barplot([0,0.5],[npos,nn],hue=[0,1],dodge=False)
    sns.despine(trim=True,ax=axi)
    #plt.xlim(-0.5,1)
    for bar in axi.patches:
        width=0.75
        bar.set_width(width)
    axins.set_xticklabels([])
    axins.legend_.remove()


    f_ax = fig.add_subplot(spec[2, Counter_Plt])

    plt.ylim(0.2, 0.9)
    plt.xlim(0.2, 0.9)

    Lim_x_up=plt.gca().get_xlim()                                    
    Lim_y_up=plt.gca().get_ylim()
    f_ax.plot (Lim_x_up,Lim_y_up,color='grey')

    x= dfdr[d]['Performance_x']
    y= dfdr[d]['Performance_y']
    z1=x/y # do opposite because of colors red
    z2=[1 if z>1 else 0 for z in z1 ]
    nn=sum(z2)/len(z2) #x= clustering y = MMSBM
    npos=1-nn
    ax1=sns.scatterplot(x,y, s=6,alpha=1,hue=z2,legend=False)
    plt.ylabel('')#, fontsize=15)
    plt.xlabel('')#, fontsize=15)
    if Counter_Plt==0:                                              
        plt.ylabel('MMSBM')#, fontsize=15)
    if Counter_Plt==2:  
        plt.xlabel('Clustering')#, fontsize=15)                        
    if Counter_Plt==4:
        plt.ylabel('RECALL',labelpad=10,fontsize=20)
        ax1.yaxis.set_label_position("right")
       ##Inset
   
    axins = inset_axes(f_ax, width="20%", height="40%", loc=4,borderpad=1)
    axi=sns.barplot([0,0.5],[npos,nn],hue=[0,1],dodge=False)
    sns.despine(trim=True,ax=axi)
    #plt.xlim(-0.5,1)
    for bar in axi.patches:
        width=0.75
        bar.set_width(width)
    axins.set_xticklabels([])
    axins.legend_.remove()
 

 
    Counter_Plt+=1
plt.tight_layout()
plt.savefig('metrics_summary.png',dpi=300)
plt.show()
