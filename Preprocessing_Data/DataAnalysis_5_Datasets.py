import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import math

#0.FUNCTIONS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#0.1. Read Data
#===============================================
def read_file(inFileName):
    lines = open(inFileName).readlines()
    d = [line.strip().split() for line in lines]
    return d
#===============================================


#0.2. Bin data for histogram
#===========================================================================
def Categories(Vector):

    #Convert to numpy array
    Vector=np.asarray(Vector)

    #Check max and min orders of magnitude
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    nonzero_min=np.min( Vector[np.nonzero(Vector)] )
    OrderOfMagnitude_Min=math.floor(math.log10( nonzero_min  ))

    nonzero_max=np.max(Vector)
    OrderOfMagnitude_Max=math.floor(math.log10( nonzero_max  ))

    Range_of_Orders= int(-OrderOfMagnitude_Min - (-OrderOfMagnitude_Max))
    print Range_of_Orders
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #Count all orders of magnitude and zeros
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    Abundances=[0]*(Range_of_Orders+1+1) #---> All orders of magnitude + 0

    for element in Vector:
        
        try:
            OrderOfMagnitude= math.floor(math.log10( element))    
            Index_in_Bins = int(OrderOfMagnitude - OrderOfMagnitude_Max - 1)
            Abundances[Index_in_Bins]+=1

        except ValueError:
            Abundances[0]+=1
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    return Abundances, int(OrderOfMagnitude_Max), int(OrderOfMagnitude_Min)

#===========================================================================

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#1.MAIN
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
path='../../Input_Data/'

Datasets=['S-8_Tot', 'V-10_Tot', 'V-22_Tot', 'V-23_24_Tot', 'V-25_Tot']
        
fig=plt.figure(figsize=(22,7))
gs = gridspec.GridSpec(2, 5)
gs.update(left=0.05, right=0.99, wspace=0.25,hspace=0.04)

column_plot=0

for Dataset in Datasets:

    #Read data and format to vector
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#    Data=read_file(path + "Leave_One_Out_" + Dataset + '/' + 'Adyacencia_' + Dataset + '.txt')
    Data=read_file(path +  Dataset + '.txt')

    Vector=[]
    print len(Data), len(Data[0])

    for line in Data:
        for column in line:
            Vector.append(float(column))

    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #Get Information about datasets
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    Abundances, OrderOfMagnitude_Max, OrderOfMagnitude_Min=Categories(Vector)
    print Abundances
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #Parameters for plot
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    bins=\
    [10**(Order) for Order in \
     range(OrderOfMagnitude_Min-1, OrderOfMagnitude_Max+3)]

    print bins
    print Abundances, sum(Abundances)
    print "zeros ratio", Abundances[0]/float(sum(Abundances))
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
    #TOP PLOT
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ax=plt.subplot(gs[0,column_plot])
    ax.set_title( Dataset + ": " + str(len(Data[0])) + ' Patients, '+ str( len(Data) ) + ' microbes',fontsize=11)
    plt.xscale('log')

    plt.hist(Vector,bins=bins)
    plt.axhline(y=Abundances[0], xmin=0, xmax=1, hold=None, color='k')
    ax.set_ylim(bottom=0.97*max(Abundances))
    
    
    #Diagonal line to show cut
    #-----------------------------------------------------------------
    d = .030  # how big to make the diagonal lines in axes coordinates
    # arguments to pass to plot, just so we don't keep repeating them
    kwargs = dict(transform=ax.transAxes, color='r', clip_on=False)
    ax.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
    ax.plot((1-d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal
    #-----------------------------------------------------------------
    
    #Fijamos xticks y posiciones
    #------------------------------------------------
    xticks=[tick for tick in bins[:-1] ]
    xticks_position=[2.6*tick for tick in bins[:-1] ]
    plt.xticks(xticks_position,([]))

    ax.xaxis.set_ticks_position('none') 
    plt.yticks(fontsize=9)
    #------------------------------------------------

    legend_x=1e-7;legend_y=Abundances[0] + 100
    ax.text(legend_x, legend_y, 'Rate of zeros: %f' %(Abundances[0]/float(sum(Abundances)) ), style='normal', bbox={'facecolor':'black', 'alpha':0.5, 'pad':5})
    

    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #BOTTOM PLOT
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ax1=plt.subplot(gs[1,column_plot])
    plt.xscale('log')
    ax1.set_ylim(0, max(Abundances[1:])*1.5 )

    plt.hist(Vector,bins=bins)
    plt.axhline(y=Abundances[0], xmin=0, xmax=1, hold=None, color='k')
    
    #Diagonal line to show cut
    #-----------------------------------------------------------------
    kwargs.update(transform=ax1.transAxes)  # switch to the bottom axes
    ax1.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
    ax1.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal
    #-----------------------------------------------------------------

    #Fijamos xticks y posiciones
    #------------------------------------------------
    xticks=[tick for tick in bins[:-1] ]
    xticks_position=[2.8*tick for tick in bins[:-1] ]
    plt.xticks(fontsize=9)
    plt.yticks(fontsize=9)
    plt.xticks(xticks_position,xticks,rotation='-30')
    #------------------------------------------------


    
    
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    column_plot+=1
    
plt.savefig("../../Plots/First_Path/Second_Histograms_" + ".pdf")
plt.show()
