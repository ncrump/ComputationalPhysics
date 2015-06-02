"""
Created on Fri Nov 1 20:15:48 2013
PHYS 613, Assignment 7 
Nick Crump
"""

# Problem 2 (EBB9): Identify Percolation Clusters
"""
This routine reads a text file containing the x,y positions of occupied sites
on a 2D square lattice to identify and label clusters for visualization and 
determining if the lattice contains a spanning (percolating) cluster.  
"""

import numpy as np
import matplotlib.pyplot as plt


# Cluster Find Routine
#*******************************************************************
# called as ClusterFind(filename, N, iters, txt)
#-------------------------------------------------------------------
# filename = name of input text file (entered as: 'filename.txt')
# N = length of square lattice (ex: N=40 for 40x40 grid)
# iters = number of correction scans for labeling (ex: iters = 2)
# txt = turn cluster numbers on/off (entered as: 'on' or 'off')
# vis = turn plotting on/off (only use for small clusters of N<100)
#------------------------------------------------------------------- 
def ClusterFind(filename, N, iters, txt, vis):
    
    # read text file containing x,y positions of occupied sites
    Ox, Oy = np.loadtxt(filename, dtype=int, skiprows=0, unpack=True) 
    Nsites = len(Ox)    # number of occupied sites
    lattice = range(N)  # length of square lattice
    
    # generate [x,y] position arrays of occupied and unoccupied sites
    Osites = [[Ox[i],Oy[i]] for i in range(Nsites)]
    Psites = [[i,j] for i in lattice for j in lattice]
    
    # define array for storing labels, initialized to integer zeros
    labels = np.zeros(Nsites) 
    labels = labels.astype(int)
    label = 1        
    
    # enter main loop to find clusters and assign cluster labels
    # ----------------------------------------------------------
    for i in range(Nsites):
        # get occupied site
        site = Osites[i]
        
        # if site has no label, assign new integer label
        if labels[i] == 0: 
            labels[i] = label
            label = label + 1  # increment label
            iLabl = label - 1  # unincremented label for checks below
        
        # get neighboring sites left, right, up, down
        lt = [site[0]-1, site[1]]
        rt = [site[0]+1, site[1]]
        up = [site[0], site[1]+1]
        dn = [site[0], site[1]-1]
        pos = [lt,rt,up,dn]
        
        # store site index and label to temp arrays for checks below       
        posIndx = [i]
        posLabl = [iLabl]
        
        # loop through neighbor sites to see if cluster sites
        for j in pos:
            if j in Osites:
                # get neighbor site index and label
                indx = Osites.index(j)
                jLabl = labels[indx]
                
                # if cluster site already labeled, store to array
                if jLabl != 0:
                    posIndx.append(indx)
                    posLabl.append(jLabl)
                
                # if cluster site not labeled, give it new label
                if jLabl == 0:
                    labels[indx] = labels[i]
        
        # get lowest label of neighbor cluster sites
        # relabel neighbor cluster sites to lowest label
        if len(posLabl) > 0:
            minLabl = min(posLabl)
            labels[(posIndx)] = minLabl
    # ----------------------------------------------------------
            
                    
    # enter correction loop to refine cluster labels
    # ----------------------------------------------------------
    # loop through lattice as many times as user input
    for chk in range(iters+1):
        for site in Osites:
            # get neighbor sites again
            lt = [site[0]-1, site[1]]
            rt = [site[0]+1, site[1]]
            up = [site[0], site[1]+1]
            dn = [site[0], site[1]-1]
            pos = [lt,rt,up,dn]
            
            # get lowest label of neighbor cluster sites
            chkIndx = [Osites.index(z) for z in pos if z in Osites]
            chkLabl = [labels[z] for z in chkIndx]
            
            # relabel neighbor cluster sites to lowest label 
            # this refines cluster identification
            if len(chkLabl) > 0:
                minLabl = min(chkLabl)
                labels[(chkIndx)] = minLabl            
    # ----------------------------------------------------------
                

    # scale labels to lowest sequential numbering of clusters
    Lsort = list(set(labels))
    mx = len(Lsort) + 1      
    for i in range(1,mx):
        indx1 = Lsort[i-1]
        indx2 = np.where(labels == indx1)[0]
        labels[indx2] = i
    
    # store x,y locations of unoccupied sites
    Px = [i[0] for i in Psites]
    Py = [i[1] for i in Psites]
        
    # enter color and plotting routine
    # ----------------------------------------------------------
    if vis == 'on':
        # define RGB color array for plotting clusters
        nColors = max(labels) + 1           
        R = np.random.uniform(0.1,1,nColors)  # random R value
        G = np.random.uniform(0.1,1,nColors)  # random G value
        B = np.random.uniform(0.1,1,nColors)  # random B value 
        colors = [[R[i],G[i],B[i]] for i in range(nColors)]
        
        # set marker size from derived visual scaling
        msize = int(128325*(N**-2.179))
        
        # loop to plot cluster sites with custom colors
        for i in range(Nsites):
            siteColor = colors[labels[i]]
            plt.scatter(Ox[i],Oy[i],s=msize,marker='s',c=siteColor,edgecolor='k')
            # plot cluster numbers if user input = 'on'
            if txt == 'on':
                plt.text(Ox[i],Oy[i],str(labels[i]),horizontalalignment='center',verticalalignment='center')
        
        # plot remaining unoccupied sites of lattice
        plt.scatter(Px,Py,s=msize,marker='s',c='none',edgecolor='k')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.xlim(-1.0,N)
        plt.ylim(-1.0,N)
    # ----------------------------------------------------------
        
    # print total clusters found
    print '\nclusters = ', len(set(labels)) 
#*******************************************************************

# sample function call
ClusterFind('PercolationCluster_40n60p.txt', 40, 2, 'off', 'on')