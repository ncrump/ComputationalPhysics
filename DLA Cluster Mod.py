"""
Created on Sun Oct 27 17:50:45 2013
"""

# DLA Cluster Mod
"""
This generates a DLA particle cluster.
"""

import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
    

# DLA Cluster Mod
#*******************************************************************
#-------------------------------------------------------------------
# N = number of sites to occupy
# seed = start location of occupied site (integers only - ex: [0,0])
# R = radius of circle from outtermost perimeter site to start walks
#------------------------------------------------------------------- 
def DLAmod(N, seed, R):

    # define initial perimeter sites
    p1 = [seed[0]-1,seed[1]]    
    p2 = [seed[0]+1,seed[1]]     
    p3 = [seed[0],seed[1]-1] 
    p4 = [seed[0],seed[1]+1]
    
    # store occupied sites and perimeter sites
    Osites = [seed]
    Psites = [p1,p2,p3,p4]
    perim0 = [p1,p2,p3,p4]
    
    # keep track of rate that perimeter grows
    dist0 = 0
    step = [0]
    size = [0]
    
    # start while loop counter
    n = 1
    
    # loop until N sites occupied
    while n < N:
        
        # get max distance of perimeter sites for updating radius
        Pdist = [(i[0]**2 + i[1]**2)**0.5 for i in Psites]
        Pmax = int(max(np.round(Pdist)))
        
        # calculate growth rate of perimeter
        if Pmax > dist0:
            dist0 = Pmax
            size.append(Pmax)
            step.append(n)
        
        # update radius and generate random angle along circle
        # this becomes start point for random walk from outter radius into cluster
        r = Pmax + R
        theta = np.random.uniform(0,2*np.pi,1)
        
        # convert randomly generated point along circle to cartesian coords
        # this is start point for random walk from outter radius into cluster
        x0 = int(r*np.cos(theta))
        y0 = int(r*np.sin(theta))
        xyRand = [x0,y0]
        
        # set distance for loop below
        dist = r
        
        
        # step through random walk as long as particle is inside circle
        while dist <= r:
            
            # perform random walk
            # generate random uniform number between [0,1)
            rand = np.random.uniform(0,1,1)
            
            # depending on random number go right, left, up, or down by 1 unit
            if 0.00 <= rand < 0.25: xyRand[0] = xyRand[0]+1
            if 0.25 <= rand < 0.50: xyRand[0] = xyRand[0]-1
            if 0.50 <= rand < 0.75: xyRand[1] = xyRand[1]+1
            if 0.75 <= rand < 1.00: xyRand[1] = xyRand[1]-1
            
            # set distance to check if walk goes outside radius
            dist = int(round((xyRand[0]**2 + xyRand[1]**2)**0.5))
            
            # if random walk hits a perimeter site it becomes an occupied site
            if xyRand in Psites:
                Osites.append(xyRand)
                Psites.remove(xyRand)
                
                # get new perimeter sites from new occupied site
                p1 = [xyRand[0]-1,xyRand[1]]    
                p2 = [xyRand[0]+1,xyRand[1]]     
                p3 = [xyRand[0],xyRand[1]-1] 
                p4 = [xyRand[0],xyRand[1]+1]
                perims = [p1,p2,p3,p4]
                
                # check that new perimeter sites not already used
                for i in perims:
                    if i not in Psites and i not in Osites:
                        Psites.append(i)
                        perim0.append(i)
                break
            

        # increment counter only if site was occupied
        n = len(Osites)
        
    # store and return x,y locations of occupied sites
    Ox = [i[0] for i in Osites]
    Oy = [i[1] for i in Osites] 
    
    # store and return x,y locations of perimeter sites
    Px = [i[0] for i in perim0]
    Py = [i[1] for i in perim0]
    
    # returns occupied sites of the cluster and unoccupied lattice sites
    return Ox,Oy,Px,Py,step,size
#*******************************************************************


t1 = datetime.now()

#-------------------------------------------------------------------
# initial parameters
N = 1000
seed = [0,0]
R = 4

# call function
x,y,Px,Py,step,size = DLAmod(N, seed, R)

# get logs of both sides
logStep = np.log(step)
logSize = np.log(size)

# calculate growth rate
rate = logSize[len(logSize)-1]/logStep[len(logStep)-1] 
print '\n','growth rate = ', round(rate,3)

# get max point for setting axis limits
mx = [max(np.abs(x)), max(np.abs(y))]
lim = max(mx)+5

# plot DLA growth cluster
plt.plot(x,y,'b.')
plt.plot(Px,Py,'m.',label='DLA Model')
plt.xlabel('x')
plt.ylabel('y')
plt.xlim(-lim,lim)
plt.ylim(-lim,lim)
plt.annotate('N = '+str(N),fontsize=14,xy=(0.14,0.84),xycoords='figure fraction')
plt.legend()

# plot growth rate of cluster
plt.figure()
plt.subplot(2,1,1)
plt.plot(step,size,'b.', label='Growth Rate')
plt.plot(step,step**rate,'r-')
plt.xlabel('Time Step')
plt.ylabel('Radius')
plt.legend(loc=2)

# plot log-log of growth rate
plt.subplot(2,1,2)
plt.plot(logStep,logSize,'g.', label='Log-Log')
plt.plot(logStep,logStep*rate,'r-')
plt.xlabel('log Time Step')
plt.ylabel('log Radius')
#plt.annotate('$slope=0.471$',fontsize=18,xy=(0.14,0.33),xycoords='figure fraction')
plt.legend(loc=2)
#-------------------------------------------------------------------

t2 = datetime.now()
print 'runtime = ',t2-t1