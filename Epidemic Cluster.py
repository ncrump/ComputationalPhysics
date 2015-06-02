"""
Created on Sun Oct 27 17:50:45 2013
"""

# Epidemic Cluster
"""
This generates an Epidemic particle cluster.
"""

import numpy as np
import matplotlib.pyplot as plt


# Epidemic Cluster Model
#*******************************************************************
#-------------------------------------------------------------------
# N = number of sites to occupy
# seed = start location of occupied site (integers only - ex: [0,0])
# p = probability to occupy neighbor site (kill site with 1-p) 
#------------------------------------------------------------------- 
def Epidemic(N, seed, p):

    # define initial perimeter sites
    p1 = [seed[0]-1,seed[1]]    
    p2 = [seed[0]+1,seed[1]]     
    p3 = [seed[0],seed[1]-1] 
    p4 = [seed[0],seed[1]+1]
    
    # store occupied sites and perimeter sites
    Osites = [seed]
    Psites = [p1,p2,p3,p4]
    
    # keep track of rate that perimeter grows
    dist0 = 0
    step = [0]
    size = [0]
    
    # start while loop counter
    n = 1    
    
    # loop until N sites occupied
    while n < N:
        # randomly choose from perimeter sites to occupy
        rand = len(Psites)
        indx = int(np.random.randint(0,rand,1)) 
        
        # get and store new occupied site or remove it
        randGen = np.random.uniform(0,1,1)[0]
        
        if randGen >= 0 and randGen <= p:
            site = Psites.pop(indx)
            Osites.append(site)
            
            # calculate growth rate of perimeter
            dist = ((seed[0]-site[0])**2 + (seed[1]-site[1])**2)**0.5
            if dist > dist0:
                dist0 = dist
                size.append(dist)
                step.append(n)
            
            # get new perimeter sites from new occupied site
            p1 = [site[0]-1,site[1]]    
            p2 = [site[0]+1,site[1]]     
            p3 = [site[0],site[1]-1] 
            p4 = [site[0],site[1]+1]
            perims = [p1,p2,p3,p4]
    
            # check that new perimeter sites not already used
            for i in perims:
                if i not in Psites and i not in Osites:
                    Psites.append(i)
                    
            # increment counter
            n = n+1
        
        else:
            Psites.remove(Psites[indx])
        
    # store and return x,y locations of occupied sites
    Ox = [i[0] for i in Osites]
    Oy = [i[1] for i in Osites]    
    
    # returns occupied sites of the cluster
    return Ox,Oy,step,size
#*******************************************************************

# call function for Eden growth model
#-------------------------------------------------------------------
# initial parameters
N = 10000
seed = [0,0]
p = 0.65

# call function
x,y,step,size = Epidemic(N,seed,p)

# get logs of both sides
logStep = np.log(step)
logSize = np.log(size)

# calculate growth rate
rate = logSize[len(logSize)-1]/logStep[len(logStep)-1] 
print '\n','growth rate = ', round(rate,3)

# get max point for setting axis limits
mx = [max(np.abs(x)), max(np.abs(y))]
lim = max(mx)+5

# plot Eden growth cluster
plt.figure()
plt.plot(x,y,'b.', label='Epidemic Model')
plt.xlabel('x')
plt.ylabel('y')
plt.xlim(-lim,lim)
plt.ylim(-lim,lim)
plt.annotate('N = '+str(N),fontsize=14,xy=(0.14,0.84),xycoords='figure fraction')
plt.annotate('p = '+str(p),fontsize=14,xy=(0.14,0.80),xycoords='figure fraction')
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