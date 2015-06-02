"""
Created on Mon Oct 14 22:02:21 2013
PHYS 613, Project 2
Nick Crump
"""

# Project 2: Genetic Algorithm
"""
Writing a genetic algorithm to find the global minimum of the
Lennard-Jones interatomic potential function of a molecule. 
"""

import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt


# define function for Lennard-Jones potential
# returns value of multivariable function of 3N variables for N atoms    
#*******************************************************************
# ------------------------------------------------------------------
# input argument 'pos0' is an array of initial positions 
# 'pos0' is input as [x1,y1,z1,x2,y2,z2 ... xn,yn,zn]
# ------------------------------------------------------------------
def LennardJones(pos0):
    
    # get number of atoms for loop from input positions
    n = int(len(pos0)/3)
    
    # term for iteration
    potential = 0
    
    # build array of initial positions of atoms in (x0,y0,z0)
    x0 = [pos0[k] for k in range(0, 3*n, 3)]
    y0 = [pos0[k] for k in range(1, 3*n, 3)]
    z0 = [pos0[k] for k in range(2, 3*n, 3)]
    
    # iterate through double sum of energy function
    # want i from i = 1 to n-1 atoms
    for i in range(n-1):
        # get initial xi positions
        xi,yi,zi = x0[i],y0[i],z0[i]
        
        # want j > i atoms
        for j in range(i+1,n):
            # get initial xj positions
            xj,yj,zj = x0[j],y0[j],z0[j]  
            
            # construct 3-D distance equation
            rij = ((xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2)**0.5

            
            # calculate terms of Lennard-Jones potential
            potential = potential + (rij**-12.0) - (rij**-6.0) 
            
    # return value for potential energy
    return potential
#******************************************************************* 


## use scipy.optimize to find min of function using conjugate gradient 
##*******************************************************************
#Natoms = 5
#Nrandom = 10
#Pmin = 0
#Pvar = 0
#for i in range(Nrandom):
#    pos0 = np.random.randint(0,10,3*Natoms)
#    minPos,minFun,fc,gc,wf = optimize.fmin_cg(LennardJones,pos0,full_output=1,disp=0)
#    if minFun < Pmin:
#        Pmin = minFun
#        Pvar = minPos
#print '\n','conjugate gradient minimum = ', round(Pmin,3)
##*******************************************************************


# performs genetic algorithm to find minimum of Lennard Jones function
#*******************************************************************
# ------------------------------------------------------------------
# nAtoms = number of atoms in simulation (# of coordinates = 3*nAtoms)
# sd = side length of bounding box to generate random atom positions
# nPop = initial population size (# of arrays to initialize)
# Pcross = crossover probability (percentage of population selected for mating)
# Pmute = mutation probability (percentage of children selected for mutation)
# Pperm = permutation probability (percentage of children selected for permutation)
# nIters = max number of iterations to run
# ------------------------------------------------------------------
def GA(nAtoms, sd, nPop, Pcross, Pmute, Pperm, nIters):
    # ------------------------------------------------------------------
    # define initial parameters
    # ------------------------------------------------------------------
    iters = 0             # number of iterations
    gMinFun = 0           # global minimum found
    gMinVar = 0           # coordinates at global minimum
    gMinNum = 0           # iteration number of global minimum
    
    nCoords = 3*nAtoms    # number of coordinates in each individual
    nVals = nCoords*nPop  # number of values to generate randomly
    fitAv = []            # stores average fitness of each population
    # ------------------------------------------------------------------
    
    
    # ------------------------------------------------------------------
    # define initial population
    # ------------------------------------------------------------------
    pop = np.random.uniform(0,sd,nVals)   # randomly generate x,y,z values
    pop = pop.reshape(nPop,nCoords)       # make array of individuals
    # ------------------------------------------------------------------
    
    
    while iters < nIters: 
        # ------------------------------------------------------------------
        # determine population fitness
        # ------------------------------------------------------------------
        # evaluate potential function for each individual in population
        fVals = [LennardJones(i) for i in pop]  
        # only keep negative potential values to evaluate fitness
        fits = [i for i in fVals if i < 0]
        # determine fitness rankings from 0 to 1 (1 is best)
        fmax = max(fits)
        fmin = min(fits)
        fitness = (fmax-fits)/(fmax-fmin)
        # keep average fitness of each generation to compare
        fitAv.append(round(np.average(fitness),3))
        # ------------------------------------------------------------------

                
        # ------------------------------------------------------------------
        # determine population selection from fitness/function values
        # ------------------------------------------------------------------
        # sort potential values and keep best individuals to make new generation
        fSort = np.sort(fVals)
        # get number of individuals for mating from input crossover probability
        # if number selected for mating is less than 10 print error and break
        nCross = int(nPop*Pcross)
        if nCross < 10:
            print '\n','WARNING: POPULATION CRITICALLY LOW'
            print 'Adjust Population Size or Crossover Probability','\n'
            iters = iters+1
            break
        # get index values of best individuals for mating
        fitIndx = [fVals.index(i) for i in fSort[0:nCross]]
        # use index values to get best individuals from population
        pop = np.array(pop)
        pop = pop[fitIndx]  # coords of minimum end up at top of array
        # ------------------------------------------------------------------

       
        # ------------------------------------------------------------------
        # create new generation by GA operators
        # ------------------------------------------------------------------
        
        # ------------------------------------------------------------------
        # two-point CROSSOVER at fixed points based on input probability
        # ------------------------------------------------------------------
        popNew = []  # stores new population
        
        # get number of children each pair must have to maintain population size
        # actually half since each pair produces 2 children
        nChild = int((nPop/nCross)/2.0)
        if nChild < 1: nChild = 1
        
        # split parents at 0.25 and 0.75 for crossover
        pt1 = int(nCoords*0.25)
        pt2 = int(nCoords*0.75)
        
        # this mates selected parents to produce next generation
        for i in range(nCross):
            # list of selected parents excluding parent i
            mates = [k for k in range(nCross) if k != i]
            
            # randomly generate second parent for selection
            for j in range(nChild):
                mate = mates.pop(np.random.randint(0,len(mates),1))
                
                # select candidate parents for mating
                parent1 = list(pop[i])
                parent2 = list(pop[mate])        
                
                # create 2 childred from mating crossover of parents selected
                child1 = parent1[0:pt1] + parent2[pt1:pt2] + parent1[pt2:nCoords]
                child2 = parent2[0:pt1] + parent1[pt1:pt2] + parent2[pt2:nCoords]
                
                # append children to new population
                popNew.append(child1)
                popNew.append(child2)
                    
        # rename new population for use in loop
        pop = popNew  
        # ------------------------------------------------------------------

        
        # ------------------------------------------------------------------
        # create random MUTATIONS on children based on input probability
        # ------------------------------------------------------------------
        # get number of children for random mutations
        # if number selected for mutation is less than 1 print error and break
        nMute = int(nPop*Pmute)
        if nMute < 1:
            print '\n','WARNING: MUTATION PROBABILITY TOO LOW'
            print 'Adjust Population Size or Mutation Probability','\n'
            iters = iters+1
            break
        # select and mutate 
        for i in range(nMute):
            mutChild = int(np.random.randint(0,nPop,1))     # random selection of child
            mutCoord = int(np.random.randint(0,nCoords,1))  # random selection of gene
            mut = float(np.random.uniform(0,sd,1))          # random mutation
            
            pop[mutChild][mutCoord] = mut
        # ------------------------------------------------------------------
            
            
        # ------------------------------------------------------------------
        # create random double PERMUTATION on children based on input probability
        # ------------------------------------------------------------------
        # get number of children for random permutation
        # if number selected for permation is less than 2 print error and break
        nPerm = int(nPop*Pperm)
        if nMute < 2:
            print '\n','WARNING: PERMUTATION PROBABILITY TOO LOW'
            print 'Adjust Population Size or Permutation Probability','\n'
            iters = iters+1
            break
        # select and permutate
        for i in range(nPerm):
            permChild = np.random.randint(0,nPop,2)      # random selection of children
            permCoord = np.random.randint(0,nCoords,2)   # random selection of genes  
            
            # get children and genes for permutation
            child1 = int(permChild[0])
            child2 = int(permChild[1])
            gene1 = int(permCoord[0])
            gene2 = int(permCoord[1])
            permGenes = [pop[child1][gene1], pop[child1][gene2], pop[child2][gene1], pop[child2][gene2]]
            
            # swap genes of children 1 and 2
            pop[child1][gene1] = permGenes[3]
            pop[child1][gene2] = permGenes[2]
            pop[child2][gene1] = permGenes[1]
            pop[child2][gene2] = permGenes[0]
        # ------------------------------------------------------------------
            
        
        # ------------------------------------------------------------------
        # compute minimum to store and iterate 
        # ------------------------------------------------------------------ 
        # gets minimum value
        gMin0 = fSort[0] 
        
        # store global minimum if less than previously found minimum
        if gMin0 < gMinFun: 
            gMinFun = gMin0   # get minimum value
            gMinVar = pop[0]  # get coordinates at minimum
            gMinNum = iters   # get iteration number at minimum
        
        # iterate counter
        iters = iters+1
        # ------------------------------------------------------------------
        
            
    return gMinFun, gMinVar, gMinNum, iters, fitAv
#*******************************************************************

# define initial parameters
Natoms = 5      # number of atoms in simulation
Bbox = 10       # side length of bounding box for initial atom positions
Npop = 100      # size of initial population
Pcross = 0.25   # crossover probability
Pmute = 0.10    # mutation probability
Pperm = 0.10    # permutation probability
Nmax = 1000     # max generations to run

# call genetic algorithm above
gMin, coords, n, iters, fitAv = GA(Natoms,Bbox,Npop,Pcross,Pmute,Pperm,Nmax)

# print global minimum and generation number
print 'genetic algorithm minimum = ', round(gMin,3)
print 'iterations to convergence = ', n,'\n'

# plot trend of average fitness vs generation number
plt.plot(range(iters),fitAv)
plt.vlines(fitAv.index(max(fitAv)),0,1,color='g',linewidth=1.5,label='Maximum Fitness')
plt.vlines(n,0,1,color='r',linewidth=1.5,label='Convergence Point')
plt.xlabel('Generation',fontsize=14)
plt.ylabel('Average Fitness',fontsize=14)
plt.ylim(0,1)
plt.legend(loc=8)
