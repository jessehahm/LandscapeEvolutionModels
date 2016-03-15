# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 14:30:17 2016

@author: wjh

Under guidance of Jean Braun/Bill Dietrich
EPS 117, UC Berkeley, Spring 2016
"""
###  Import required libraries
import numpy as np
import matplotlib.pyplot as plt
import time
#from numba import autojit
start = time.time()


####################################
#USER DEFINED LANDSCAPE VARIABLES

# Scale of the grid; km
xl = 10**2
yl = 10**2
#Resolution of the grid
nx = 10**2
ny = 10**2 #number of nodes
nn = nx*ny
printFreq = 20

# Create array that stores elevation scalar
# Convention wherein single number describes node
h = np.random.rand((nn))
#h = np.arange(nn)
"""
################
#  TEST ARRAY  #
################
xl = 10.
yl = 5.
nx = 10
ny = 5
nn = nx*ny
h = np.array([9,0,0,0,6,6,6,5,4,3,
              2,2,2,2,5,5,5,4,4,2,
              3,3,3,3,5,4,3,2,1,0,
              2,2,2,2,5,5,5,4,4,2,
              0,0,0,0,6,6,6,5,4,3], dtype=float)
#################
"""


k = 1.0*10**(-4)
U = 2.0*10**(-3) #(m/yr)
delta_t = 10.0**3 #yrs; timestep
num_timesteps = 300 #number of timesteps
n = 1.0 #Slope exponent
m = 0.4 # drainage area exponent
################

# DERIVED VARIABLES
dx = xl/(nx)
dy = yl/(ny) 
indexVector = np.arange(nn)
reshaped_index = indexVector.reshape(ny,nx)

def plot_mesh(oneD_in, cmap = 'gist_earth'):
    """Take 1-d array, plot as color mesh grid
    optional colormap"""
    grid = np.reshape(oneD_in,(ny,nx))
    plt.pcolormesh(grid, cmap=cmap)
    plt.colorbar()
    #plt.gca().invert_yaxis()
    #plt.show()


def long_profile():
    #get index of node with biggest drainage area
    node = np.argmax(area[indexVector])
    
    #declare python lists for length, elevations
    profile_distance = [0]
    profile_h = [h[node]]
    
    #ask donors who has largest area...
    #get the donor with the biggest drainage area... 
    while (ndon[node] > 0):
        node = donor[np.argmax(area[donor[0:ndon[node], node]]), node] 
        profile_h.append(h[node])
        profile_distance.append(delta_x[node]+profile_distance[-1])
    plt.plot(profile_distance, profile_h)
    plt.show()
     
reshaped_h = h.reshape(ny,nx)


diag_dist = np.sqrt((dx**2) + (dy**2))
receiver = np.arange(nn)
slope = np.zeros(nn)
delta_x = np.zeros(nn)
direction = np.full(nn, 4)
twoD_index = indexVector.reshape(ny,nx)
twoD_noBoundary = twoD_index[1:-1,1:-1]
oneD_noBoundary = twoD_noBoundary.ravel()
oneD_Boundary = np.delete(indexVector,oneD_noBoundary)
# to get boundaries; not this!



#%% The loop
# Add U
# Then erode
# 1: Build receiver array
# 2: Build donor array

#Loop over timesteps
for istep in range(num_timesteps):    
    if istep > 150:
        k = 0.5*10**(-4) #(m/yr)

    #Build receiver array
# For cyclic: ---> include as your neighbor the cell on opposite side

    #receiver array stores each node's lowest neighbor
    for ij in oneD_noBoundary:
        dist_neighbors = np.array([diag_dist, dy, diag_dist,
                           dx,        1,  dx,
                           diag_dist, dy, diag_dist])

        # if not on boundary:
        ij_neighbors =np.array([ij-nx-1, ij-nx, ij-nx+1,
                                ij-1,    ij,    ij+1,
                                ij+nx-1, ij+nx, ij+nx+1])  
        
        h_neighbors = h[ij_neighbors]
    
        delta_h = h[ij] - h_neighbors    
        
        slope_neighbors = delta_h/dist_neighbors
    
        steepest_descent = max(slope_neighbors)
        steepest_descent_index = np.argmax(slope_neighbors)
        receiver[ij] = ij_neighbors[steepest_descent_index]
        slope[ij] = steepest_descent
        direction[ij] = steepest_descent_index    
        delta_x[ij] = dist_neighbors[steepest_descent_index]
   #Receiver array for bottom wall
    dist_neighbors = np.array([dx,        1,  dx,
                           diag_dist, dy, diag_dist])
    for ij in np.arange(1,nx-1):
        ij_neighbors =np.array([ij-1,    ij,    ij+1,
                                ij+nx-1, ij+nx, ij+nx+1])  
        
        h_neighbors = h[ij_neighbors]
    
        delta_h = h[ij] - h_neighbors    
        
        slope_neighbors = delta_h/dist_neighbors
    
        steepest_descent = max(slope_neighbors)
        steepest_descent_index = np.argmax(slope_neighbors)
        receiver[ij] = ij_neighbors[steepest_descent_index]
        slope[ij] = steepest_descent
        direction[ij] = steepest_descent_index    
        delta_x[ij] = dist_neighbors[steepest_descent_index]
 
    
    reshaped_receiver = receiver.reshape(ny,nx)
    
    
    ## ndon = total number of donors to a node
    ndon = np.zeros(nn,int)
    # donors = indices of donors to a node 
    donor = np.full([8,nn],-1,int)
    
    #Build donor and ndon array
    #donor array: list of nodes that drain to you
    for ij in indexVector:
        if receiver[ij] != ij:   #if not myself
            recij = receiver[ij]
            donor[ndon[recij], recij] = ij
            #Increment number of donors by one for this receiver
            ndon[recij] = ndon[recij] + 1 
    
    reshaped_ndon = ndon.reshape(ny,nx)
    
   
    
    #%% ####### MAKE THE STACK
    
    # START WITH BASE-LEVELS / PITS 
    # (Where receiver[ij] = ij)
    
    baseLevels = receiver[receiver == indexVector]
    stack = np.empty(nn,int)
    
    #Create array of length nn, 'catchment'
    #where value = name of catchment
    #which is actually the baselevel!
    catchment = np.empty(nn,int)
    nstack = 0
    
    def add_to_stack(ij, cc):
        global catchment 
        global nstack
        stack[nstack] = ij
        catchment[ij] = cc    
        nstack = nstack+1
        for k in donor[0:ndon[ij],ij]:
                  add_to_stack(k, cc)
    
    for ij in baseLevels:
        cc = ij
        add_to_stack(ij, cc)
        
    
    #%% Drainage area
    area = np.ones(nn)*dx*dy
    reversed_stack = stack[::-1]
    
    for ij in reversed_stack:
        if receiver[ij] != ij:
            area[receiver[ij]] = area[receiver[ij]] + area[ij] 
    
 #   print 'area:'
#    print area.reshape(ny,nx)    
    
    #%% Calculate new heights
    
    #Two exceptional cases
    # Boundaries: Don't add uplift, don't change h
    # Local minima: Just add uplift
    
    #Add uplift to all non-boundary nodes
    #for ij in oneD_noBoundary:
      #  h[ij] = h[ij] + U*delta_t
    
    #Add uplift to all nodes away from boundary
    for ij in oneD_noBoundary:
        h[ij] = h[ij] + U*delta_t
    #One-sided boundary condition (reflective):
    for ij in np.arange(nx):
        h[ij] = h[ij] + U*delta_t
    
        
    #calculate erosion for all in stack
    for ij in stack:
        if (receiver[ij] != ij):
            C = k*area[ij]**m*delta_t/delta_x[ij]
            h[ij] = (h[ij] + C*h[receiver[ij]])/(1 + C) 
            #note, for stability, C should be between 0 and 1
            #where stability = erosion not too fast or too slow

    #Print the landscape
    if (istep % printFreq == 0):
        plot_mesh(h)
        plt.show()
        long_profile()
        print 'timestep: ' + str(istep)
#%% Plotting
print 'It took', time.time()-start, 'seconds.'


# Plot grid-coded by height 
plot_mesh(h)

# Plot arrows of steepest descent
#Build arrow index vector arrays
U = np.zeros(nn)
V = np.zeros(nn)
for ij in indexVector:
    if direction[ij] == 0:
        U[ij] = -1
        V[ij] = 1
    if direction[ij] == 1:
        U[ij] = 0
        V[ij] = 1
    if direction[ij] == 2:
        U[ij] = 1
        V[ij] = 1
    if direction[ij] == 3:
        U[ij] = -1
        V[ij] = 0
    if direction[ij] == 4:
        U[ij] = 0
        V[ij] = 0
    if direction[ij] == 5:
        U[ij] = 1
        V[ij] = 0
    if direction[ij] == 6:
        U[ij] = -1
        V[ij] = -1
    if direction[ij] == 7:
        U[ij] = 0
        V[ij] = -1
    if direction[ij] == 8:
        U[ij] = 1
        V[ij] = -1



qx = np.arange(nx)
qy = np.arange(ny)
qU = U.reshape(ny,nx)
qV = V.reshape(ny,nx)
Q = plt.quiver(qx+(dx/2.0),qy+(dy/2.0),qU,qV)

plt.gca().invert_yaxis()
plt.show()
plot_mesh(catchment)
plt.show()

