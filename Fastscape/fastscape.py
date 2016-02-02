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
"""
####################################
#USER DEFINED LANDSCAPE VARIABLES

# Scale of the grid; km
xl = 10**1
yl = 10**1
#Resolution of the grid
nx = 10**1
ny = 10**1 #number of nodes
nn = nx*ny

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
              0,0,0,0,6,6,6,5,4,3])
#################



dt = 1000 #yrs; timestep
nstep = 1000 #number of timesteps
n = 1 #Slope exponent
m = 0.4 # drainage area exponent
################

# DERIVED VARIABLES
dx = xl/(nx)
dy = yl/(ny) 
indexVector = np.arange(nn)
reshaped_index = indexVector.reshape(ny,nx)
def plot_mesh(oneD_in):
    """Take 1-d array, plot as color mesh grid"""
    grid = np.reshape(oneD_in,(ny,nx))
    plt.pcolormesh(grid, cmap='RdBu')
    #plt.gca().invert_yaxis()
    #plt.show()

reshaped_h = h.reshape(ny,nx)


diag_dist = np.sqrt((dx**2) + (dy**2))
receiver = np.arange(nn)
slope = np.zeros(nn)
direction = np.full(nn, 4)
twoD_index = indexVector.reshape(ny,nx)
twoD_noBoundary = twoD_index[1:-1,1:-1]
oneD_noBoundary = twoD_noBoundary.ravel()
# to get boundaries; not this!

dist_neighbors = np.array([diag_dist, dy, diag_dist,
                           dx,        1,  dx,
                           diag_dist, dy, diag_dist])
#Build receiver array
for ij in oneD_noBoundary:
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
    
reshaped_receiver = receiver.reshape(ny,nx)

## ndon = total number of donors to a node
ndon = np.zeros(nn,int)
# donors = indices of donors to a node 
donor = np.full([8,nn],-1,int)

#Build donor and ndon array
for ij in indexVector:
    if receiver[ij] != ij:   #if not myself
        recij = receiver[ij]
        donor[ndon[recij], recij] = ij
        #Increment number of donors by one for this receiver
        ndon[recij] = ndon[recij] + 1 

reshaped_ndon = ndon.reshape(ny,nx)

print 'h'
print reshaped_h
print 'index'
print reshaped_index
print 'receiver'
print reshaped_receiver
print 'ndon'
print reshaped_ndon

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

plot_mesh(h)
qx = np.arange(nx)
qy = np.arange(ny)
qU = U.reshape(ny,nx)
qV = V.reshape(ny,nx)
Q = plt.quiver(qx+(dx/2.0),qy+(dy/2.0),qU,qV)

plt.gca().invert_yaxis()
plt.show()


######## MAKE THE STACK

# START WITH BASE-LEVELS / PITS 
# (Where receiver[ij] = ij)

baseLevels = receiver[receiver == indexVector]
stack = np.empty(nn,int)
i = 0
def add_to_stack(ij):
    global i
    stack[i] = ij
    i = i+1
    for k in donor[0:ndon[ij],ij]:
              add_to_stack(k)
abc = 0

for ij in baseLevels:
    add_to_stack(ij)
    

    

    