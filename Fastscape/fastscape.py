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


# Scale of the grid; km
xl = 100 * 10**3
yl = 100 * 10**3
#Resolution of the grid
nx = 101
<<<<<<< HEAD
ny = 101 #number of nodes
nn = nx*ny

=======
ny = 101 #number of cells
nn = nx*ny
>>>>>>> origin/master
dx = xl/(nx-1)
dy = yl/(ny-1) 

dt = 1000 #yrs; timestep
nstep = 1000 #number of timesteps
n = 1 #Slope exponent
m = 0.4 # drainage area exponent

def plot_h(h_in):
    """Take 1-d array, plot as color mesh grid"""
    grid = np.reshape(h_in,(nx,ny))
    plt.pcolormesh(grid)


# Create array that stores elevation scalar
# Convention wherein single number describes node
h = np.random.rand((nn))


plot_h(h)

diag_dist = np.sqrt((dx**2) + (dy**2))
receiver = range(len(h))
slope = np.zeros(nn)

<<<<<<< HEAD
for ij in range(len(h)):
=======
indexVector = range(nn)


# to translate between 2D i,j matrix and
# 1D ij vector, ij = i + j*nx 
# Boundary rows:
# i = 0; i = ny-1
# j = 0; j = nx-1
# in python, j: [0,ny-1]; i: [0, nx-1]
# Need to take subset of indexVector away from boundaries

#i = 0
boundary1 = range(ny)*nx
#i = ny-1
boundary2 = (ny-1) + range(ny)*nx
#j = 0
boundary3 = range(nx)
#j = nx-1
boundary4 = range(nx) + (nx-1)*nx

#Now, create subset of indexVector away from boundaries

for ij in range(nn):
>>>>>>> origin/master
    # if not on boundary:
    ij_neighbors =[[ij-nx-1], [ij-nx], [ij-nx+1],
                   [ij-1],    [ij],    [ij+1],
                   [ij+nx+1], [ij+nx], [ij+nx+1]]  
    
    h_neighbors = h[ij_neighbors]
    
    dist_neighbors = [diag_dist, dy, diag_dist,
                      dx,        1,  dx,
                      diag_dist, dy, diag_dist]

    delta_h = h[ij] - h_neighbors    
    
    slope_neighbors = delta_h/dist_neighbors

    steepest_descent = max(slope_neighbors)
    steepest_descent_index = np.argmax(slope_neighbors)
    receiver[ij] = ij_neighbors[steepest_descent_index][0]
    slope[ij] = steepest_descent
    