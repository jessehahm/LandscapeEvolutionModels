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
xl = 100 * 10**-2
yl = 100 * 10**-2
#Resolution of the grid
nx = 4
ny = 5 #number of nodes
nn = nx*ny
indexVector = np.arange(nn)
dx = xl/(nx)
dy = yl/(ny) 

dt = 1000 #yrs; timestep
nstep = 1000 #number of timesteps
n = 1 #Slope exponent
m = 0.4 # drainage area exponent

def plot_h(h_in):
    """Take 1-d array, plot as color mesh grid"""
    grid = np.reshape(h_in,(nx,ny))
    plt.pcolormesh(grid)
    plt.gca().invert_yaxis()


# Create array that stores elevation scalar
# Convention wherein single number describes node
#h = np.random.rand((nn))
h = np.arange(nn)

plot_h(h)

diag_dist = np.sqrt((dx**2) + (dy**2))
receiver = np.arange(nn)
slope = np.zeros(nn)

twoD_index = indexVector.reshape(ny,nx)
twoD_noBoundary = twoD_index[1:-1,1:-1]
oneD_noBoundary = twoD_noBoundary.ravel()

for ij in oneD_noBoundary:
    # if not on boundary:
    ij_neighbors =np.array([ij-nx-1, ij-nx, ij-nx+1,
                            ij-1,    ij,    ij+1,
                            ij+nx+1, ij+nx, ij+nx+1])  
    
    h_neighbors = h[ij_neighbors]
    
    dist_neighbors = np.array([diag_dist, dy, diag_dist,
                               dx,        1,  dx,
                               diag_dist, dy, diag_dist])

    delta_h = h[ij] - h_neighbors    
    
    slope_neighbors = delta_h/dist_neighbors

    steepest_descent = max(slope_neighbors)
    steepest_descent_index = np.argmax(slope_neighbors)
    receiver[ij] = ij_neighbors[steepest_descent_index]
    slope[ij] = steepest_descent
    