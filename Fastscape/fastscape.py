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

def plot_h(h_in):
    """Take 1-d array, plot as color mesh grid"""
    grid = np.reshape(h_in,(nx,ny))
    plt.pcolormesh(grid)

# Scale of the grid; km
xl = 100 * 10**3
yl = 100 * 10**3
nn = nx*ny
#Resolution of the grid
nx = 101
ny = 101 #number of cells
dx = xl/(nx-1)
dy = yl/(ny-1) 

dt = 1000 #yrs; timestep
nstep = 1000 #number of timesteps
n = 1 #Slope exponent
m = 0.4 # drainage area exponent

# Create array that stores elevation scalar
# Convention wherein single number describes node
h = np.random.rand((nn))


plot_h(h)

diag_dist = np.sqrt((dx**2) + (dy**2))
receiver = np.zeros(nn)
slope = np.zeros(nn)
for ij in range(len(h)):
    # if not on boundary:
    ij_neighbors =[[ij-nx-1], [ij-nx], [ij-nx+1],
                   [ij-1],    [ij],    [ij+1],
                   [ij+nx+1], [ij+nx], [ij+nx+1]]  
    
    h_neighbors = [h[ij-nx-1], h[ij-nx], h[ij-nx+1],
                   h[ij-1],    h[ij],    h[ij+1],
                   h[ij+nx+1], h[ij+nx], h[ij+nx+1]]
                   
    dist_neighbors = [diag_dist, dy, diag_dist,
                      dx,        1,  dx,
                      diag_dist, dy, diag_dist]

    delta_h = h[ij] - h_neighbors    
    
    slope_neighbors = delta_h/dist_neighbors

    steepest_descent = max(slope_neighbors)
    steepest_descent_index = np.argmax(slope_neighbors)
    receiver[ij] = ij_neighbors[steepest_descent_index][0]
    slope[ij] = steepest_descent
    