#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  2 12:22:02 2026

@author: rbarcrosspbar
"""

import numpy as np
import TDE_units

rng = np.random.default_rng()
# R = float(input("Radius of star:"))
# n = int(input("Number of particles:"))
n = 10000
particle_grid = rng.uniform(low=-TDE_units.r_star , high=+TDE_units.r_star ,size=(n,3))

print(particle_grid)

#%%
ind_to_del = []
for i in range(len(particle_grid)):
    # print(np.sqrt((particle_grid[i][0]**2+particle_grid[i][1]**2+particle_grid[i][2]**2)))
    if(np.sqrt((particle_grid[i][0]**2+particle_grid[i][1]**2+particle_grid[i][2]**2))>TDE_units.r_star):
        ind_to_del.append(i)

particle_grid = np.delete(particle_grid, ind_to_del, axis=0)

print(len(particle_grid))
#%%
x = particle_grid[:,0]
y = particle_grid[:,1]
z = particle_grid[:,2]

#%%
import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter(x, y, z, cmap='viridis', marker='o')