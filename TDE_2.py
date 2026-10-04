#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  2 14:11:00 2026

@author: rbarcrosspbar
"""

import numpy as np
import matplotlib.pyplot as plt
import TDE_units

def uniform_sphere(n, radius, rng):
    points = np.empty((0,3))
    while len(points) < n:
        trial = rng.uniform(-radius, radius, size=(2*n,3))
        r = np.linalg.norm(trial, axis=1)
        points = np.vstack([points, trial[r<=radius]])
    return points[:n]

if __name__ == "__main__":
    rng = np.random.default_rng(42)
    R = TDE_units.r_star
    pos = uniform_sphere(5000, R, rng)
    print(pos.shape)
    
    r = np.linalg.norm(pos, axis=1)
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5))
    
    ax1.scatter(pos[:, 0], pos[:, 1], s=1)
    ax1.set_aspect("equal")
    ax1.set_xlabel("x")
    ax1.set_ylabel("y")
    
    ax2.hist(r, bins=40, density=True, alpha=0.6, label="sampled")
    rr = np.linspace(0, R, 200)
    ax2.plot(rr, 3*rr**2/ R**3, "k--", label=r"$3r^2/R^3$")
    ax2.set_xlabel("r")
    ax2.set_ylabel("probability density")
    ax2.legend()
    
    plt.tight_layout()
    plt.show()
    
    
    
    
    