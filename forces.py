#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  2 20:28:00 2026

@author: rbarcrosspbar
"""

import numpy as np
import TDE_units as u

def acceleration(pos):
    """Point-mass gravity of the black hole at the origin. pos is (N, 3)."""
    r = np.linalg.norm(pos, axis=1, keepdims=True)
    return -u.G * u.M_bh * pos / r**3

def specific_energy(pos, vel):
    r = np.linalg.norm(pos, axis=1)
    return 0.5 * np.sum(vel**2, axis=1) - u.G * u.M_bh / r