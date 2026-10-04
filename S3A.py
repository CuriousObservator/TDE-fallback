#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  2 20:11:07 2026

@author: rbarcrosspbar
"""

import numpy as np
import TDE_units as u
from LaneEmden import polytrope_sphere

def place_at_pericentre(offsets):
    """Shift star to pericentre; give every particle the centre-of-mass velocity."""
    pos = offsets + np.array([u.r_p, 0.0, 0.0])
    vel = np.tile([0.0, u.v_p, 0.0], (len(offsets), 1))
    return pos, vel

if __name__ == "__main__":
    offsets = polytrope_sphere(5000, u.r_star, 1.5, np.random.default_rng(42))
    pos, vel = place_at_pericentre(offsets)
    r = np.linalg.norm(pos, axis=1)
    E = 0.5 * np.sum(vel**2, axis=1) - u.G * u.M_bh / r

    dE = u.G * u.M_bh * u.r_star / u.r_t**2
    print("bound fraction:", np.mean(E < 0))        # expect ~0.5
    print("E range:", E.min(), E.max())             # expect about -dE to +dE
    print("dE estimate:", dE)