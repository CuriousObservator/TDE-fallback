#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  1 23:53:15 2026

@author: rbarcrosspbar
"""

import numpy as np

# --- Code units: G = 1, mass in M_sun, length in R_sun ---
G = 1.0
M_star = 1.0
r_star = 1.0
M_bh = 1e6
beta = 1.0

# --- Derived scales (code units) ---
r_t = r_star * (M_bh / M_star)**(1/3)      # tidal radius
r_p = r_t / beta                            # pericentre distance
t_p = np.sqrt(r_p**3 / (G * M_bh))          # pericentre passage time
v_p = np.sqrt(2 * G * M_bh / r_p)           # parabolic speed at pericentre

# --- Conversions back to physical units ---
G_SI = 6.6743e-11        # m^3 kg^-1 s^-2
M_SUN_KG = 1.98841e30
R_SUN_M = 6.957e8

T_UNIT_S = np.sqrt(R_SUN_M**3 / (G_SI * M_SUN_KG))   # ~1593 s
T_UNIT_DAYS = T_UNIT_S / 86400
V_UNIT_KMS = R_SUN_M / T_UNIT_S / 1e3                # ~437 km/s