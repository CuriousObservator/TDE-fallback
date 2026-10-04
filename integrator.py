#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  2 20:28:15 2026

@author: rbarcrosspbar
"""

from forces import acceleration

def leapfrog_step(pos, vel, acc, dt):
    """One kick-drift-kick step. Returns updated pos, vel, acc."""
    vel_half = vel + 0.5 * dt * acc
    pos = pos + dt * vel_half
    acc = acceleration(pos)
    vel = vel_half + 0.5 * dt * acc
    return pos, vel, acc