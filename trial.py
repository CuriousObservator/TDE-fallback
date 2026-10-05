#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  2 20:28:36 2026

@author: rbarcrosspbar
"""

import numpy as np
import TDE_units as u
from forces import acceleration, specific_energy
from integrator import leapfrog_step
from LaneEmden import polytrope_sphere
from S3A import place_at_pericentre
import matplotlib.pyplot as plt


def run_tests():
    dt = 0.01 * u.t_p
    
    # Test 1: circular orbit at r_p. Period is 2*pi*t_p.
    pos = np.array([[u.r_p, 0.0, 0.0]])
    vel = np.array([[0.0, np.sqrt(u.G * u.M_bh / u.r_p), 0.0]])
    acc = acceleration(pos)
    E0 = specific_energy(pos, vel)[0]
    
    n_steps = int(round(2 * np.pi * u.t_p / dt))
    for _ in range(n_steps):
        pos, vel, acc = leapfrog_step(pos, vel, acc, dt)
    
    print("circular: final position", pos[0], "(expect ~[100, 0, 0])")
    print("circular: relative energy error", abs(specific_energy(pos, vel)[0] / E0 - 1))
    
    # Test 2: parabolic orbit from pericentre. Energy should stay ~0.
    pos = np.array([[u.r_p, 0.0, 0.0]])
    vel = np.array([[0.0, u.v_p, 0.0]])
    acc = acceleration(pos)
    for _ in range(10000):
        pos, vel, acc = leapfrog_step(pos, vel, acc, dt)
    
    E_test = specific_energy(pos, vel)[0]
    print("parabolic: E / (GM/r_p) =", E_test / (u.G * u.M_bh / u.r_p), "(expect ~0)")
    print("parabolic: distance", np.linalg.norm(pos[0]), "(should be far beyond r_p)")

if __name__ == "__main__":
    run_tests()

#%%


def fallback_from_energies(E, n_bins=40):
    """Return bin centres (days) and dM/dt (M_sun/yr) from specific energies."""
    bound = E[E < 0]
    T = 2 * np.pi * u.G * u.M_bh / (2 * np.abs(bound))**1.5   # code units
    T_days = T * u.T_UNIT_DAYS

    bins = np.logspace(np.log10(T_days.min()), np.log10(T_days.min() * 100), n_bins)
    counts, edges = np.histogram(T_days, bins=bins)
    m_particle = u.M_star / len(E)
    dMdt = counts * m_particle / np.diff(edges) * 365.25      # M_sun / yr
    centres = np.sqrt(edges[:-1] * edges[1:])
    return centres, dMdt

offsets = polytrope_sphere(5000, u.r_star, 1.5, np.random.default_rng(42))
pos, vel = place_at_pericentre(offsets)
E = specific_energy(pos, vel)

t, dMdt = fallback_from_energies(E)
plt.loglog(t, dMdt, "o-", label="n = 1.5 polytrope")

# ref = t > 3 * t.min()
# plt.loglog(t[ref], dMdt[ref][0] * (t[ref] / t[ref][0])**(-5/3), "k--", label=r"$t^{-5/3}$")

late = t > 10 * t.min()
A = np.mean(dMdt[late] * t[late]**(5/3))
plt.loglog(t[late], A * t[late]**(-5/3), "k--", label=r"$t^{-5/3}$")

plt.xlabel("time since disruption (days)")
plt.ylabel(r"$dM/dt$ ($M_\odot$/yr)")
plt.legend()
plt.show()
#%%
from TDE_2 import uniform_sphere
rng = np.random.default_rng(42)
N = 100000
stars = {
    "uniform": uniform_sphere(N, u.r_star, rng),
    "n = 1.5": polytrope_sphere(N, u.r_star, 1.5, rng),
    "n = 3":   polytrope_sphere(N, u.r_star, 3, rng),
}

plt.figure()
for label, offsets in stars.items():
    pos, vel = place_at_pericentre(offsets)
    t, dMdt = fallback_from_energies(specific_energy(pos, vel))
    plt.loglog(t, dMdt, label=label)

late = t > 10 * t.min()                      # anchored to the last curve (n = 3)
A = np.mean(dMdt[late] * t[late]**(5/3))
plt.loglog(t[late], 2 * A * t[late]**(-5/3), "k--", label=r"$t^{-5/3}$")

plt.xlabel("time since disruption (days)")
plt.ylabel(r"$dM/dt$ ($M_\odot$/yr)")
plt.legend()
plt.savefig("figures/fallback_comparison.png", dpi=200)
plt.show()