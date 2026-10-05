#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  2 20:30:50 2026

@author: rbarcrosspbar
"""

import numpy as np
import TDE_units as u
from forces import acceleration, specific_energy
from integrator import leapfrog_step
from LaneEmden import polytrope_sphere
from S3A import place_at_pericentre
import matplotlib.pyplot as plt
import os
  
os.makedirs("figures", exist_ok=True)




def run_direct(pos, vel, t_end, dt):
    """Integrate bound debris; return each particle's energy and return time."""
    E = specific_energy(pos, vel)
    bound = E < 0
    pos, vel, E = pos[bound], vel[bound], E[bound]

    acc = acceleration(pos)
    t_return = np.full(len(pos), np.nan)
    rv_old = np.sum(pos * vel, axis=1)
    t = 0.0

    for _ in range(int(t_end / dt)):
        pos, vel, acc = leapfrog_step(pos, vel, acc, dt)
        t += dt
        rv = np.sum(pos * vel, axis=1)
        crossed = (rv_old < 0) & (rv >= 0) & np.isnan(t_return) & (t > 10 * u.t_p)
        t_return[crossed] = t
        rv_old = rv

    return E, t_return


if __name__ == "__main__":
    
    offsets = polytrope_sphere(5000, u.r_star, 1.5, np.random.default_rng(42))
    pos, vel = place_at_pericentre(offsets)
    E, t_return = run_direct(pos, vel, t_end=3000, dt=0.01 * u.t_p)
    
    T_pred = 2 * np.pi * u.G * u.M_bh / (2 * np.abs(E))**1.5
    done = ~np.isnan(t_return)
    print("returned:", done.sum())
    print("max relative difference:", np.max(np.abs(t_return[done] / T_pred[done] - 1)))
    
    plt.loglog(T_pred[done] * u.T_UNIT_DAYS, t_return[done] * u.T_UNIT_DAYS, ".")
    plt.xlabel("predicted return time (days)")
    plt.ylabel("measured return time (days)")
    plt.savefig("figures/return_time.png", dpi=300)
    plt.show()

    
    N = 5000
    offsets = polytrope_sphere(N, u.r_star, 1.5, np.random.default_rng(42))
    pos, vel = place_at_pericentre(offsets)
    
    t_end = 20000
    E, t_return = run_direct(pos, vel, t_end=t_end, dt=0.01 * u.t_p)
    T_pred = 2 * np.pi * u.G * u.M_bh / (2 * np.abs(E))**1.5
    done = ~np.isnan(t_return)
    print("returned:", done.sum(), "of", len(E))
    
    bins = np.logspace(np.log10(45), np.log10(t_end * u.T_UNIT_DAYS), 20)
    centres = np.sqrt(bins[:-1] * bins[1:])
    m_particle = u.M_star / N
    
    plt.figure()
    for times, style, label in [(T_pred, "-", "from energies"),
                                (t_return[done], "o", "direct simulation")]:
        counts, _ = np.histogram(times * u.T_UNIT_DAYS, bins=bins)
        plt.loglog(centres, counts * m_particle / np.diff(bins) * 365.25, style, label=label)
    
    plt.xlabel("time since disruption (days)")
    plt.ylabel(r"$dM/dt$ ($M_\odot$/yr)")
    plt.legend()
    plt.savefig("figures/dM-dt_plot.png", dpi=300)
    plt.show()

    E_all = specific_energy(pos, vel)            # all particles, bound and unbound
    dE = u.G * u.M_bh * u.r_star / u.r_t**2      # 100 in code units
    
    plt.figure()
    plt.hist(E_all / dE, bins=60, density=True)
    plt.xlabel(r"$E / \Delta E$")
    plt.ylabel(r"$dM/dE$ (normalised)")
    plt.savefig("figures/dM-dE_plot.png", dpi=300)
    plt.show()