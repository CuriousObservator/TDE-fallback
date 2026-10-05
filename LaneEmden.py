#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  2 15:51:07 2026

@author: rbarcrosspbar
"""

import numpy as np
import matplotlib.pyplot as plt



def f(x, y, n):
    theta, phi = y
    theta = max(theta, 0.0)          # avoids NaN from negative**1.5
    return np.array([phi, -theta**n - 2*phi/x])


def solve_lane_emden(n, h=1e-3, x0=1e-4):
    y0 = np.array([1 - x0**2/6, -x0/3])
    xval, theta, phi = [x0], [y0[0]], [y0[1]]

    while True:
        k1 = f(x0, y0, n)
        k2 = f(x0 + h/2, y0 + k1*h/2, n)
        k3 = f(x0 + h/2, y0 + k2*h/2, n)
        k4 = f(x0 + h, y0 + k3*h, n)
        y = y0 + (k1 + 2*k2 + 2*k3 + k4)*h/6

        if y[0] <= 0:
            break

        x0 += h
        y0 = y
        xval.append(x0); theta.append(y0[0]); phi.append(y0[1])

    return np.array(xval), np.array(theta), np.array(phi)

#%%
if __name__ == "__main__":
    # Test against the analytic n = 1 solution first
    xi, theta, phi = solve_lane_emden(1)
    print("n = 1 surface:", xi[-1], "(expect pi)")
    plt.figure()
    plt.plot(xi, theta, label="RK4, n = 1")
    plt.plot(xi, np.sin(xi)/xi, "k--", label=r"$\sin\xi/\xi$")

    xi, theta, phi = solve_lane_emden(1.5)
    print("n = 1.5 surface:", xi[-1], "(expect ~3.654)")
    print("-xi^2 phi at surface:", -xi[-1]**2 * phi[-1], "(expect ~2.714)")
    plt.plot(xi, theta, label="RK4, n = 1.5")

    plt.xlabel(r"$\xi$")
    plt.ylabel(r"$\theta$")
    plt.legend()
    plt.show()
#%%
def polytrope_sphere(n_particles, radius, n_poly, rng):
    """Return an (N, 3) array of points following a polytrope density profile."""
    xi, theta, phi = solve_lane_emden(n_poly)

    m_frac = -xi**2 * phi
    m_frac /= m_frac[-1]                      # 0 at centre, 1 at surface

    # Radii: invert the cumulative mass distribution
    u = rng.uniform(0, 1, n_particles)
    r = np.interp(u, m_frac, xi) / xi[-1] * radius

    # Directions: uniform on the sphere
    cos_t = rng.uniform(-1, 1, n_particles)
    sin_t = np.sqrt(1 - cos_t**2)
    az = rng.uniform(0, 2*np.pi, n_particles)

    return np.column_stack([r * sin_t * np.cos(az),
                            r * sin_t * np.sin(az),
                            r * cos_t])  
#%%
if __name__ == "__main__":
    import TDE_units

    rng = np.random.default_rng(42)
    R = TDE_units.r_star
    pos = polytrope_sphere(5000, R, 1.5, rng)
    r = np.linalg.norm(pos, axis=1)

    xi, theta, phi = solve_lane_emden(1.5)
    expected = xi**2 * theta**1.5
    rr = xi / xi[-1] * R
    expected /= np.sum(expected) * (rr[1] - rr[0])   # normalise to unit area

    plt.figure()
    plt.hist(r, bins=40, density=True, alpha=0.6, label="sampled")
    plt.plot(rr, expected, "k--", label=r"$\xi^2\theta^n$")
    plt.xlabel("r")
    plt.legend()
    plt.show()
#%%

