#!/usr/bin/env python3
"""Compare a fishstep Sod shock tube profile with the exact Riemann solution.

Usage: python3 tests/shocktube.py shocktube.dat [output.png]

The simulation box is periodic with the high pressure state in
x_L <= x < x_R, so there are two mirrored Riemann problems: one at x_L
(high pressure on the right) and one at x_R (high pressure on the left).
The data file is written by fishstep built with -DSHOCK_TUBE=1, density
and pressure are normalized to the low pressure state.
"""

import sys

import numpy as np


def read_profile(path):
    with open(path) as f:
        header = f.readline().split()
    meta = {header[i]: float(header[i + 1]) for i in range(1, len(header) - 1, 2)}
    data = np.loadtxt(path)
    return meta, data


def exact_riemann(rho_l, u_l, p_l, rho_r, u_r, p_r, gamma, xi):
    """Exact Riemann solution (Toro, ch. 4) sampled at xi = x/t."""
    c_l = np.sqrt(gamma * p_l / rho_l)
    c_r = np.sqrt(gamma * p_r / rho_r)

    def f_k(p, rho_k, p_k, c_k):
        if p > p_k:
            a_k = 2.0 / ((gamma + 1.0) * rho_k)
            b_k = (gamma - 1.0) / (gamma + 1.0) * p_k
            f = (p - p_k) * np.sqrt(a_k / (p + b_k))
            df = np.sqrt(a_k / (b_k + p)) * (1.0 - 0.5 * (p - p_k) / (b_k + p))
        else:
            f = 2.0 * c_k / (gamma - 1.0) * ((p / p_k) ** ((gamma - 1.0) / (2.0 * gamma)) - 1.0)
            df = 1.0 / (rho_k * c_k) * (p / p_k) ** (-(gamma + 1.0) / (2.0 * gamma))
        return f, df

    # Newton iteration for the star region pressure
    p = 0.5 * (p_l + p_r)
    for _ in range(100):
        fl, dfl = f_k(p, rho_l, p_l, c_l)
        fr, dfr = f_k(p, rho_r, p_r, c_r)
        p_new = max(p - (fl + fr + u_r - u_l) / (dfl + dfr), 1.0e-12)
        if abs(p_new - p) < 1.0e-12 * p:
            p = p_new
            break
        p = p_new
    fl, _ = f_k(p, rho_l, p_l, c_l)
    fr, _ = f_k(p, rho_r, p_r, c_r)
    u_star = 0.5 * (u_l + u_r) + 0.5 * (fr - fl)
    p_star = p

    g = gamma
    rho = np.empty_like(xi)
    vel = np.empty_like(xi)
    prs = np.empty_like(xi)

    for i, s in enumerate(xi):
        if s <= u_star:
            # left of contact
            if p_star > p_l:
                s_l = u_l - c_l * np.sqrt((g + 1) / (2 * g) * p_star / p_l + (g - 1) / (2 * g))
                if s <= s_l:
                    rho[i], vel[i], prs[i] = rho_l, u_l, p_l
                else:
                    r = p_star / p_l
                    rho[i] = rho_l * (r + (g - 1) / (g + 1)) / ((g - 1) / (g + 1) * r + 1)
                    vel[i], prs[i] = u_star, p_star
            else:
                c_star_l = c_l * (p_star / p_l) ** ((g - 1) / (2 * g))
                if s <= u_l - c_l:
                    rho[i], vel[i], prs[i] = rho_l, u_l, p_l
                elif s >= u_star - c_star_l:
                    rho[i] = rho_l * (p_star / p_l) ** (1 / g)
                    vel[i], prs[i] = u_star, p_star
                else:
                    fac = 2 / (g + 1) + (g - 1) / ((g + 1) * c_l) * (u_l - s)
                    rho[i] = rho_l * fac ** (2 / (g - 1))
                    vel[i] = 2 / (g + 1) * (c_l + (g - 1) / 2 * u_l + s)
                    prs[i] = p_l * fac ** (2 * g / (g - 1))
        else:
            # right of contact
            if p_star > p_r:
                s_r = u_r + c_r * np.sqrt((g + 1) / (2 * g) * p_star / p_r + (g - 1) / (2 * g))
                if s >= s_r:
                    rho[i], vel[i], prs[i] = rho_r, u_r, p_r
                else:
                    r = p_star / p_r
                    rho[i] = rho_r * (r + (g - 1) / (g + 1)) / ((g - 1) / (g + 1) * r + 1)
                    vel[i], prs[i] = u_star, p_star
            else:
                c_star_r = c_r * (p_star / p_r) ** ((g - 1) / (2 * g))
                if s >= u_r + c_r:
                    rho[i], vel[i], prs[i] = rho_r, u_r, p_r
                elif s <= u_star + c_star_r:
                    rho[i] = rho_r * (p_star / p_r) ** (1 / g)
                    vel[i], prs[i] = u_star, p_star
                else:
                    fac = 2 / (g + 1) - (g - 1) / ((g + 1) * c_r) * (u_r - s)
                    rho[i] = rho_r * fac ** (2 / (g - 1))
                    vel[i] = 2 / (g + 1) * (-c_r + (g - 1) / 2 * u_r + s)
                    prs[i] = p_r * fac ** (2 * g / (g - 1))

    return rho, vel, prs, p_star, u_star


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "shocktube.dat"
    png = sys.argv[2] if len(sys.argv) > 2 else path.rsplit(".", 1)[0] + ".png"

    meta, data = read_profile(path)
    t, g = meta["t"], meta["gamma"]
    rho_l, p_l = meta["rho_L"], meta["P_L"]
    x_l, x_r = meta["x_L"], meta["x_R"]
    n = len(data)
    x, rho, vel, prs = data[:, 0], data[:, 1], data[:, 2], data[:, 3]

    # exact solution with rho_R = P_R = 1, velocities scale with sqrt(P_R/rho_R)
    p_unit = meta["p_over_rho_R"]
    u_left = p_l / rho_l * p_unit / (g - 1.0)

    mid = 0.5 * (x_l + x_r)
    rho_ex = np.empty(n)
    vel_ex = np.empty(n)
    prs_ex = np.empty(n)

    # problem at x_R: high pressure on the left
    right = x >= mid
    xi = (x[right] - x_r) / t / np.sqrt(p_unit)
    r_, v_, p_, p_star, u_star = exact_riemann(rho_l, 0.0, p_l, 1.0, 0.0, 1.0, g, xi)
    rho_ex[right], vel_ex[right], prs_ex[right] = r_, v_ * np.sqrt(p_unit), p_

    # problem at x_L: mirror image, high pressure on the right
    left = ~right
    xi = -(x[left] - x_l) / t / np.sqrt(p_unit)
    r_, v_, p_, _, _ = exact_riemann(rho_l, 0.0, p_l, 1.0, 0.0, 1.0, g, xi)
    rho_ex[left], vel_ex[left], prs_ex[left] = r_, -v_ * np.sqrt(p_unit), p_

    # L1 errors relative to the mean absolute exact value
    def l1(sim, ex):
        return np.mean(np.abs(sim - ex)) / np.mean(np.abs(ex))

    # shock position on the x_R side, taken where v_x falls through u*/2 behind
    # the shock, velocity is smooth in PIC while density carries particle noise
    r = p_star / 1.0
    rho_post = (r + (g - 1) / (g + 1)) / ((g - 1) / (g + 1) * r + 1)
    s_shock = np.sqrt(g * p_unit) * np.sqrt((g + 1) / (2 * g) * p_star + (g - 1) / (2 * g))
    x_shock_exact = x_r + s_shock * t
    level = 0.5 * u_star * np.sqrt(p_unit)
    idx = np.where(right & (x > x_r) & (vel > level))[0]
    x_shock_sim = np.nan
    if len(idx):
        i = idx[-1]
        x_shock_sim = x[i] + (vel[i] - level) / (vel[i] - vel[i + 1]) * (x[i + 1] - x[i])
    x_contact = x_r + u_star * np.sqrt(p_unit) * t

    print(f"t = {t:g}, gamma = {g:.4f}, rho_L/rho_R = {rho_l:g}, P_L/P_R = {p_l:g}")
    print(f"exact star state: P*/P_R = {p_star:.4f}, u* = {u_star*np.sqrt(p_unit):.4f} cells/time")
    print(f"L1 relative error: rho {l1(rho, rho_ex):.4f}  P {l1(prs, prs_ex):.4f}  "
          f"v_x {np.mean(np.abs(vel - vel_ex)) / (u_star*np.sqrt(p_unit)):.4f} (of u*)")
    print(f"shock position: simulation {x_shock_sim:.2f}, exact {x_shock_exact:.2f}, "
          f"difference {x_shock_sim - x_shock_exact:+.2f} cells")
    plateau = (x > x_contact + 0.25 * (x_shock_exact - x_contact)) & (x < x_shock_exact - 0.25 * (x_shock_exact - x_contact))
    print(f"post shock plateau density: simulation mean {np.mean(rho[plateau]):.4f} "
          f"rms scatter {np.std(rho[plateau]):.4f}, exact {rho_post:.4f}")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available, no plot written")
        return

    fig, axes = plt.subplots(2, 2, figsize=(11, 7), sharex=True)
    for ax, sim, ex, label in zip(axes.flat,
                                  (rho, vel, prs, data[:, 4] / u_left),
                                  (rho_ex, vel_ex, prs_ex, prs_ex / rho_ex * rho_l / p_l),
                                  ("density / rho_R", "v_x [cells/time]", "pressure / P_R", "u / u_L")):
        ax.plot(x, ex, color="0.2", lw=1.2, label="exact")
        ax.plot(x, sim, ".", ms=2.5, color="tab:orange", label="fishstep")
        ax.set_ylabel(label)
        ax.grid(alpha=0.3)
    for ax in axes[1]:
        ax.set_xlabel("x [cells]")
    axes[0, 0].legend()
    fig.suptitle(f"Sod shock tube, t = {t:g}, gamma = {g:.3f}")
    fig.tight_layout()
    fig.savefig(png, dpi=110)
    print(f"plot written to {png}")


if __name__ == "__main__":
    main()
