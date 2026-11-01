#!/usr/bin/env python3
"""Compare fishstep self-gravitating gas cylinder profiles with the exact
gamma = 2 polytropic equilibrium rho(r) = rho_c J0(k r), k = j01/R.

Usage: python3 tests/gascloud.py gascloud.dat [output.png]

The data file is written by fishstep built with -DGAS_CLOUD=1: one block per
output time with a header line holding the time, gas center of mass,
Lagrangian radii containing 10, 25, 50, 75 and 90% of the mass and, after
the first mesh solve, the energy budget, followed by the radial density
profile normalized to the exact central density.
"""

import sys

import numpy as np

J01 = 2.404825557695773


def bessel_j0(x):
    # series is accurate to ~1e-15 for the 0 <= x <= j01 range used here
    s, term, k = 0.0, 1.0, 0
    while abs(term) > 1e-17 or k < 3:
        s += term
        k += 1
        term *= -(x / 2.0) ** 2 / (k * k)
    return s


def bessel_j1(x):
    s, term, k = 0.0, x / 2.0, 0
    while abs(term) > 1e-17 or k < 3:
        s += term
        k += 1
        term *= -(x / 2.0) ** 2 / (k * (k + 1))
    return s


def read_blocks(path):
    blocks = []
    meta, rows = None, []
    for line in open(path):
        if line.startswith("#"):
            if meta is not None:
                blocks.append((meta, np.array(rows)))
            tok = line[1:].split()
            meta = {tok[i]: float(tok[i + 1]) for i in range(0, len(tok) - 1, 2)}
            rows = []
        elif line.strip():
            rows.append([float(v) for v in line.split()])
    if meta is not None:
        blocks.append((meta, np.array(rows)))
    return blocks


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "gascloud.dat"
    png = sys.argv[2] if len(sys.argv) > 2 else path.rsplit(".", 1)[0] + ".png"

    blocks = read_blocks(path)
    big_r = blocks[0][0]["R"]
    k = J01 / big_r
    j1_edge = bessel_j1(J01)

    def rho_exact(r):
        return np.array([bessel_j0(k * x) if x < big_r else 0.0 for x in r])

    def mass_fraction(r):
        return r * bessel_j1(k * r) / (big_r * j1_edge)

    # exact Lagrangian radii by bisection on the enclosed mass fraction
    fracs = (0.1, 0.25, 0.5, 0.75, 0.9)
    names = ("r10", "r25", "r50", "r75", "r90")
    exact_lagr = []
    for f in fracs:
        lo, hi = 0.0, big_r
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if mass_fraction(mid) < f else (lo, mid)
        exact_lagr.append(mid)

    times = np.array([m["t"] for m, _ in blocks])
    lagr = np.array([[m[n] for n in names] for m, _ in blocks])
    rel = lagr / np.array(exact_lagr) - 1.0

    print(f"polytropic cylinder R = {big_r:g} cells, {len(blocks)} outputs to t = {times[-1]:g}")
    print("Lagrangian radius   exact    t=0 error   max |error| over run   mean error over run")
    for i, n in enumerate(names):
        print(f"  {n}            {exact_lagr[i]:7.3f}   {100*rel[0, i]:+7.2f}%     {100*np.max(np.abs(rel[:, i])):6.2f}%"
              f"                {100*np.mean(rel[:, i]):+6.2f}%")

    # profile errors inside 0.9 R, central density from the innermost 4 cells
    print("time   L1 profile error (r < 0.9 R)   central density / exact   center drift [cells]")
    for meta, prof in blocks:
        r, rho = prof[:, 0], prof[:, 1]
        inner = r < 0.9 * big_r
        ex = rho_exact(r)
        l1 = np.mean(np.abs(rho[inner] - ex[inner])) / np.mean(ex[inner])
        core = r < 4.0
        rho_core = np.sum(rho[core] * r[core]) / np.sum(ex[core] * r[core])
        drift = np.hypot(meta["cm_x"] - blocks[0][0]["cm_x"], meta["cm_y"] - blocks[0][0]["cm_y"])
        print(f"{meta['t']:5.1f}  {100*l1:6.2f}%                         {rho_core:.4f}                    {drift:.4f}")

    # total energy relative to the first output with an energy budget, the
    # potential needs a mesh solve so t = 0 has none
    eb = [m for m, _ in blocks if "E_kin" in m]
    if eb:
        e0 = eb[0]["E_kin"] + eb[0]["E_th"] + eb[0]["E_pot"]
        print("time   E_kin         E_th          E_pot          (E - E0)/|E0|   FLIP loss      work residual")
        for m in eb:
            e = m["E_kin"] + m["E_th"] + m["E_pot"]
            print(f"{m['t']:5.1f}  {m['E_kin']:.5e}  {m['E_th']:.5e}  {m['E_pot']:.6e}  {(e - e0)/abs(e0):+.3e}      "
                  f"{m['flip_loss']:.2e}       {m['work_res']:+.2e}")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available, no plot written")
        return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.8))
    rr = np.linspace(0, 1.3 * big_r, 400)
    ax1.plot(rr, rho_exact(rr), color="0.2", lw=1.5, label="exact")
    picks = sorted(set([0, len(blocks) // 2, len(blocks) - 1]))
    colors = ("tab:blue", "tab:orange", "tab:red")
    for c, i in zip(colors, picks):
        meta, prof = blocks[i]
        ax1.plot(prof[:, 0], prof[:, 1], ".", ms=3, color=c, label=f"t = {meta['t']:g}")
    ax1.set_xlim(0, 1.3 * big_r)
    ax1.set_xlabel("r [cells]")
    ax1.set_ylabel("density / rho_c")
    ax1.set_title("radial density profile")
    ax1.grid(alpha=0.3)
    ax1.legend()

    for i, n in enumerate(names):
        ax2.plot(times, 100 * rel[:, i], "o-", ms=3, label=n)
    ax2.axhline(0, color="0.2", lw=1)
    ax2.set_xlabel("t")
    ax2.set_ylabel("Lagrangian radius error [%]")
    ax2.set_title("Lagrangian radii vs exact")
    ax2.grid(alpha=0.3)
    ax2.legend(ncol=5, fontsize=8)

    fig.suptitle("Self-gravitating gamma = 2 polytropic gas cylinder")
    fig.tight_layout()
    fig.savefig(png, dpi=110)
    print(f"plot written to {png}")


if __name__ == "__main__":
    main()
