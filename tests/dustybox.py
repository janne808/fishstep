#!/usr/bin/env python3
"""Check the fishstep gas-dust drag against the DUSTYBOX solution.

Usage: python3 tests/dustybox.py dustybox.dat [dustybox2.dat ...] [output.png]

Each data file is written by the DUSTY_BOX setup: uniform gas at rest and
uniform dust streaming through it with v0 along x (Laibe and Price 2011). Only
the drag acts, so the center of mass velocity stays at v0 eps/(1 + eps) and
the relative velocity decays as

    v_gas - v_dust = -v0 exp(-(1 + eps) t/t_s)

exactly for any time step, as the drag is solved analytically. The kinetic
energy lost heats the gas. The test fails if the relative velocity deviates
from the solution by more than TOL_V of v0, the momentum by more than TOL_P
or the total energy by more than TOL_E of the initial kinetic energy.
"""

import sys

import numpy as np

TOL_V = 1.0e-6
TOL_P = 1.0e-10
TOL_E = 1.0e-8


def read(path):
    meta = {}
    for line in open(path):
        if line.startswith("#"):
            tok = line[1:].split()
            if tok and tok[0] == "eps":
                meta = {tok[i]: float(tok[i + 1]) for i in range(0, len(tok) - 1, 2)}
    return meta, np.loadtxt(path)


def main():
    args = sys.argv[1:] or ["dustybox.dat"]
    png = args.pop() if args[-1].endswith(".png") else args[0].rsplit(".", 1)[0] + ".png"
    ok = True
    runs = []
    print("    t_s      dt/tau    eps   max |dw - exact|/v0   momentum drift   energy drift   heat residual")
    for path in args:
        meta, d = read(path)
        eps, ts, v0 = meta["eps"], meta["t_s"], meta["v0"]
        t, vg, vd, px, ek, eth, res = d[:, 0], d[:, 1], d[:, 3], d[:, 5], d[:, 6], d[:, 7], d[:, 8]
        tau = ts / (1.0 + eps)
        w = vg - vd
        exact = -v0 * np.exp(-t / tau)
        err_v = np.max(np.abs(w - exact)) / v0
        err_p = np.max(np.abs(px - px[0])) / abs(px[0])
        e_tot = ek + eth
        err_e = np.max(np.abs(e_tot - e_tot[0])) / ek[0]
        dt = meta["dt"]
        print(f"{ts:9.4g}  {dt / tau:9.3g}  {eps:5.2g}   {err_v:14.3e}      {err_p:12.3e}   {err_e:12.3e}   {res[-1]:12.3e}")
        ok &= err_v < TOL_V and err_p < TOL_P and err_e < TOL_E
        runs.append((meta, t, w, exact, tau))

    print("PASS" if ok else "FAIL")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available, no plot written")
        return 0 if ok else 1

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for meta, t, w, exact, tau in runs:
        label = f"t_s = {meta['t_s']:g}, dt/tau = {meta['dt'] / tau:.3g}"
        line, = ax[0].plot(t, -w / meta["v0"], "o", ms=3, label=label)
        tt = np.linspace(0, t[-1], 400)
        ax[0].plot(tt, np.exp(-tt / tau), "-", color=line.get_color(), lw=1)
        ax[1].semilogy(t[1:], np.maximum(np.abs(w - exact)[1:] / meta["v0"], 1e-17), "o-", ms=3,
                       color=line.get_color(), label=label)
    ax[0].set_xlabel("t"); ax[0].set_ylabel("(v_dust - v_gas)/v0"); ax[0].set_yscale("log")
    ax[0].set_ylim(1e-12, 2); ax[0].legend(fontsize=8); ax[0].grid(alpha=0.3)
    ax[0].set_title("relative velocity, points: fishstep, lines: exact")
    ax[1].set_xlabel("t"); ax[1].set_ylabel("|error|/v0"); ax[1].grid(alpha=0.3)
    ax[1].set_title("deviation from the exact solution")
    fig.suptitle(f"DUSTYBOX, dust to gas ratio {runs[0][0]['eps']:g}")
    fig.tight_layout()
    fig.savefig(png, dpi=100)
    print(f"plot written to {png}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
