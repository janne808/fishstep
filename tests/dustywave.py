#!/usr/bin/env python3
"""Check the fishstep gas-dust coupling against the linear DUSTYWAVE solution.

Usage: python3 tests/dustywave.py dustywave.dat [dustywave2.dat ...] [output.png]

Each data file is written by the DUSTY_BOX setup with DB_WAVE > 0: a sound
wave of wave number k in gas and dust with dust to gas ratio eps (Laibe and
Price 2011). The file holds per step the Fourier coefficients at k of the gas
and dust density and momentum, normalized by the species mass,

    d = sum_p m exp(-i k x_p)/M,   p = sum_p m v_p exp(-i k x_p)/M.

To first order in the amplitude they obey

    dd_g/dt = -i k p_g,   dp_g/dt = -i k cs^2 d_g - eps (p_g - p_d)/t_s
    dd_d/dt = -i k p_d,   dp_d/dt = (p_g - p_d)/t_s

which is solved exactly from the measured initial coefficients. Without
coupling it is a sound wave at cs, for t_s -> 0 the dusty gas carries sound
at cs/sqrt(1 + eps), in between the drag damps the wave, most strongly for
k cs t_s ~ 1. Reported is the deviation of the gas density coefficient from
the linear solution relative to its initial value. The run without coupling
(the largest t_s) gives the error of the gas dynamics alone.
"""

import sys

import numpy as np


def read(path):
    meta = {}
    for line in open(path):
        if line.startswith("#"):
            tok = line[1:].split()
            if tok and tok[0] == "eps":
                meta = {tok[i]: float(tok[i + 1]) for i in range(0, len(tok) - 1, 2)}
    return meta, np.loadtxt(path)


def linear(meta, t, x0):
    eps, ts, cs, k = meta["eps"], meta["t_s"], meta["cs"], meta["k"]
    A = np.array([[0, -1j * k, 0, 0],
                  [-1j * k * cs ** 2, -eps / ts, 0, eps / ts],
                  [0, 0, 0, -1j * k],
                  [0, 1 / ts, 0, -1 / ts]], dtype=complex)
    lam, vec = np.linalg.eig(A)
    c = np.linalg.solve(vec, x0)
    return np.real_if_close((vec[None, :, :] * (c * np.exp(lam[None, :] * t[:, None]))[:, None, :]).sum(-1)), lam


def main():
    args = sys.argv[1:] or ["dustywave.dat"]
    png = args.pop() if args[-1].endswith(".png") else args[0].rsplit(".", 1)[0] + ".png"
    runs = []
    print("     t_s   k cs t_s    dt/tau    wave speed/cs  damping time   max |d_gas - linear|/|d_gas(0)|   at t_end")
    for path in args:
        meta, d = read(path)
        t = d[:, 0]
        x = d[:, 1:9:2] + 1j * d[:, 2:9:2]
        lin, lam = linear(meta, t, x[0])
        err = np.abs(x[:, 0] - lin[:, 0]) / np.abs(x[0, 0])
        # the gas-like mode: largest oscillation frequency
        mode = lam[np.argmax(np.abs(lam.imag))]
        speed = abs(mode.imag) / meta["k"] / meta["cs"]
        damp = -1.0 / mode.real if mode.real < 0 else np.inf
        tau = meta["t_s"] / (1 + meta["eps"])
        print(f"{meta['t_s']:9.3g}  {meta['k'] * meta['cs'] * meta['t_s']:8.3g}  {meta['dt'] / tau:9.3g}   "
              f"{speed:10.4f}    {damp:10.4g}       {err.max():14.4e}               {err[-1]:.4e}")
        runs.append((meta, t, x, lin))

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available, no plot written")
        return

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for meta, t, x, lin in runs:
        label = f"t_s = {meta['t_s']:g}"
        line, = ax[0].plot(t, x[:, 0].real / abs(x[0, 0]), "o", ms=2, label=label)
        ax[0].plot(t, lin[:, 0].real / abs(x[0, 0]), "-", color=line.get_color(), lw=1)
        ax[1].plot(t, x[:, 2].real / abs(x[0, 2]), "o", ms=2, color=line.get_color(), label=label)
        ax[1].plot(t, lin[:, 2].real / abs(x[0, 2]), "-", color=line.get_color(), lw=1)
    ax[0].set_title("gas density coefficient, points: fishstep, lines: linear")
    ax[1].set_title("dust density coefficient")
    for a in ax:
        a.set_xlabel("t"); a.grid(alpha=0.3); a.legend(fontsize=8)
    fig.suptitle(f"DUSTYWAVE, dust to gas ratio {runs[0][0]['eps']:g}")
    fig.tight_layout()
    fig.savefig(png, dpi=100)
    print(f"plot written to {png}")


if __name__ == "__main__":
    main()
