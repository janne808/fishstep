#!/usr/bin/env python3
"""Check that the fishstep protostellar gas disk stays in equilibrium rotation.

Usage: python3 tests/stardisk.py stardisk.dat [output.png]

The data file is written by the default fishstep setup (STAR_DISK): one block
per output time with a header line holding the time, the star position and
velocity, the center of mass, the total angular momentum Lz and the disk
parameters, followed by radial profiles around the star in bins of one cell:
r, surface density, mean azimuthal velocity, mean radial velocity and rms
radial velocity. The t = 0 block is written right after the rotation has been
set from the initial forces.
"""

import sys

import numpy as np


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


def sigma_shape(r, r_in, r_out, w):
    return 0.25 * (1 + np.tanh((r - r_in) / w)) * (1 - np.tanh((r - r_out) / w)) / r


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "stardisk.dat"
    png = sys.argv[2] if len(sys.argv) > 2 else path.rsplit(".", 1)[0] + ".png"

    blocks = read_blocks(path)
    m0, p0 = blocks[0]
    r_in, r_out, v_k, cs = m0["R_in"], m0["R_out"], m0["v_K"], m0["cs"]
    t_orb = 2 * np.pi * r_out / v_k

    # compare inside the disk body, away from the tapered edges
    def body(r):
        return (r > r_in + 6) & (r < r_out - 6)

    # exact shape of the initial surface density, normalized on the body
    r0 = p0[:, 0]
    shape = sigma_shape(r0, r_in, r_out, 3.0)
    b0 = body(r0)
    norm = np.sum(p0[b0, 1]) / np.sum(shape[b0])

    print(f"protostellar disk: disk/star mass {m0['disk_to_star']:g}, R {r_in:g} - {r_out:g} cells, "
          f"v_K {v_k:.3f}, cs {cs:.3f} cells/time, outer orbit {t_orb:.2f}")
    print(f"t = 0 surface density vs exact profile: L1 {np.mean(np.abs(p0[b0, 1] - norm * shape[b0])) / np.mean(norm * shape[b0]):.4f}")
    print("   t   orbits(R_out)  sigma L1 vs t=0  v_phi change  <v_r>/cs   rms v_r/cs   star offset   Lz/Lz0 - 1")
    lz0 = m0["Lz"]
    rows = []
    for meta, prof in blocks:
        r = prof[:, 0]
        common = np.intersect1d(r[body(r)], r0[b0])
        i = np.searchsorted(r, common)
        j = np.searchsorted(r0, common)
        sig_l1 = np.mean(np.abs(prof[i, 1] - p0[j, 1])) / np.mean(p0[j, 1])
        dvphi = np.mean(np.abs(prof[i, 2] - p0[j, 2])) / np.mean(p0[j, 2])
        w = prof[i, 1] * common
        vr_mean = np.sum(w * prof[i, 3]) / np.sum(w) / cs
        vr_rms = np.sqrt(np.sum(w * prof[i, 4] ** 2) / np.sum(w)) / cs
        off = np.hypot(meta["star_x"] - meta["cm_x"], meta["star_y"] - meta["cm_y"])
        dlz = meta["Lz"] / lz0 - 1
        rows.append((meta["t"], sig_l1, dvphi, vr_rms))
        print(f"{meta['t']:5.1f}  {meta['t'] / t_orb:8.2f}      {100 * sig_l1:7.2f}%       {100 * dvphi:7.3f}%   "
              f"{vr_mean:+8.4f}   {vr_rms:8.4f}      {off:8.4f}     {dlz:+.2e}")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available, no plot written")
        return

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 4.5))
    rr = np.linspace(1, r_out + 18, 400)
    ax1.plot(rr, norm * sigma_shape(rr, r_in, r_out, 3.0), color="0.2", lw=1.5, label="initial profile")
    picks = sorted(set([0, len(blocks) // 2, len(blocks) - 1]))
    colors = ("tab:blue", "tab:orange", "tab:red")
    for c, k in zip(colors, picks):
        meta, prof = blocks[k]
        ax1.plot(prof[:, 0], prof[:, 1], ".", ms=3, color=c, label=f"t = {meta['t']:g}")
        ax2.plot(prof[:, 0], prof[:, 2], ".", ms=3, color=c, label=f"t = {meta['t']:g}")
    ax2.axhline(v_k, color="0.2", lw=1, ls="--", label="v_K of the star")
    ax1.set_yscale("log")
    ax1.set_xlabel("r [cells]"); ax1.set_ylabel("surface density"); ax1.legend(); ax1.grid(alpha=0.3)
    ax1.set_title("surface density")
    ax2.set_xlabel("r [cells]"); ax2.set_ylabel("v_phi [cells/time]"); ax2.legend(); ax2.grid(alpha=0.3)
    ax2.set_title("rotation curve")
    rows = np.array(rows)
    ax3.plot(rows[:, 0] / t_orb, 100 * rows[:, 1], "o-", ms=3, label="sigma L1 change [%]")
    ax3.plot(rows[:, 0] / t_orb, 100 * rows[:, 2], "o-", ms=3, label="v_phi change [%]")
    ax3.plot(rows[:, 0] / t_orb, rows[:, 3], "o-", ms=3, label="rms v_r / cs")
    ax3.set_xlabel("outer orbits"); ax3.legend(); ax3.grid(alpha=0.3)
    ax3.set_title("departure from the initial equilibrium")
    fig.suptitle("Protostellar gas disk in equilibrium rotation")
    fig.tight_layout()
    fig.savefig(png, dpi=110)
    print(f"plot written to {png}")


if __name__ == "__main__":
    main()
