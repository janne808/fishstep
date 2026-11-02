#!/usr/bin/env python3
"""Follow the fishstep protostellar gas disk as it cools and fragments.

Usage: python3 tests/stardisk.py stardisk.dat [output.png]

The data file is written by the default fishstep setup (STAR_DISK): one block
per output time with a header line holding the time, the star position and
velocity, the center of mass, the total angular momentum Lz, the disk and
cooling parameters, the gas thermal energy E_th, the energy radiated so far
E_rad and the mass accreted by the star M_acc, followed by radial profiles
around the star in bins of one cell: r, surface density, mean azimuthal
velocity, mean radial velocity, rms radial velocity and mass weighted sound
speed. The t = 0 block is written right after the rotation has been set from
the initial forces. stardisk_rho.bin next to it holds one gas surface density
map (NUM_CELLS^2 float32) per block.

With dust (SD_DUST 1) the profiles have three more columns, dust surface
density, mean azimuthal and radial dust velocity, the header holds the dust to
gas ratio, the Stokes number and the accreted dust mass, and
stardisk_dust.bin holds the dust surface density maps. The script then also
reports the fraction of the dust in the gas clumps and their dust to gas
ratio relative to the initial one.

Per output the script reports the Toomre parameter Q = cs kappa/(pi G sigma)
of the disk body, the thermal and radiated energy, the accreted mass and a
census of gas clumps in the maps: local density maxima more than CONTRAST
times the azimuthal mean density at their radius holding at least
CLUMP_MIN_MASS of the disk mass within CLUMP_RADIUS cells.
"""

import os
import sys

import numpy as np

CONTRAST = 30.0
CLUMP_RADIUS = 3.0
CLUMP_MIN_MASS = 1.0e-4
R_ACC = 8.0


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


def toomre_q(prof, G):
    """Q(r) from the binned profile, kappa^2 = (2 v/r^2) d(r v)/dr."""
    r, sigma, vphi, cs = prof[:, 0], prof[:, 1], prof[:, 2], prof[:, 5]
    v = np.convolve(vphi, np.ones(3) / 3, mode="same")
    kap2 = 2.0 * v / r ** 2 * np.gradient(r * v, r)
    with np.errstate(invalid="ignore", divide="ignore"):
        q = cs * np.sqrt(np.where(kap2 > 0.0, kap2, np.nan)) / (np.pi * G * sigma)
    return q


def find_clumps(rho, cx, cy, m_disk, dust=None):
    n = rho.shape[0]
    yy, xx = np.mgrid[0:n, 0:n].astype(float)
    r = np.hypot(xx - cx, yy - cy)
    rb = r.astype(int)
    mean = np.bincount(rb.ravel(), rho.ravel()) / np.maximum(np.bincount(rb.ravel()), 1)
    peak = rho > CONTRAST * mean[rb]
    for dy in range(-2, 3):
        for dx in range(-2, 3):
            if dx or dy:
                peak &= rho >= np.roll(np.roll(rho, dy, 0), dx, 1)
    peak &= r > R_ACC + 2
    clumps = []
    for y0, x0 in zip(*np.nonzero(peak)):
        sl = (slice(max(y0 - 4, 0), y0 + 5), slice(max(x0 - 4, 0), x0 + 5))
        near = np.hypot(xx[sl] - x0, yy[sl] - y0) <= CLUMP_RADIUS
        m = rho[sl][near].sum()
        if m >= CLUMP_MIN_MASS * m_disk:
            md = dust[sl][near].sum() if dust is not None else 0.0
            clumps.append((x0, y0, r[y0, x0], m, md))
    return clumps


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "stardisk.dat"
    png = sys.argv[2] if len(sys.argv) > 2 else path.rsplit(".", 1)[0] + ".png"
    map_path = path.rsplit(".", 1)[0] + "_rho.bin"
    dust_path = path.rsplit(".", 1)[0] + "_dust.bin"

    blocks = read_blocks(path)
    m0, p0 = blocks[0]
    r_in, r_out, v_k, G = m0["R_in"], m0["R_out"], m0["v_K"], m0["G"]
    edge_out = m0.get("edge_out", 6.0)
    t_orb = 2 * np.pi * r_out / v_k
    m_disk = m0["disk_to_star"] / (1 + m0["disk_to_star"])
    e_th0 = m0["E_th"]

    maps = None
    if os.path.exists(map_path):
        raw = np.fromfile(map_path, dtype=np.float32)
        n = int(round(np.sqrt(raw.size / len(blocks))))
        if n * n * len(blocks) == raw.size:
            maps = raw.reshape(len(blocks), n, n)
        else:
            print(f"{map_path}: {raw.size} values do not match {len(blocks)} frames, no clump census")
    dmaps = None
    eps0 = m0.get("dust_to_gas", 0.0)
    if maps is not None and eps0 > 0 and os.path.exists(dust_path):
        raw = np.fromfile(dust_path, dtype=np.float32)
        if raw.size == maps.size:
            dmaps = raw.reshape(maps.shape)
        else:
            print(f"{dust_path}: {raw.size} values do not match the gas maps, no dust census")

    # the disk body inside the tapered edges
    def body(r):
        return (r > r_in + 6) & (r < r_out - max(edge_out, 6))

    # Toomre Q away from the pressure supported edges
    def q_range(r):
        return (r > r_in + 8) & (r < r_out - max(edge_out, 12))

    print(f"protostellar disk: disk/star mass {m0['disk_to_star']:g}, R {r_in:g} - {r_out:g} cells, "
          f"outer orbit {t_orb:.2f}, gamma {m0['gamma']:.3f}, Q_init {m0['Q_init']:g}, "
          + (f"beta cooling {m0['beta']:g} to the irradiation floor, Q_floor {m0['Q_floor']:g}" if m0["cooling"] else "no cooling"))
    if eps0 > 0:
        print(f"dust: dust/gas mass {eps0:g}, St {m0['St']:g} at the initial gas surface density in the middle of the disk")
    print("   t  orbits(R_out)  sigma L1   <v_r>/cs   Q min   Q median   E_th/E_th0  E_rad/E_th0  "
          "M_acc/M_disk  Lz/Lz0-1   clumps  clump mass  largest"
          + ("   dust in clumps  clump dust/gas / eps0  dust accreted" if dmaps is not None else ""))
    lz0 = m0["Lz"]
    r0 = p0[:, 0]
    b0 = body(r0)
    series = []
    clump_sets = []
    for k, (meta, prof) in enumerate(blocks):
        r = prof[:, 0]
        common = np.intersect1d(r[body(r)], r0[b0])
        i = np.searchsorted(r, common)
        j = np.searchsorted(r0, common)
        sig_l1 = np.mean(np.abs(prof[i, 1] - p0[j, 1])) / np.mean(p0[j, 1])
        w = prof[i, 1] * common
        vr_mean = np.sum(w * prof[i, 3]) / np.sum(w * prof[i, 5])
        q = toomre_q(prof, G)[q_range(r)]
        clumps = find_clumps(maps[k], meta["star_x"], meta["star_y"], m_disk,
                             dmaps[k] if dmaps is not None else None) if maps is not None else []
        clump_sets.append(clumps)
        cm = sum(c[3] for c in clumps) / m_disk
        cmax = max((c[3] for c in clumps), default=0.0) / m_disk
        series.append((meta["t"], np.nanmin(q), np.nanmedian(q), meta["E_th"] / e_th0, meta["E_rad"] / e_th0,
                       len(clumps), cm))
        print(f"{meta['t']:5.1f}  {meta['t'] / t_orb:8.2f}     {100 * sig_l1:6.2f}%   {vr_mean:+8.4f}  "
              f"{np.nanmin(q):6.2f}  {np.nanmedian(q):7.2f}     {meta['E_th'] / e_th0:8.4f}    "
              f"{meta['E_rad'] / e_th0:8.4f}      {meta['M_acc'] / m_disk:.2e}    {meta['Lz'] / lz0 - 1:+.2e}   "
              f"{len(clumps):4d}    {cm:8.4f}   {cmax:.4f}"
              + (f"      {sum(c[4] for c in clumps) / (eps0 * m_disk):8.4f}       "
                 f"{(sum(c[4] for c in clumps) / max(sum(c[3] for c in clumps), 1e-300)) / eps0:10.3f}           "
                 f"{meta['M_acc_dust'] / (eps0 * m_disk):.3e}" if dmaps is not None else ""))
    if clump_sets and clump_sets[-1]:
        last = sorted(clump_sets[-1], key=lambda c: -c[3])
        print("clumps at the last output, radius and mass in disk masses:",
              ", ".join(f"r {c[2]:.0f} m {c[3] / m_disk:.4f}" for c in last[:16]))

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available, no plot written")
        return

    picks = sorted(set(np.linspace(0, len(blocks) - 1, 4).astype(int)))
    fig = plt.figure(figsize=(16, 8.5))
    win = r_out + 2 * edge_out + 10
    for col, k in enumerate(picks):
        ax = fig.add_subplot(2, 4, col + 1)
        meta = blocks[k][0]
        if maps is not None:
            cx, cy = meta["star_x"], meta["star_y"]
            n = maps.shape[1]
            lo_x, lo_y = int(cx - win), int(cy - win)
            img = maps[k][max(lo_y, 0):int(cy + win), max(lo_x, 0):int(cx + win)]
            floor = np.percentile(maps[0][maps[0] > 0], 20)
            ax.imshow(np.log10(np.maximum(img, floor)), origin="lower", cmap="magma",
                      extent=(lo_x - cx, lo_x - cx + img.shape[1], lo_y - cy, lo_y - cy + img.shape[0]))
            for (x0, y0, rr, m, *_) in clump_sets[k]:
                ax.plot(x0 - cx, y0 - cy, "o", mfc="none", mec="cyan", ms=8, mew=1.2)
        ax.set_title(f"t = {meta['t']:g}, clumps {len(clump_sets[k])}")
        ax.set_xlabel("x - x_star [cells]")
        if col == 0:
            ax.set_ylabel("log gas surface density")

    colors = plt.cm.viridis(np.linspace(0, 0.9, len(picks)))
    ax = fig.add_subplot(2, 4, 5)
    for c, k in zip(colors, picks):
        meta, prof = blocks[k]
        ax.plot(prof[:, 0], prof[:, 1], color=c, lw=1, label=f"t = {meta['t']:g}")
        if eps0 > 0 and prof.shape[1] > 6:
            ax.plot(prof[:, 0], prof[:, 6] / eps0, color=c, lw=1, ls="--")
    ax.set_yscale("log"); ax.set_xlabel("r [cells]")
    ax.set_title("azimuthal mean surface density" + (", dashed: dust/eps0" if eps0 > 0 else ""))
    ax.set_xlim(0, r_out + 20); ax.legend(fontsize=8); ax.grid(alpha=0.3)

    ax = fig.add_subplot(2, 4, 6)
    for c, k in zip(colors, picks):
        meta, prof = blocks[k]
        sel = q_range(prof[:, 0])
        ax.plot(prof[sel, 0], toomre_q(prof, G)[sel], color=c, lw=1, label=f"t = {meta['t']:g}")
    ax.axhline(1.0, color="0.3", ls="--", lw=1)
    ax.set_yscale("log"); ax.set_xlabel("r [cells]"); ax.set_title("Toomre Q"); ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    s = np.array(series)
    ax = fig.add_subplot(2, 4, 7)
    ax.plot(s[:, 0] / t_orb, s[:, 1], "o-", ms=3, label="Q min")
    ax.plot(s[:, 0] / t_orb, s[:, 2], "o-", ms=3, label="Q median")
    ax.axhline(1.0, color="0.3", ls="--", lw=1)
    ax2 = ax.twinx()
    ax2.plot(s[:, 0] / t_orb, s[:, 5], "s-", ms=3, color="tab:red", label="clumps")
    ax2.set_ylabel("number of clumps", color="tab:red")
    ax.set_xlabel("outer orbits"); ax.set_yscale("log"); ax.legend(fontsize=8, loc="upper left")
    ax.set_title("disk stability and fragments"); ax.grid(alpha=0.3)

    ax = fig.add_subplot(2, 4, 8)
    ax.plot(s[:, 0] / t_orb, s[:, 3], "o-", ms=3, label="E_th / E_th0")
    ax.plot(s[:, 0] / t_orb, s[:, 4], "o-", ms=3, label="E_rad / E_th0")
    ax.plot(s[:, 0] / t_orb, 100 * s[:, 6], "o-", ms=3, label="clump mass [% of disk]")
    ax.set_xlabel("outer orbits"); ax.legend(fontsize=8); ax.grid(alpha=0.3)
    ax.set_title("thermal budget and clump mass")

    fig.suptitle("Protostellar gas disk: " + (f"beta = {m0['beta']:g} cooling" if m0["cooling"] else "no cooling"))
    fig.tight_layout()
    fig.savefig(png, dpi=100)
    print(f"plot written to {png}")


if __name__ == "__main__":
    main()
