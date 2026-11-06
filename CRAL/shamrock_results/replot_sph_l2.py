"""Re-plot the L2-error-vs-time curve of Shamrock's examples/sph/run_sod.py.

The example draws this plot onto the figure left over from its GIF animation; with a
non-interactive backend (MPLBACKEND=Agg) that figure is never cleared, so the saved
l2_error.png comes out garbled. The underlying data (l2_error_*.npy, written by the
example's AnalysisHelper) are fine -- this script only re-draws them.

usage: python replot_sph_l2.py <sim_folder, e.g. .../_to_trash/sod_128> <out.png>
"""
import glob
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

folder, out = sys.argv[1], sys.argv[2]
files = sorted(glob.glob(f"{folder}/l2_error_*.npy"))
data = [np.load(f, allow_pickle=True).item() for f in files]
t = np.array([d["time"] for d in data])
l2_rho = np.array([d["l2_rho"] for d in data])
l2_v = np.array([d["l2_v"] for d in data])
l2_P = np.array([d["l2_P"] for d in data])

fig, ax = plt.subplots(figsize=(6, 4), dpi=150)
ax.plot(t, l2_rho, "o-", label=r"$L_2(\rho)$")
ax.plot(t, l2_v[:, 0], "o-", label=r"$L_2(v_x)$")
ax.plot(t, l2_v[:, 1], "o-", label=r"$L_2(v_y)$")
ax.plot(t, l2_v[:, 2], "o-", label=r"$L_2(v_z)$")
ax.plot(t, l2_P, "o-", label=r"$L_2(P)$")
ax.set_yscale("log")
ax.set_xlabel("t")
ax.set_ylabel("L2 distance to analytic Sod (Shamrock compute_L2_dist)")
ax.set_title(f"SPH Sod tube, {len(files)} snapshots ({folder.rstrip('/').split('/')[-1]})")
ax.legend(fontsize=8)
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(out)
for ti, r, v, p in zip(t, l2_rho, l2_v[:, 0], l2_P):
    print(f"t={ti:.4f}  L2_rho={r:.3e}  L2_vx={v:.3e}  L2_P={p:.3e}")
print("saved", out)
