"""Run 20 (Malata): every validated layout (diagnostics/swgotrig/validate_count_layouts.txt: baselines, Stage 4 runs with
and without the cable cost) on the detector region. Each panel: the 266 faces (East-North, coloured by Up), the
detectors, the 3-D minimum spanning tree linking them (modules/optimize/cable.py) and its length, and the validated
share of triggered showers within 5 deg for the drawn-count recon (valcnt_count score files, all seeds and chunks):
no background, and 10 us background with the kNN gate. Writes <run 20>/diagnostics/opt_layouts.png.

    source slurm/run20_guard.sh; python plots/analysis/run20/plot_opt_layouts.py
"""
import glob
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import torch  # noqa: E402
from matplotlib.collections import LineCollection, PolyCollection  # noqa: E402

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, "plots", "analysis", "run19", "swgotrig"))
import modules.constants as C  # noqa: E402
from modules.geometry import load_tr_mountain, SurfaceUpMap  # noqa: E402
from modules.geometry.mountain import _read_detector_region, _ecef_to_enu  # noqa: E402
from modules.optimize.cable import cable_length, mst_edges, points_3d  # noqa: E402
from layout_names import layouts  # noqa: E402

D = os.path.join(C.RUN_LOCATION, "diagnostics", "swgotrig")
OUT = os.path.join(C.RUN_LOCATION, "diagnostics", "opt_layouts.png")
PATHS = dict(l.strip().split("=", 1) for l in open(os.path.join(D, "validate_count_layouts.txt")) if "=" in l)


def within5(rd, w):
    """{layout: share of triggered showers within 5 deg}, drawn-count recon, all seeds and chunks pooled."""
    num, den = {}, {}
    for f in glob.glob(os.path.join(D, "valcnt_count_s*_*.json")):
        for k, by in json.load(open(f))["scores"].items():
            r = by[rd][w]
            num[k] = num.get(k, 0.0) + r["n_showers"] * r["all_5"]
            den[k] = den.get(k, 0.0) + r["n_showers"] * r["triggered"]
    return {k: num[k] / den[k] for k in num if den[k] > 0}


clean, noisy = within5("none", "0"), within5("knn", "10000")
m = load_tr_mountain(C.GEOMETRY_PATH_RESOLVED, C.GEOMETRY_GROUP, C.DET_KEY, east_entry=C.EAST_ENTRY,
                     layer_east_dx=C.LAYER_EAST_DX, n_planes=C.N_PLANES)
surface = SurfaceUpMap.from_mountain(m, grid_h=256, grid_w=256)
verts, faces, det_idx, loc = _read_detector_region(C.GEOMETRY_PATH_RESOLVED, C.GEOMETRY_GROUP, C.DET_KEY)
tri = np.stack([_ecef_to_enu(verts[:, faces[k, det_idx]], loc[0], loc[1]).T for k in range(3)], 1)   # (F, 3, 3) ENU

panels = layouts(D)
ncol = 5
nrow = -(-len(panels) // ncol)
fig, axs = plt.subplots(nrow, ncol, figsize=(3.9 * ncol, 4.5 * nrow), constrained_layout=True)
axs = np.atleast_2d(axs)
xl = (tri[..., 0].min() - 60, tri[..., 0].max() + 60)
yl = (tri[..., 1].min() - 60, tri[..., 1].max() + 60)
for ax, (key, name) in zip(axs.flat, panels):
    lay = torch.load(PATHS[key], map_location="cpu")
    x, y = lay["x"].float(), lay["y"].float()
    with torch.no_grad():
        L = float(cable_length(x, y, surface)[0]) / 1e3
        i, j = mst_edges(points_3d(x, y, surface))
    ax.add_collection(PolyCollection(tri[:, :, :2], array=tri[:, :, 2].mean(1), cmap="Greys", edgecolors="white",
                                     linewidths=0.3))
    xe, ye = x.numpy(), y.numpy()
    ax.add_collection(LineCollection(np.stack([np.c_[xe[i], ye[i]], np.c_[xe[j], ye[j]]], 1),
                                     colors="#e07b39", linewidths=0.9))
    ax.scatter(xe, ye, s=10, c="#2a6fb0", lw=0, zorder=3)
    ax.set_xlim(*xl); ax.set_ylim(*yl); ax.set_aspect("equal")
    ax.set_title(f"{name}\ncable {L:.2f} km   5°: {100 * clean.get(key, np.nan):.1f}% / "
                 f"{100 * noisy.get(key, np.nan):.1f}%", fontsize=9)
    ax.tick_params(labelsize=7)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    print(f"{name}: cable {L:.2f} km, 5 deg of triggered {100 * clean.get(key, np.nan):.1f}% (no bg) / "
          f"{100 * noisy.get(key, np.nan):.1f}% (10 us, kNN)")
for ax in axs.flat[len(panels):]:
    ax.axis("off")
for ax in axs[:, 0]:
    ax.set_ylabel("North [m]")
for ax in axs[-1]:
    ax.set_xlabel("East [m]")
fig.suptitle("Run 20, Malata: validated layouts (blue) and their 3-D cable tree (orange). 5° = triggered showers within "
             "5°, drawn-count recon: no background / 10 µs background with the kNN gate", fontsize=11)
fig.savefig(OUT, dpi=170, facecolor="white")
print("[write]", OUT)
