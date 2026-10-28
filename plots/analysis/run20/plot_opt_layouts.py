"""Run 20 (Malata): every validated layout (diagnostics/swgotrig/validate_count_layouts.txt: baselines, Stage 4 runs with
and without the cable cost) on the detector region. Each panel: the 266 faces (East-North, coloured by Up), the
detectors, the 3-D minimum spanning tree linking them (modules/optimize/cable.py) and its length, and the validated
number of triggered showers within 1 / 2.5 / 5 deg for the drawn-count recon (valcnt_count score files):
no background, all read out and with the kNN gate. Writes <run 20>/diagnostics/opt_layouts.png.

    source slurm/run20_guard.sh; python plots/analysis/run20/plot_opt_layouts.py
"""
import glob
import json
import os
import re
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


def within(rd, w):
    """{layout: number of triggered showers within 1 / 2.5 / 5 deg}, drawn-count recon, summed over the held-out chunks
    and averaged over the count seeds (all_X = triggered and within X deg, as a share of all showers)."""
    files = glob.glob(os.path.join(D, "valcnt_count_s*_*.json"))
    n_seeds = len({re.search(r"_s(\d+)_\d+\.json$", f).group(1) for f in files})
    num = {}
    for f in files:
        for k, by in json.load(open(f))["scores"].items():
            r = by[rd][w]
            num[k] = num.get(k, 0.0) + r["n_showers"] * np.array([r["all_1"], r["all_2.5"], r["all_5"]])
    return {k: v / n_seeds for k, v in num.items()}


N_HELD = sum(json.load(open(f))["scores"]["grid_n99"]["none"]["0"]["n_showers"]
             for f in glob.glob(os.path.join(D, "valcnt_count_s1_*.json")))
clean, knn = within("none", "0"), within("knn", "0")
fmt = lambda v: "–" if v is None else "  ".join(f"{x:,.0f}" for x in v)                   # noqa: E731
m = load_tr_mountain(C.GEOMETRY_PATH_RESOLVED, C.GEOMETRY_GROUP, C.DET_KEY, east_entry=C.EAST_ENTRY,
                     layer_east_dx=C.LAYER_EAST_DX, n_planes=C.N_PLANES)
surface = SurfaceUpMap.from_mountain(m, grid_h=256, grid_w=256)
verts, faces, det_idx, loc = _read_detector_region(C.GEOMETRY_PATH_RESOLVED, C.GEOMETRY_GROUP, C.DET_KEY)
tri = np.stack([_ecef_to_enu(verts[:, faces[k, det_idx]], loc[0], loc[1]).T for k in range(3)], 1)   # (F, 3, 3) ENU

panels = layouts(D)
ncol = 5
nrow = -(-len(panels) // ncol)
fig, axs = plt.subplots(nrow, ncol, figsize=(3.9 * ncol, 4.9 * nrow), constrained_layout=True)
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
    ax.set_title(f"{name}, cable {L:.2f} km\nall read out {fmt(clean.get(key))}\nkNN gate {fmt(knn.get(key))}", fontsize=8.5)
    ax.tick_params(labelsize=7)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    print(f"{name}: cable {L:.2f} km, showers within 1/2.5/5 deg: {fmt(clean.get(key))} (all read out), "
          f"{fmt(knn.get(key))} (kNN gate)")
for ax in axs.flat[len(panels):]:
    ax.axis("off")
for ax in axs[:, 0]:
    ax.set_ylabel("North [m]")
for ax in axs[-1]:
    ax.set_xlabel("East [m]")
fig.suptitle(f"Run 20, Malata: validated layouts (blue) and their 3-D cable tree (orange). Triggered showers within "
             f"1° / 2.5° / 5° (of {N_HELD:,} held-out, mean of the count seeds), drawn-count recon, no background: "
             "all read out and kNN gate", fontsize=11)
fig.savefig(OUT, dpi=170, facecolor="white")
print("[write]", OUT)
