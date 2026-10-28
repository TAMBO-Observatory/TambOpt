"""Run 20 (Malata): the Stage 4 layouts (layout_best.pt of every swgo_trigger run, with and without the cable cost) and
the grid_n99 baseline on the detector region. Each panel: the 266 faces (East-North, coloured by Up), the detectors,
the 3-D minimum spanning tree linking them (modules/optimize/cable.py) and its length; the title gives the best
scoring-set U of the run (with the cable runs, U minus the cable penalty, as optimised).
Writes <run 20>/diagnostics/opt_layouts.png.

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
import modules.constants as C  # noqa: E402
from modules.geometry import load_tr_mountain, SurfaceUpMap  # noqa: E402
from modules.geometry.mountain import _read_detector_region, _ecef_to_enu  # noqa: E402
from modules.optimize.cable import cable_length, mst_edges, points_3d  # noqa: E402

OUT = os.path.join(C.RUN_LOCATION, "diagnostics", "opt_layouts.png")
PRE = os.path.join(C.RUN_LOCATION, "test_v6_run_04_optimize_lbfgs_kernel_swgotrig_l0p3_k3_")
RUNS = [("a2_lr10_grid_b2048", "Poisson 2 m², grid start"), ("a2_lr10_center_b2048", "Poisson 2 m², centre start"),
        ("a20_lr10_grid_b2048", "Poisson 20 m², grid start"), ("a2_step_lr10_grid_b2048", "Step 2 m², grid start"),
        ("a2_step_lr10_center_b2048", "Step 2 m², centre start"),
        ("a2_lr10_grid_b2048_cab0p05", "Cable λ 0.05, grid start"), ("a2_lr10_grid_b2048_cab0p1", "Cable λ 0.1, grid start"),
        ("a2_lr10_grid_b2048_cab0p2", "Cable λ 0.2, grid start"), ("a2_lr10_center_b2048_cab0p1", "Cable λ 0.1, centre start")]

m = load_tr_mountain(C.GEOMETRY_PATH_RESOLVED, C.GEOMETRY_GROUP, C.DET_KEY, east_entry=C.EAST_ENTRY,
                     layer_east_dx=C.LAYER_EAST_DX, n_planes=C.N_PLANES)
surface = SurfaceUpMap.from_mountain(m, grid_h=256, grid_w=256)
verts, faces, det_idx, loc = _read_detector_region(C.GEOMETRY_PATH_RESOLVED, C.GEOMETRY_GROUP, C.DET_KEY)
tri = np.stack([_ecef_to_enu(verts[:, faces[k, det_idx]], loc[0], loc[1]).T for k in range(3)], 1)   # (F, 3, 3) ENU

grid = glob.glob(os.path.join(C.RUN_LOCATION, "diagnostics", "triplets", "layouts", "grid_n99.pt"))
panels = [(grid[0], "Grid (baseline)", None)] if grid else []
for sfx, name in RUNS:
    d = PRE + sfx
    if os.path.exists(os.path.join(d, "layout_best.pt")):
        u = json.load(open(os.path.join(d, "optimize_log.json")))["best_u"]
        panels.append((os.path.join(d, "layout_best.pt"), name, u))

ncol = 5
nrow = -(-len(panels) // ncol)
fig, axs = plt.subplots(nrow, ncol, figsize=(3.9 * ncol, 4.3 * nrow), constrained_layout=True)
axs = np.atleast_2d(axs)
xl = (tri[..., 0].min() - 60, tri[..., 0].max() + 60)
yl = (tri[..., 1].min() - 60, tri[..., 1].max() + 60)
for ax, (path, name, u) in zip(axs.flat, panels):
    lay = torch.load(path, map_location="cpu")
    x, y = lay["x"].float(), lay["y"].float()
    with torch.no_grad():
        L = float(cable_length(x, y, surface)[0]) / 1e3
        i, j = mst_edges(points_3d(x, y, surface))
    pc = PolyCollection(tri[:, :, :2], array=tri[:, :, 2].mean(1), cmap="Greys", edgecolors="white", linewidths=0.3)
    ax.add_collection(pc)
    xe, ye = x.numpy(), y.numpy()
    ax.add_collection(LineCollection(np.stack([np.c_[xe[i], ye[i]], np.c_[xe[j], ye[j]]], 1),
                                     colors="#e07b39", linewidths=0.9))
    ax.scatter(xe, ye, s=10, c="#2a6fb0", lw=0, zorder=3)
    ax.set_xlim(*xl); ax.set_ylim(*yl); ax.set_aspect("equal")
    ax.set_title(f"{name}\ncable {L:.2f} km" + (f"   U {u:.2f}" if u is not None else ""), fontsize=9.5)
    ax.tick_params(labelsize=7)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    print(f"{name}: {len(xe)} detectors, cable {L:.2f} km" + (f", best U {u:.3f}" if u is not None else ""))
for ax in axs.flat[len(panels):]:
    ax.axis("off")
for ax in axs[:, 0]:
    ax.set_ylabel("North [m]")
for ax in axs[-1]:
    ax.set_xlabel("East [m]")
fig.suptitle("Run 20, Malata: optimised layouts (blue) and their 3-D cable tree (orange). "
             "U = best scoring-set utility (cable runs: minus the cable penalty)", fontsize=11)
fig.savefig(OUT, dpi=170, facecolor="white")
print("[write]", OUT)
