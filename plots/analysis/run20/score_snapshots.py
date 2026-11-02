"""Run 20: score the layout snapshots of the sur5snap_<name> runs (slurm/run20_snap_runs.sh) on the kernel: New U with
drawn counts, no background, all read out, on one held-out chunk of 5,000 showers (score_knn.py), against SGD step;
next to the run's own surrogate scoring-set U. Step 1 (--write): write every snapshot as a layout file and the list
<run 20>/diagnostics/snaps/layouts.txt; step 2 (default): plot from <run 20>/diagnostics/snaps/scores.json.
Writes <run 20>/diagnostics/r20b/utility_curves_kernel_snapshots.png.

    source slurm/run20_guard.sh; python plots/analysis/run20/score_snapshots.py --write
    python plots/analysis/run19/triplets/score_knn.py --layouts <snaps>/layouts.txt --windows 0 ... --out <snaps>/scores.json
    python plots/analysis/run20/score_snapshots.py
"""
import argparse
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import torch  # noqa: E402

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, REPO)
import modules.constants as C  # noqa: E402

PRE = os.path.join(C.RUN_LOCATION, "test_v6_run_04_optimize_lbfgs_kernel_")
D = os.path.join(C.RUN_LOCATION, "diagnostics", "snaps")
RUNS = [("eall", "all E"), ("e500_567", "E 5.00–5.67"), ("pads", "+ pads"), ("contour", "+ contour rows")]
ap = argparse.ArgumentParser()
ap.add_argument("--write", action="store_true")
args = ap.parse_args()
os.makedirs(D, exist_ok=True)

if args.write:
    lines = []
    for k, _ in RUNS:
        st = torch.load(PRE + f"sur5snap_{k}/opt_state.pt", map_location="cpu", weights_only=False)
        for step, (x, y) in sorted(st["snaps"].items()):
            f = os.path.join(D, f"{k}_{step:04d}.pt")
            torch.save(dict(x=x.float(), y=y.float()), f)
            lines.append(f"{k}_{step:04d}={f}")
    open(os.path.join(D, "layouts.txt"), "w").write("\n".join(lines) + "\n")
    print("[write]", os.path.join(D, "layouts.txt"), len(lines), "snapshots")
    sys.exit(0)

S = json.load(open(os.path.join(D, "scores.json")))["scores"]
fig, axs = plt.subplots(1, len(RUNS), figsize=(5.2 * len(RUNS), 4.6), constrained_layout=True)
for ax, (k, name) in zip(axs, RUNS):
    pts = sorted((int(n.split("_")[-1]), v["none"]["0"]["u_swgo_drawn"]) for n, v in S.items() if n.startswith(k + "_"))
    ax.plot([p[0] for p in pts], [p[1] for p in pts], color="#e07b39", marker="o", lw=1.8, label="kernel New U")
    ax.axhline(pts[0][1], color="#999999", lw=0.8, ls=":")
    ax.set_xlabel("SGD step")
    ax.set_title(name, fontsize=11)
    ax.grid(alpha=0.25, lw=0.5)
    st = torch.load(PRE + f"sur5snap_{k}/opt_state.pt", map_location="cpu", weights_only=False)
    a2 = ax.twinx()
    a2.plot([s["epoch"] for s in st["score_log"]], [s["U_score"] for s in st["score_log"]], color="#2a6fb0", lw=1.2,
            alpha=0.8, label="surrogate scoring U (− penalty)")
    a2.tick_params(labelsize=8, colors="#2a6fb0")
    if k == RUNS[0][0]:
        ax.set_ylabel("kernel New U (5,000 held-out showers)")
        h1, l1 = ax.get_legend_handles_labels(); h2, l2 = a2.get_legend_handles_labels()
        ax.legend(h1 + h2, l1 + l2, fontsize=8, frameon=False, loc="lower left")
fig.suptitle("Surrogate optimisation (5° utility, SGD from the 100-detector grid) measured on the kernel along the path: "
             "New U (drawn counts, no background, all read out) every 100 steps; dotted = start", fontsize=11)
out = os.path.join(C.RUN_LOCATION, "diagnostics", "r20b", "utility_curves_kernel_snapshots.png")
fig.savefig(out, dpi=150, facecolor="white")
print("[write]", out)
