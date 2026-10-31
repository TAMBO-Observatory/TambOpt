"""Run 20 quick tests: each test's last layout (opt_state.pt -> layout_final.pt next to it; the final, not the best, layout,
since the best is often still the start) and the validation layout list <run 20>/diagnostics/quick/layouts.txt:
grid_n99, the earlier full-range SGD test (sgdtest_nobg_lr1), the energy-band tests and the construction tests.

    source slurm/run20_guard.sh; [QTAG=quick5] python plots/analysis/run20/export_final_layouts.py
QTAG=quick5: the 5 deg tests (THETA=5 in run20_quick_tests.sh), list in diagnostics/quick5.
"""
import os
import sys

import torch

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, REPO)
import modules.constants as C  # noqa: E402

PRE = os.path.join(C.RUN_LOCATION, "test_v6_run_04_optimize_lbfgs_kernel_")
TAG = os.environ.get("QTAG", "quick")
RUNS = [("sgd_nobg_lr1", "sgdtest_nobg_lr1")] + [
    (f"{TAG}_{k}", f"{TAG}_{k}") for k in ("e500_567", "e567_633", "e633_700", "eall", "slope", "access", "pads",
                                            "contour")]
OUT = os.path.join(C.RUN_LOCATION, "diagnostics", TAG)
os.makedirs(OUT, exist_ok=True)
lines = ["grid_n99=" + os.path.join(C.RUN_LOCATION, "diagnostics", "triplets", "layouts", "grid_n99.pt")]
for key, sfx in RUNS:
    st = torch.load(os.path.join(PRE + sfx, "opt_state.pt"), map_location="cpu", weights_only=False)
    f = os.path.join(PRE + sfx, "layout_final.pt")
    torch.save(dict(x=st["xy"]["x"].detach().float().cpu(), y=st["xy"]["y"].detach().float().cpu(),
                    epoch=st["epoch"]), f)
    lines.append(f"{key}={f}")
    print(f"[write] {f} (step {st['epoch']})")
open(os.path.join(OUT, "layouts.txt"), "w").write("\n".join(lines) + "\n")
print("[write]", os.path.join(OUT, "layouts.txt"))
