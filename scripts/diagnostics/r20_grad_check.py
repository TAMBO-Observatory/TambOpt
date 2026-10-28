"""Run 20 Stage 4 diagnostic: is the layout gradient signal or batch noise, and where is the utility lost?

For each layout (the grid start of Stage 4, and the optimised Poisson-grid layout_best.pt): the swgo_trigger utility and
its (East, North) gradient on --k independent pool batches of --batch showers (Stage 4's pool_utility, noise-free
kernel, no background, as the runs optimised). Reports
  * cosine similarity between batch gradients (pairwise mean, and each batch vs the mean of the others): ~0 = noise;
  * |mean gradient| / mean |batch gradient| (pure noise gives ~1/sqrt(k)) and the per-detector SNR
    |mean| / (std / sqrt(k)) (median and share of detectors above 2);
  * the utility headroom: U / 30 = mean reward; lost to showers that do not trigger (1 - trigger rate) and lost among
    triggered showers (trigger rate - U / 30).
Writes <run>/diagnostics/grad_check.json.

    TAMBOPT_UTILITY=swgo_trigger ... python scripts/diagnostics/r20_grad_check.py [--k 8] [--batch 2048]
"""
import argparse
import importlib.util
import json
import os
import sys

import numpy as np
import torch

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, REPO)
spec = importlib.util.spec_from_file_location("s4", os.path.join(REPO, "scripts", "04_optimize_lbfgs_kernel.py"))
S4 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(S4)
import modules.constants as C  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--k", type=int, default=8)
ap.add_argument("--batch", type=int, default=2048)
ap.add_argument("--chunk", type=int, default=32)
ap.add_argument("--opt", default=os.path.join(C.RUN_LOCATION, "test_v6_run_04_optimize_lbfgs_kernel_swgotrig_l0p3_k3_a2_"
                                              "lr10_grid_b2048", "layout_best.pt"))
args = ap.parse_args()
dev = S4.DEVICE

mountain = S4.load_tr_mountain(C.GEOMETRY_PATH_RESOLVED, C.GEOMETRY_GROUP, C.DET_KEY, east_entry=C.EAST_ENTRY,
                               layer_east_dx=C.LAYER_EAST_DX, n_planes=C.N_PLANES)
surface = S4.SurfaceUpMap.from_mountain(mountain, grid_h=256, grid_w=256).to(dev)
recon = S4.load_recon(C.RECON_FOLDER + "_deepsets", dev)
pool = S4.MainPool(mountain)
n_sweep = pool.n - 1024                                         # Stage 4's sweep part (the scoring set is the tail)
starts = np.random.default_rng(7).choice(n_sweep // args.batch, args.k, replace=False) * args.batch
print(f"[pool] {pool.n} events; {args.k} disjoint batches of {args.batch} at {sorted(starts.tolist())}", flush=True)
batches = [pool.load(int(a), int(a) + args.batch) for a in starts]

e_g, n_g = S4.init_layout(mountain, "grid", dev, C.N_DETECTORS)
opt = torch.load(args.opt, map_location="cpu")
LAYOUTS = {"grid start": (e_g, n_g), "opt Poisson grid": (opt["x"].float().to(dev), opt["y"].float().to(dev))}

res = {}
for name, (e0, n0) in LAYOUTS.items():
    G, rows = [], []
    for b in batches:
        xy = S4.LearnableXY(e0.clone(), n0.clone(), device=str(dev)).to(dev)
        tot = S4.pool_utility(xy, b, surface, recon, args.chunk, backward=True)
        G.append(torch.cat([xy.x.grad, xy.y.grad]).detach().cpu().double())
        rows.append(dict(U=tot["U"], trig=tot["u_phi"] / 30, good_trig=tot["u_e"] / 30))
    G = torch.stack(G)                                            # (k, 2 n_det)
    k = len(G)
    Gn = G / G.norm(dim=1, keepdim=True)
    cos = Gn @ Gn.T
    pair = cos[~torch.eye(k, dtype=bool)]
    loo = [float(torch.nn.functional.cosine_similarity(G[i], (G.sum(0) - G[i]) / (k - 1), dim=0)) for i in range(k)]
    mean = G.mean(0)
    nd = G.shape[1] // 2
    gx, gy = G[:, :nd], G[:, nd:]                                  # per-detector 2-vectors
    m_det = torch.stack([gx.mean(0), gy.mean(0)], 1).norm(dim=1)
    s_det = torch.sqrt((gx.var(0) + gy.var(0)) / 2) / np.sqrt(k)
    snr = (m_det / s_det.clamp_min(1e-30)).numpy()
    U = np.array([r["U"] for r in rows]); trig = np.array([r["trig"] for r in rows])
    res[name] = dict(
        U_mean=float(U.mean()), U_std_between_batches=float(U.std(ddof=1)),
        trigger_rate=float(trig.mean()), mean_reward=float(U.mean() / 30),
        lost_untriggered=float(1 - trig.mean()), lost_among_triggered=float(trig.mean() - U.mean() / 30),
        cos_pairwise_mean=float(pair.mean()), cos_pairwise_std=float(pair.std()),
        cos_leave_one_out=[round(c, 3) for c in loo],
        mean_over_batch_norm=float(mean.norm() / G.norm(dim=1).mean()), pure_noise_expect=float(1 / np.sqrt(k)),
        snr_median=float(np.median(snr)), snr_share_above_2=float((snr > 2).mean()),
        grad_norm_mean=float(G.norm(dim=1).mean()))
    r = res[name]
    print(f"\n[{name}] U = {r['U_mean']:.3f} ± {r['U_std_between_batches']:.3f} (between batches); trigger rate "
          f"{100 * r['trigger_rate']:.1f}%, mean reward {r['mean_reward']:.3f}: lost untriggered "
          f"{100 * r['lost_untriggered']:.1f}%, lost among triggered {100 * r['lost_among_triggered']:.1f}%")
    print(f"  gradient cosine between batches {r['cos_pairwise_mean']:+.3f} ± {r['cos_pairwise_std']:.3f}; "
          f"each vs mean of the others {r['cos_leave_one_out']}")
    print(f"  |mean grad| / mean |grad| = {r['mean_over_batch_norm']:.3f} (pure noise ≈ {r['pure_noise_expect']:.3f}); "
          f"per-detector SNR median {r['snr_median']:.2f}, share > 2: {100 * r['snr_share_above_2']:.0f}%", flush=True)

out = os.path.join(C.RUN_LOCATION, "diagnostics", "grad_check.json")
json.dump(dict(config=vars(args), results=res), open(out, "w"), indent=1)
print("[write]", out)
