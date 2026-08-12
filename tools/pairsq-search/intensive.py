#!/usr/bin/env python3
"""Intensive verification runs on the most promising 8-shapes, plus controls.

Jobs:
  - top 8-shapes x 3 T-seeds, B=192, steps=28000 (norm-one T, nT=12)
  - control: known-feasible 9-shape (3,6,0) under identical settings x 3 seeds
  - c-alone subproblem (targets c0,c1 only) at totals 5 and 6
"""
import os, sys, json, time
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
import multiprocessing as mp

def work(args):
    label, k, m, j, seed, subset = args
    import torch
    torch.set_num_threads(1)
    torch.set_default_dtype(torch.float64)
    import search2
    from search2 import BatchModel, target_coeffs, sample_T, run_shape
    nT = 12
    T0, T1 = sample_T(nT, normone=True, seed=seed)
    tgt = torch.stack([target_coeffs(float(a), float(b)) for a, b in zip(T0, T1)])
    if subset == "c":
        tgt = tgt[:, :2, :]  # c0, c1 only
        # monkeypatch: W has 4 output rows; restrict by comparing only 2 rows
        # simplest: pad targets back with zeros and zero-weight d rows -> instead
        # just use 2-row W via subclass: easier to slice model output
    t0 = time.time()
    if subset == "c":
        best, mx = run_c_only(k, m, j, tgt, nT, seed)
        return (label, k, m, j, seed, best, mx, time.time() - t0)
    best, mx, model, bi = run_shape(k, m, j, tgt, nT, B=192, steps=28000)
    if best < 1e-10:
        d = {n: p[bi].detach().tolist() for n, p in model.named_parameters()}
        with open(f"sol_intensive_{k}_{m}_{j}_s{seed}.json", "w") as f:
            json.dump(d, f)
    return (label, k, m, j, seed, best, mx, time.time() - t0)

def run_c_only(k, m, j, tgt, nT, seed):
    import torch
    from search2 import BatchModel, polish
    torch.manual_seed(seed)
    B = 128
    model = BatchModel(B, k, m, j, nT)
    # shrink W to 2 outputs
    with torch.no_grad():
        model.W = torch.nn.Parameter(model.W[:, :2, :].clone())
    opt = torch.optim.Adam(model.parameters(), lr=0.05)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=20000, eta_min=1e-4)
    tt = tgt.unsqueeze(0)
    for step in range(20000):
        opt.zero_grad()
        res = (model() - tt).pow(2).sum(dim=(1, 2, 3))
        res.sum().backward()
        opt.step(); sched.step()
    with torch.no_grad():
        res = (model() - tt).pow(2).sum(dim=(1, 2, 3))
    order = res.argsort()
    import math
    best = math.inf
    for cand in order[:3]:
        bi = int(cand)
        one = BatchModel(1, k, m, j, nT)
        with torch.no_grad():
            one.W = torch.nn.Parameter(one.W[:, :2, :].clone())
            for (n, p), (n2, q) in zip(one.named_parameters(), model.named_parameters()):
                p.copy_(q[bi:bi + 1])
        r = polish(one, tgt)
        best = min(best, r)
        if best < 1e-20:
            break
    mx = 0.0
    return best, mx

if __name__ == "__main__":
    jobs = []
    for (k, m, j) in [(2, 5, 1), (3, 5, 0), (2, 6, 0), (1, 5, 2), (3, 4, 1)]:
        for seed in (99, 7, 31337):
            jobs.append(("eight", k, m, j, seed, "full"))
    for seed in (99, 7, 31337):
        jobs.append(("nine-ctrl", 3, 6, 0, seed, "full"))
    # c-alone: shapes with totals 5 and 6
    for (k, m, j) in [(2, 3, 0), (2, 2, 1), (1, 3, 1), (3, 2, 0), (1, 4, 0), (0, 3, 2)]:
        jobs.append(("c5" if k + m + j == 5 else "c6", k, m, j, 99, "c"))
    for (k, m, j) in [(3, 3, 0), (2, 4, 0), (2, 3, 1), (1, 4, 1), (3, 2, 1), (0, 4, 2)]:
        jobs.append(("c6", k, m, j, 99, "c"))
    print(f"{len(jobs)} jobs", flush=True)
    with mp.Pool(22) as pool:
        for (label, k, m, j, seed, best, mx, el) in pool.imap_unordered(work, jobs):
            status = "ZERO" if best < 1e-15 else ("NEAR" if best < 1e-8 else "no")
            print(f"[{label}] ({k},{m},{j}) seed={seed}: {best:.3e} [{status}] ({el:.0f}s)", flush=True)
