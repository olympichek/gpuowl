#!/usr/bin/env python3
"""Parallel runner: one process per shape, torch single-threaded."""
import os, sys, json, time
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
import multiprocessing as mp

def work(args):
    which, total, k, m, j, normone = args
    import torch
    torch.set_num_threads(1)
    torch.set_default_dtype(torch.float64)
    from search2 import BatchModel, target_coeffs, sample_T, run_shape
    nT = 12
    T0, T1 = sample_T(nT, normone=normone, seed=99)
    tgt = torch.stack([target_coeffs(float(a), float(b)) for a, b in zip(T0, T1)])
    t0 = time.time()
    best, mx, model, bi = run_shape(k, m, j, tgt, nT, B=96, steps=14000)
    el = time.time() - t0
    hit = best < 1e-8
    if hit:
        d = {n: p[bi].detach().tolist() for n, p in model.named_parameters()}
        with open(f"sol_{which}_{k}_{m}_{j}.json", "w") as f:
            json.dump(d, f)
    return (k, m, j, best, mx, el)

if __name__ == "__main__":
    which = sys.argv[1]
    total = {"validate9": 9, "eight": 8, "eight-generic": 8, "seven": 7}[which]
    normone = which != "eight-generic"
    if which == "validate9":
        shapes = [(3, 6, 0), (3, 5, 1), (0, 6, 3), (2, 5, 2)]
    else:
        shapes = [(k, m, total - k - m) for k in range(0, total - 1)
                  for m in range(2, total - k + 1) if total - k - m >= 0]
    args = [(which, total, k, m, j, normone) for (k, m, j) in shapes]
    print(f"== {which}: {len(shapes)} shapes in parallel", flush=True)
    with mp.Pool(min(len(shapes), 22)) as pool:
        out = {}
        for (k, m, j, best, mx, el) in pool.imap_unordered(work, args):
            status = "ZERO" if best < 1e-18 else ("NEAR" if best < 1e-8 else "no")
            print(f"shape ({k},{m},{j}): residual {best:.3e}  max|coef| {mx:.1e}  [{status}] ({el:.0f}s)", flush=True)
            out[f"{k},{m},{j}"] = (best, mx)
    with open(f"results_{which}.json", "w") as f:
        json.dump(out, f, indent=1)
    print("done", flush=True)
