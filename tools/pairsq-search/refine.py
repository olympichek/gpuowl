#!/usr/bin/env python3
"""Robustness + structure probes for the (2,6,0) 8-product hit.

1. Freeze discrete structure from sol_eight_2_6_0.json; per fresh T refit ONLY
   the 2 kappas (2 DOF vs 40 target coefficients): exact fit on fresh norm-one
   T proves the scheme is generic, not a 4-sample interpolation.
   Also test generic (non-norm-one) T: does the scheme need |T| = 1?
2. Canonical-F probe: force f1 = k1*b0, f2 = k2*(b1 + gamma*f1)? simplest:
   f1 = k1*b0, f2 = k2*b1; reoptimize everything else from scratch.
"""
import torch, json, math, sys
torch.set_num_threads(8)
torch.set_default_dtype(torch.float64)
from search2 import BatchModel, target_coeffs, sample_T, polish

def load_model(path, nT):
    d = json.load(open(path))
    m = BatchModel(1, 2, 6, 0, nT)
    with torch.no_grad():
        m.F.copy_(torch.tensor(d['F']).unsqueeze(0))
        m.OL.copy_(torch.tensor(d['OL']).unsqueeze(0))
        m.OR.copy_(torch.tensor(d['OR']).unsqueeze(0))
        m.W.copy_(torch.tensor(d['W']).unsqueeze(0))
        m.kf.copy_(torch.zeros(1, nT, 2))
    return m

def fit_kappas_per_T(m, tgt, iters=6000):
    for n, p in m.named_parameters():
        p.requires_grad_(n == 'kf')
    opt = torch.optim.Adam([m.kf], lr=0.05)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=iters, eta_min=1e-5)
    tt = tgt.unsqueeze(0)
    for _ in range(iters):
        opt.zero_grad()
        loss = (m() - tt).pow(2).sum()
        loss.backward()
        opt.step(); sched.step()
    # polish kappas with LBFGS
    opt = torch.optim.LBFGS([m.kf], max_iter=300, tolerance_grad=1e-16,
                            tolerance_change=1e-20, line_search_fn='strong_wolfe')
    def cl():
        opt.zero_grad()
        loss = (m() - tt).pow(2).sum() * 1e8
        loss.backward()
        return loss
    opt.step(cl)
    with torch.no_grad():
        per_t = (m() - tt).pow(2).sum(dim=(2, 3)).squeeze(0)
    return per_t

if __name__ == "__main__":
    step = sys.argv[1] if len(sys.argv) > 1 else "fresh"
    if step == "fresh":
        for name, normone, seed in [("norm-one", True, 12345), ("generic", False, 54321)]:
            nT = 12
            T0, T1 = sample_T(nT, normone=normone, seed=seed)
            tgt = torch.stack([target_coeffs(float(a), float(b)) for a, b in zip(T0, T1)])
            m = load_model('sol_eight_2_6_0.json', nT)
            per_t = fit_kappas_per_T(m, tgt)
            print(f"fresh {name} T: per-T residuals:")
            for i in range(nT):
                print(f"  T=({float(T0[i]):+.4f},{float(T1[i]):+.4f})  res={float(per_t[i]):.3e}")
    elif step == "canonF":
        nT = 6
        T0, T1 = sample_T(nT, normone=True, seed=777)
        tgt = torch.stack([target_coeffs(float(a), float(b)) for a, b in zip(T0, T1)])
        best = math.inf
        for trial in range(24):
            m = BatchModel(1, 2, 6, 0, nT)
            with torch.no_grad():
                m.F.copy_(torch.tensor([[[0., 0., 1., 0., 0., 0.],
                                         [0., 0., 0., 1., 0., 0.]]]))
            m.F.requires_grad_(False)
            opt = torch.optim.Adam([p for n, p in m.named_parameters() if n != 'F'], lr=0.03)
            sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=14000, eta_min=1e-4)
            tt = tgt.unsqueeze(0)
            for _ in range(14000):
                opt.zero_grad()
                loss = (m() - tt).pow(2).sum()
                loss.backward()
                opt.step(); sched.step()
            r = polish(m, tgt)
            best = min(best, r)
            print(f"  trial {trial}: {r:.3e} (best {best:.3e})", flush=True)
            if best < 1e-20:
                d = {n: p.detach().tolist() for n, p in m.named_parameters()}
                json.dump(d, open('sol_canonF.json', 'w'))
                break
        print(f"canonical F (f1=k1*b0, f2=k2*b1): best {best:.3e}")
