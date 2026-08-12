#!/usr/bin/env python3
"""v2: batched continuous-feasibility search over homogeneous program shapes.

By Strassen homogenization, ANY straight-line program (adds/shifts free, every
wide product counted, arbitrary precomputed T-constants free) computing the
four homogeneous quadratic targets converts, with no more counted products,
into exactly:
  k formations  f_i = kappa_i * L_i(vars, f_<i)        [deg-0 x deg-1]
  m bilinears   p_i = L(vars, f_*) * R(vars, f_*)      [deg-1 x deg-1]
  j scalings    s_i = kappa'_i * C_i(p_*, s_<i)        [deg-0 x deg-2]
  outputs       = combos of (p_*, s_*)
(deg0xdeg0 products become free precomputed constants.)

Phase 1 here: CONTINUOUS relaxation (all combo coefficients free reals,
kappas free per T sample). If no shape with k+m+j = 8 reaches residual 0,
no 8-product program exists over the reals with generic norm-one T --
a genuine algebraic infeasibility argument (modulo optimizer completeness,
mitigated by many restarts + validation at 9 where zeros must exist).

Batched over restarts. Tracks coefficient blowup to flag border-rank escapes.
"""
import torch, math, sys, time, json

torch.set_default_dtype(torch.float64)
torch.manual_seed(0)

NV = 4
MON = [(i, j) for i in range(NV) for j in range(i, NV)]
NM = len(MON)

def target_coeffs(T0, T1):
    t = torch.zeros(4, NM)
    def setm(k, i, j, v):
        t[k, MON.index((min(i, j), max(i, j)))] += v
    setm(0, 0, 0, 1); setm(0, 1, 1, -1); setm(0, 2, 2, -T0); setm(0, 3, 3, T0); setm(0, 2, 3, 2 * T1)
    setm(1, 0, 1, 2); setm(1, 2, 2, -T1); setm(1, 3, 3, T1); setm(1, 2, 3, -2 * T0)
    setm(2, 0, 2, 2); setm(2, 1, 3, -2)
    setm(3, 0, 3, 2); setm(3, 1, 2, 2)
    return t

def sample_T(n, normone=True, seed=0):
    g = torch.Generator().manual_seed(seed)
    if normone:
        s = torch.rand(n, generator=g) * 4 - 2
        T0 = (1 - s * s) / (1 + s * s); T1 = 2 * s / (1 + s * s)
    else:
        T0 = torch.rand(n, generator=g) * 4 - 2
        T1 = torch.rand(n, generator=g) * 4 - 2
    return T0, T1

class BatchModel(torch.nn.Module):
    """B independent restarts; nT T-samples shared."""
    def __init__(self, B, k, m, j, nT, scale=0.7):
        super().__init__()
        def P(*shape):
            return torch.nn.Parameter(torch.randn(*shape) * scale)
        self.B, self.k, self.m, self.j, self.nT = B, k, m, j, nT
        # formation i combo over [vars, f_<i]: ragged; use square mask
        self.F = P(B, k, NV + k) if k else None      # mask applied for f_<i
        self.kf = P(B, nT, k) if k else None
        self.OL = P(B, m, NV + k)
        self.OR = P(B, m, NV + k)
        self.S = P(B, j, m + j) if j else None       # combo over [p_*, s_<i], masked
        self.ks = P(B, nT, j) if j else None
        self.W = P(B, 4, m + j)
        if k:
            fm = torch.zeros(k, NV + k)
            fm[:, :NV] = 1
            for i in range(k):
                fm[i, NV:NV + i] = 1
            self.register_buffer('fmask', fm)
        if j:
            sm = torch.zeros(j, m + j)
            sm[:, :m] = 1
            for i in range(j):
                sm[i, m:m + i] = 1
            self.register_buffer('smask', sm)

    def forward(self):
        B, k, m, j, nT = self.B, self.k, self.m, self.j, self.nT
        # build linear quantities: (B, nT, NV + k, NV)
        base = torch.eye(NV).view(1, 1, NV, NV).expand(B, nT, NV, NV)
        if k:
            lins = [base[:, :, i] for i in range(NV)]  # each (B, nT, NV)
            F = self.F * self.fmask                    # (B, k, NV+k)
            for i in range(k):
                combo = sum(F[:, i, q].view(B, 1, 1) * lins[q] for q in range(NV + i))
                fi = self.kf[:, :, i].unsqueeze(-1) * combo
                lins.append(fi)
            lin = torch.stack(lins, dim=2)             # (B, nT, NV+k, NV)
        else:
            lin = base
        L = torch.einsum('bmq,btqv->btmv', self.OL, lin)
        R = torch.einsum('bmq,btqv->btmv', self.OR, lin)
        P = torch.zeros(B, nT, m, NM)
        for idx, (i, jj) in enumerate(MON):
            if i == jj:
                P[:, :, :, idx] = L[:, :, :, i] * R[:, :, :, i]
            else:
                P[:, :, :, idx] = L[:, :, :, i] * R[:, :, :, jj] + L[:, :, :, jj] * R[:, :, :, i]
        quads = [P[:, :, i] for i in range(m)]         # each (B, nT, NM)
        if j:
            S = self.S * self.smask
            for i in range(j):
                combo = sum(S[:, i, q].view(B, 1, 1) * quads[q] for q in range(m + i))
                si = self.ks[:, :, i].unsqueeze(-1) * combo
                quads.append(si)
        Q = torch.stack(quads, dim=2)                  # (B, nT, m+j, NM)
        out = torch.einsum('bkp,btpn->btkn', self.W, Q)
        return out

def extract_restart(model, bi):
    one = BatchModel(1, model.k, model.m, model.j, model.nT)
    with torch.no_grad():
        for (n, p), (n2, q) in zip(one.named_parameters(), model.named_parameters()):
            p.copy_(q[bi:bi + 1])
    return one

def polish(model, tgt, iters=800):
    tt = tgt.unsqueeze(0)
    # second-stage Adam at low lr first
    opt = torch.optim.Adam(model.parameters(), lr=2e-3)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=5000, eta_min=1e-6)
    for _ in range(5000):
        opt.zero_grad()
        loss = (model() - tt).pow(2).sum()
        loss.backward()
        opt.step(); sched.step()
    # scaled-loss LBFGS to sharpen true zeros
    opt = torch.optim.LBFGS(model.parameters(), max_iter=iters, tolerance_grad=1e-14,
                            tolerance_change=1e-20, history_size=50, line_search_fn='strong_wolfe')
    def closure():
        opt.zero_grad()
        loss = (model() - tt).pow(2).sum() * 1e8
        loss.backward()
        return loss
    opt.step(closure)
    with torch.no_grad():
        return (model() - tt).pow(2).sum().item()

def run_shape(k, m, j, tgt, nT, B=48, steps=24000, lr0=0.05, tag=""):
    model = BatchModel(B, k, m, j, nT)
    opt = torch.optim.Adam(model.parameters(), lr=lr0)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=steps, eta_min=1e-4)
    tt = tgt.unsqueeze(0)
    for step in range(steps):
        opt.zero_grad()
        res = (model() - tt).pow(2).sum(dim=(1, 2, 3))   # per restart
        res.sum().backward()
        opt.step(); sched.step()
    with torch.no_grad():
        res = (model() - tt).pow(2).sum(dim=(1, 2, 3))
    order = res.argsort()
    best, bestm, bi = math.inf, None, order[0]
    for cand in order[:3]:
        one = extract_restart(model, int(cand))
        r = polish(one, tgt)
        if r < best:
            best, bestm, bi = r, one, int(cand)
        if best < 1e-20:
            break
    mx = max(p.abs().max().item() for p in bestm.parameters())
    return best, mx, bestm, 0

if __name__ == "__main__":
    which = sys.argv[1]
    nT = 4
    T0, T1 = sample_T(nT, normone=(which != "eight-generic"), seed=99)
    tgt = torch.stack([target_coeffs(float(a), float(b)) for a, b in zip(T0, T1)])
    total = {"validate9": 9, "eight": 8, "eight-generic": 8, "seven": 7}[which]
    if which == "validate9":
        shapes = [(3, 6, 0), (3, 5, 1), (0, 6, 3), (2, 5, 2)]
    else:
        shapes = [(k, m, total - k - m) for k in range(0, total - 1)
                  for m in range(2, total - k + 1) if total - k - m >= 0]
    print(f"== {which}: total={total}, {len(shapes)} shapes, nT={nT}")
    results = {}
    for (k, m, j) in shapes:
        t0 = time.time()
        best, mx, model, bi = run_shape(k, m, j, tgt, nT)
        results[f"{k},{m},{j}"] = (best, mx)
        status = "ZERO" if best < 1e-18 else ("NEAR" if best < 1e-8 else "no")
        print(f"shape ({k},{m},{j}): residual {best:.3e}  max|coef| {mx:.1e}  [{status}] ({time.time()-t0:.0f}s)", flush=True)
        if best < 1e-8:
            d = {n: p[bi].detach().tolist() for n, p in model.named_parameters()}
            with open(f"sol_{which}_{k}_{m}_{j}.json", "w") as f:
                json.dump(d, f)
    with open(f"results_{which}.json", "w") as f:
        json.dump(results, f, indent=1)
