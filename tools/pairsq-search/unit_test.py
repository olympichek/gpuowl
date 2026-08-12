#!/usr/bin/env python3
"""Unit test: hard-code the known 9-product scheme into the Model and check
zero residual; then diagnose continuous-only optimization at shape (3,6,0)."""
import torch
from struct_search import Model, target_coeffs, sample_T, try_shape

torch.set_default_dtype(torch.float64)

nT = 4
T0, T1 = sample_T(nT, normone=True, seed=99)
tgt = torch.stack([target_coeffs(float(a), float(b)) for a, b in zip(T0, T1)])

m = Model(3, 6, 0, nT, seed=0)
with torch.no_grad():
    # formations: f1 = T0*b0, f2 = T1*b1, f3 = (T0+T1)*(b0+b1)
    m.F.copy_(torch.tensor([
        [0., 0., 1., 0.],
        [0., 0., 0., 1.],
        [0., 0., 1., 1.],
    ]))
    # half-root t = sqrt(T): T norm-one real embedding T=(cos,sin) -> t=(cos/2 angle)
    t0h = torch.sqrt((1 + T0) / 2)
    t1h = T1 / (2 * t0h)
    m.kf.copy_(torch.stack([t0h, t1h, t0h + t1h], dim=1))
    # w0 = f1 - f2, w1 = f3 - f1 - f2  (as combos over basis [a0,a1,b0,b1,f1,f2,f3])
    # bilinears:
    # K1 = (a0+w0)(a0-w0); K2 = (a1+w1)(a1-w1); K3 = (a0+a1+w0+w1)(a0+a1-w0-w1)
    # K4 = a0*b0; K5 = a1*b1; K6 = (a0+a1)(b0+b1)
    OL = torch.zeros(6, 7); OR = torch.zeros(6, 7)
    def combo(a0=0., a1=0., b0=0., b1=0., f1=0., f2=0., f3=0.):
        return torch.tensor([a0, a1, b0, b1, f1, f2, f3])
    w0 = combo(f1=1., f2=-1.); w1 = combo(f1=-1., f2=-1., f3=1.)
    a0 = combo(a0=1.); a1 = combo(a1=1.); b0 = combo(b0=1.); b1 = combo(b1=1.)
    OL[0] = a0 + w0; OR[0] = a0 - w0
    OL[1] = a1 + w1; OR[1] = a1 - w1
    OL[2] = a0 + a1 + w0 + w1; OR[2] = a0 + a1 - w0 - w1
    OL[3] = a0; OR[3] = b0
    OL[4] = a1; OR[4] = b1
    OL[5] = a0 + a1; OR[5] = b0 + b1
    m.OL.copy_(OL); m.OR.copy_(OR)
    # outputs: c0 = K1-K2; c1 = K3-K1-K2; d0 = 2K4-2K5; d1 = 2K6-2K4-2K5
    m.W.copy_(torch.tensor([
        [1., -1., 0., 0., 0., 0.],
        [-1., -1., 1., 0., 0., 0.],
        [0., 0., 0., 2., -2., 0.],
        [0., 0., 0., -2., -2., 2.],
    ]))
res = (m() - tgt).pow(2).sum().item()
print(f"hard-coded 9-scheme residual: {res:.3e}  ({'PASS' if res < 1e-20 else 'FAIL'})")

# continuous-only diagnosis
model = Model(3, 6, 0, nT, seed=42)
opt = torch.optim.Adam(model.parameters(), lr=0.03)
for step in range(6000):
    opt.zero_grad()
    loss = (model() - tgt).pow(2).sum()
    loss.backward()
    opt.step()
print(f"continuous-only (3,6,0) after 6000 steps: {loss.item():.3e}")
