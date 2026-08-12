#!/usr/bin/env python3
"""Flat bilinear-rank probe for the M61 tail pair-square system.

System (x = (a0,a1,b0,b1)), targets as quadratic forms x^T Q_k x:
  c0 = a0^2 - a1^2 - t0(b0^2-b1^2) + 2 t1 b0 b1
  c1 = 2 a0 a1 - t1(b0^2-b1^2) - 2 t0 b0 b1
  d0 = 2(a0 b0 - a1 b1)
  d1 = 2(a0 b1 + a1 b0)

Flat model ("scalar mults free"): r products P_s = (u_s.x)(v_s.x) with fully
generic real coefficients; outputs are generic combinations sum_s alpha[k,s] P_s.
Question: minimal r. Known upper bound 6. ALS with restarts.
"""
import numpy as np
import sys

rng = np.random.default_rng(12345)

def targets(t0, t1):
    Q = np.zeros((4, 4, 4))
    # c0
    Q[0][0, 0] = 1; Q[0][1, 1] = -1
    Q[0][2, 2] = -t0; Q[0][3, 3] = t0
    Q[0][2, 3] = Q[0][3, 2] = t1
    # c1
    Q[1][0, 1] = Q[1][1, 0] = 1
    Q[1][2, 2] = -t1; Q[1][3, 3] = t1
    Q[1][2, 3] = Q[1][3, 2] = -t0
    # d0
    Q[2][0, 2] = Q[2][2, 0] = 1
    Q[2][1, 3] = Q[2][3, 1] = -1
    # d1
    Q[3][0, 3] = Q[3][3, 0] = 1
    Q[3][1, 2] = Q[3][2, 1] = 1
    return Q

def sym(u, v):
    return 0.5 * (np.outer(u, v) + np.outer(v, u))

def als(Q, r, iters=400, seed=0):
    rg = np.random.default_rng(seed)
    U = rg.normal(size=(r, 4)); V = rg.normal(size=(r, 4))
    A = None
    tvec = Q.reshape(4, 16)  # targets flattened
    for it in range(iters):
        # solve alpha: least squares over span of S_s
        S = np.array([sym(U[s], V[s]).ravel() for s in range(r)])  # r x 16
        A, *_ = np.linalg.lstsq(S.T, tvec.T, rcond=None)  # r x 4
        A = A.T  # 4 x r
        # solve U (residual linear in U): for each k: sum_s a_ks * sym(u_s,v_s)
        # build M: (4*16) x (r*4); vec(sym(u,v)) = 0.5(kron(v,u)+kron(u,v)) -> linear in u
        M = np.zeros((4 * 16, r * 4))
        for k in range(4):
            for s in range(r):
                # d vec(sym(u_s,v_s))/du_s = 0.5*(I kron v + v kron I) pattern
                B = np.zeros((16, 4))
                for i in range(4):
                    for j in range(4):
                        row = i * 4 + j
                        B[row, i] += 0.5 * V[s][j]
                        B[row, j] += 0.5 * V[s][i]
                M[k * 16:(k + 1) * 16, s * 4:(s + 1) * 4] = A[k, s] * B
        sol, *_ = np.linalg.lstsq(M, tvec.ravel(), rcond=None)
        U = sol.reshape(r, 4)
        # solve V symmetrically
        M = np.zeros((4 * 16, r * 4))
        for k in range(4):
            for s in range(r):
                B = np.zeros((16, 4))
                for i in range(4):
                    for j in range(4):
                        row = i * 4 + j
                        B[row, i] += 0.5 * U[s][j]
                        B[row, j] += 0.5 * U[s][i]
                M[k * 16:(k + 1) * 16, s * 4:(s + 1) * 4] = A[k, s] * B
        sol, *_ = np.linalg.lstsq(M, tvec.ravel(), rcond=None)
        V = sol.reshape(r, 4)
    S = np.array([sym(U[s], V[s]).ravel() for s in range(r)])
    A, *_ = np.linalg.lstsq(S.T, tvec.T, rcond=None)
    R = tvec.T - S.T @ A
    return np.sqrt((R ** 2).sum()), U, V, A.T

def probe(r, t0, t1, tries=60):
    Q = targets(t0, t1)
    best = (np.inf, None)
    for s in range(tries):
        res, U, V, A = als(Q, r, iters=300, seed=s * 7919 + r)
        if res < best[0]:
            best = (res, (U, V, A))
        if res < 1e-9:
            break
    return best

if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "both"
    # generic T and norm-one T
    samples = []
    if mode in ("both", "generic"):
        samples.append(("generic", 0.7345298, -1.318809))
    if mode in ("both", "normone"):
        s = 0.4183321
        samples.append(("norm-one", (1 - s * s) / (1 + s * s), 2 * s / (1 + s * s)))
    for name, t0, t1 in samples:
        print(f"--- T {name}: t0={t0:.6f} t1={t1:.6f} (t0^2+t1^2={t0*t0+t1*t1:.6f})")
        for r in (4, 5, 6):
            res, sol = probe(r, t0, t1)
            print(f"  rank {r}: best residual {res:.3e}  -> {'ACHIEVABLE' if res < 1e-8 else 'not found'}")
