# M61 tail pair-square multiplicative-complexity search

Machine search for an `<=8`-wide-product spelling of the PRPLL tail pair-square

```
c = a^2 - T*b^2
d = 2*a*b            a, b in GF(M61^2), T = t^2 the dynamic tail trig
```

in the canonical DGT basis (registry gate in `sol-prpll-speedup-attempts.md`,
"Scaled two-multiply M61 tail rotation audit"). Cost model: every 61-bit wide
product counts (including constant*variable); adds, subtractions, shifts and
`+-i`/conjugation are free; any constant derived from T is free to precompute.

Model: by Strassen homogenization, any straight-line program for the four
homogeneous quadratic targets reduces, with no more counted products, to
k formations (kappa * linear), m bilinears (linear x linear), j scalings
(kappa * quadratic), chained; T-dependence enters only through the kappas.
`search2.py` optimizes the continuous relaxation (combo coefficients free
reals shared across nT simultaneous norm-one T samples; kappas free per T)
with batched Adam restarts + LBFGS polish; `runner.py` fans shapes across
cores; `intensive.py` is the heavy verification pass. nT >= 12 is required:
at nT = 4 the model can interpolate the T samples (see `refine.py`, which
exposed exactly that on an apparent (2,6,0) hit).

Results (2026-08-12, driver run recorded in fable ledger):
- flat bilinear rank of the system (T-dependent coefficients free) = 4,
  via the split of GF(M61^2)[j]/(j^2=-T) (-T is always a QR);
- c-subproblem alone: exactly 6 products (all 5-product shapes fail,
  (3,3,0) hits an exact zero); d alone: 3 (Karatsuba, classical);
- all 28 shapes with k+m+j = 8 fail (best residual 2.1e-4, reproducible
  across 3 independent T-seeds at B=192/28k steps), while the known
  9-product schemes are rediscovered to 1e-30 in every seed;
- verdict: 9 products is optimal in this (fully general) model; the <=8
  gate is closed as unreachable. Both 9-product spellings are already in
  `src/cl/tailsquare.cl` (ENABLE_BETTER/PARALLEL_ONEPAIRSQ) and measured
  tied-or-lost against the production 10 at 300 W.

Caveat: numerical-search evidence (continuous infeasibility), not a proof.
Requires: python3, torch (CPU), numpy.
