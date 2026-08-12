# Fable campaign: PRPLL speedup attempts

Continuation of the architectural speedup campaign recorded in
[`sol-prpll-speedup-attempts.md`](sol-prpll-speedup-attempts.md).  That file is
the authoritative registry of everything already tried; nothing recorded there
as rejected may be repeated without satisfying its stated reopen condition.
This file tracks the new campaign: ideas to try, experiments run, and lessons.

## Inherited state (verified 2026-08-12)

- Hardware: NVIDIA RTX PRO 6000 Blackwell Max-Q (GB202, sm_120), 97.9 GiB,
  128 MiB L2, 300 W enforced power limit (325 W max requires unavailable
  permissions), driver 595.71.05, CUDA 13.2.  No Nsight Compute counters
  (`ERR_NVGPUCTRPERM`).  Workload is integer-ALU/power limited, not
  bandwidth limited; INT32 and FP32 share unified pipes (no free FP32
  reservoir); no FP64 tensor MMA on sm_120.
- Production: exponent 136279841, `-fft 1:512:8:512:202` (4M M31+M61), tuned
  `-use INPLACE=1,LOADS=22042,STORES=21,TABMUL_CHAIN32=1,MODM31=2,MULTI_Q=1,L1CUDA=3,GRAPHS=0`.
  Reverified today: **201.5 us/iteration** over 100k iterations,
  iteration-2,000 residue `05d6515c416b83e2`, iteration-100,000 residue
  `52775eea4730be87` — identical to Sol's records.  Long 1M-run reference:
  205.5 us.
- Kernel medians (Nsight timeline): carryFused 71-74, M61 tail 56.1,
  M61 middle-out 28.8, M61 middle-in 27.0 (M61 bottom half 111.9);
  M31 bottom half 65.6 total, fully hidden under M61.  Ideal floor with the
  current kernel set: 111.9 + 73.8 = 185.7 us.
- Sol's acceptance gate was **exact, sustained <= 180 us/iteration over 1M
  iterations** with a verifiable residue.  ~90 experiments later, no exact
  end-to-end result beat even the 205.5-us production number.

## Registry of rejected idea classes (do not repeat; see Sol's file for details)

1. **Alternative fields/moduli** — all rejected at arithmetic, root, capacity,
   or complete-budget gates: 3x31-bit direct RNS; 2M five-prime; 2M
   M31+M61+q0+q1; sub-2^30 Riesel redundant arithmetic; packed Riesel pairs
   (M31R2); M127; M89; M31^3; M31-prime-power ring; two-61-bit-field 3M;
   Goldilocks (scalar, packed, folded, 2M); 86-91-bit single scalar field;
   near-2^53 Riesel; q53, q45, q55, q24, qC=7*2^21-1, q15/q7; M19 (3M
   M31+M61+M19 at complete-budget; radix-7 3.5M); M29/M43 (composite);
   quartic M17; composite rings M31*M61, M31*M19, M31*q53; all-Mersenne
   M31/M43/M19; scalar-prime family closed except Goldilocks; exhaustive u32
   Proth root search empty; negacyclic root-of-minus-two algebraically void.
2. **FP32/FP64 planes** — FP32+M61 wide quotient is *169.4 us but inexact at
   p150*: dense quotient errors in [-120,120] (~589k pairs/iteration err>=1),
   leaving a 10.6-us exact-repair budget that nothing fits.  Rejected repairs:
   compensated FP32 (2-3x cost), double-single, half residual, stage- and
   operation-class-selective compensation, input splitting, block-scaled Q31,
   FP32-limb M31 (3.1x), FP64 backend/8M, FP64-in-M61 FMA offload, carry-aware
   lifting (algebraic no-go: carry is normalization, not a congruence).
3. **Folded/syndrome correction sidecars** — N/8, N/16, factor-4 with selected
   parity bit (decoder proofs pass; transforms cost 14-64 us vs 10.6-us
   budget), ballot parity, ternary decoders, Z/256Z Nussbaumer sidecar,
   byte-ring, CLMAD, factor-two 24-bit Riesel, Tensor sidecars (contention).
4. **Tensor cores** — dense INT8 (42% slower), TF32, Toeplitz M61, binary
   b1 AND/POPC (41-73x slower), 2:4 sparsity (algebraic no-go), Toom tiles.
5. **Fusion/residency/clusters** — middle/width/carry fusion (311 us);
   M61-only 64-KiB resident tile (+2-3%); DSM clusters 2/4/8-block and
   line-pair (all lose, up to 2.2x); clustered width-to-middle (+8..64 us);
   width/middle factorization swap (rotation count conserved, loses 7-8 us);
   carry-to-middle ready-queue (ownership audit: only 12% overlap possible);
   one-warp 512-point tail; TMA rewrite of unchanged transform (no evidence);
   warp specialization (hardware no-go).
6. **Scheduling/host** — multi-exponent batching (-13..-35%); batch-native
   carry; stream priorities (-4 us); PDL (upstream tested, no help); CUDA
   graphs (~1-2 us only, GRAPHS=0 tuned); persisting L2 (no effect); L2
   striping; LL recurrence (~206 us, no); two workers (see tune notes).
7. **Representation micro-changes** — persistent radix-2^31 M61 limbs (2.9x);
   fused/whole-radix lazy limbs; Montgomery state; twiddle compression; full
   M61 middle-root planes (isolated +7.6-8.3%, production -0.7..-0.9%);
   quotient-domain carry; delayed normalization; alternate/parallel pair-square
   identities; scaled two-multiply rotation (algebra); spectral preweight;
   Hartley; circle FFT; Nussbaumer polynomial ring; DGT folds; radix-2 width
   ownership (NW=2: +31%); split-radix/Winograd audits (<=5.6% nominal, no).
8. **Transform shapes** — 8M two-prime; 3M/3.5M/3.875M/3.9375M/4.0625M/4.5M
   odd-radix architectures (radix 3/7/9/13/31/33/63/65 edges all fail on
   ownership-complete odd-DFT cost or capacity); PFA9 bridge; 33/32 near-M61;
   further equivalent 4M factorization sweeps closed.
9. **Hardware controls** — power limit locked at 300 W; app clocks deprecated;
   no NCU counters.  cuFFT exact-convolution lower bound rejected.

Key structural facts driving everything (REVISED by the Lead-1 measurement):
- Iteration = ~71.5-us serial fused carry + ~130-us bottom-half window in
  which the M61 chain (112 us pure) and M31 chain (65.6 us pure) co-run.
- **There is no idle execution reservoir.**  Under concurrency the M31
  kernels dilate to ~108 dilated us and the M61 chain to ~130; the window is
  bound by total work, not by the longer chain alone.  Sol's "M31 is hidden"
  held only for kernel medians in isolation.
- Consequences: (a) relocating work between queues or out of the fused carry
  is paid at nearly full price; (b) only *deleting* work (instructions or
  bytes) shortens the iteration; (c) adding any sidecar work costs
  proportionally — which retroactively explains the uniform failure of every
  sidecar/fold/batching experiment in Sol's campaign.
- The serial fused carry co-runs with nothing, and both boundary-overlap
  directions (carry->middle, middle-out->carry) are closed by ownership
  audits (~12% max overlap).

## Active leads (Fable)

### Lead 1: exact detached M31 width kernels (ACTIVE)

Sol's `DETACHED_M31_EDGE=1` timing scaffold (correctness-disabled) measured
**186.3 us vs 205.1 us matched control** on the exact production path by
removing both M31 width transforms from `carryFused` (it keeps the flat
loads/stores, so the carry-kernel traffic is already charged).  Sol rejected
it because even the free-transform bound misses the 180-us gate — but it was
never completed, and ~15-19 us is still a real end-to-end speedup if the two
external M31 width kernels can hide in the ~46-us M31-stream slack:

- Kernel A (`fftWIn31`): reads `in31` in carryFused line layout, performs
  `fft_WIDTH1` with `smallTrig31`, writes flat
  `[(line*NW+i)*G_W+lowMe]` — exactly what detached carryFused reads.
  Runs on the M31 queue after `fftMiddleOut31`, hidden under M61 tail/out.
- Kernel B (`fftWOut31`): reads flat pre-width residues written by detached
  carryFused, performs `fft_WIDTH2`, `writeCarryFusedLine` to `out31`.
  Runs on the M31 queue right after carryFused, hidden under M61 middle-in.
- Scratch: one 16-MiB GF31 buffer for A-out; carryFused' writes B-in flat
  (can reuse `in31`, dead after A consumes it — verify replay ordering).
- Extra traffic vs scaffold: A+B ~64 MiB/iteration through L2 on the hidden
  stream (~0.75 TB/s over the slack window; aggregate L2 demand stays well
  under the ~7-8 TB/s ceiling).
- Projection: 186.3 + few us contention => ~189-194 vs 201.5 control.
  Success gate: exact residues at 2k/10k/20k/40k/50k checkpoints AND
  a reproducible >= 5-us end-to-end win in matched alternating runs.
- Risks: A/B contention with M61 kernels (power + issue slots), sync/launch
  slop on the main queue, layout/padding subtleties under INPLACE=1.

### Lead 2: config.txt's suggested tune toggles — CLOSED (no gain)

Matched alternating 50k runs of
`TAIL_TRIGS31=1,UNROLL_W=0,UNROLL_H=0,ZEROHACK_H=0` (the "slightly faster in
-tune" hints) vs baseline: trial blocks 199.2-207.0 vs adjacent controls
197.9-205.6 — every trial at or slightly above its neighboring control, with
visible thermal drift across the batch.  No sustained gain; residues exact.
Do not adopt.

### Lead 3: partial-edge overlap — CLOSED by Lead 1's lesson

The premise (usable M31-stream slack) is invalid: both queues are
work-saturated under concurrency; only the ~22-us pre-merge idle exists and
mode 2 showed even that absorbs nothing usable.

### Lead 4: full -tune ntt rerun on current driver

The complete NTT option sweep on driver 595.71.05 / CUDA 13.2 confirmed every
production setting at its current value: TAIL_KERNELS=2, WMUL=2, MULTI_Q=1,
L1CUDA=3, GRAPHS=0, NOREG=0, UNROLL_W/H=1, ZEROHACK_W/H=1, and the current
IN/OUT_WG, PAD, TABMUL, MODM31 choices.  The FFT-shape table also reconfirmed
`1:512:8:512:202` as the only 4M shape covering p=136.28M (cheaper variants
cap at ~100-133M; every 4.7M+ shape is >= 2.4x slower).  The tune appended
one candidate delta: `LOADS=22042->10042, STORES=21->20, TAIL_TRIGS31=1`
(trig-once load 1 instead of 2, trig-several load 0 instead of 2, FFT store 0
instead of 1, tabulated M31 tail trigs).  Probe differences were 0.2-us
class.

Matched alternating 50k verification (all residues exact):

| Pair | Control blocks | Trial blocks | Verdict |
|---|---|---|---|
| 1 | 199.3 / 203.8 / 204.2 / 204.5 | 200.3 / 205.6 / 206.2 / 206.2 | trial +1.7 |
| 2 | 201.6 / 207.0 / 204.5 / 204.9 | 200.2 / 205.1 / 205.5 / 205.9 | trial -0.3 |

Decision: **do not adopt; the tune deltas are thermal-noise-level and do not
reproduce.**  The 2026-08-12 tune otherwise confirmed every production
setting, so the incumbent `-use` line is validated as (still) optimal on the
current driver.  This closes the systematic micro-tuning direction: there is
no latent configuration gain left in the tune space.

### 2026-08-12: closing 1M-iteration production validation

A full fresh 1M-iteration run of the production configuration completed with
the Gerbicz check passing and iteration-1,000,000 residue
`52b03a7cc55e677d`, exactly matching the known-correct value in the
registry.  100k-block averages: 202.5, 204.9, 206.3, 207.6, 208.1, 208.9,
209.3, 209.5, 209.5, 209.8 — a ~207-us whole-run mean with a 209.5-us
thermal steady state, consistent with the 205.5-us historical reference.
During steady state the GPU sits at exactly **300.02 W (the enforced cap)
with SM clocks throttled to 1552 MHz — half the 3090-MHz maximum — at 81 C
and 100% utilization.**  This quantifies the campaign boundary: the software
runs the silicon at its power wall; the ~2x clock headroom of this GPU is
locked behind the Max-Q 300-W limit (325 W permission-blocked in this
container), and no software change measured by either campaign moves the
wall.  The same workload on a full-power (600-W) GB202 board would be the
single largest available speedup.

## Campaign status and boundary (2026-08-12)

No exact end-to-end speedup over the tuned production `M31+M61` path has been
achieved by either campaign.  The Fable campaign's contributions:

1. **Completed and closed Sol's last open lead.**  The detached M31 width
   edge — the only rejected-at-180 idea that still promised a real speedup —
   is now implemented exactly (correct residues through 50k in both modes)
   and measured: mode 1 +20.5 us, mode 2 +11 us.  Nothing in the scaffold's
   apparent 19-us saving survives exact completion.
2. **Corrected the performance model.**  The M31 bottom half is not hidden;
   both queues are work-saturated (M31 kernels dilate 2-3x under
   concurrency).  The iteration is bound by total executed work plus the
   serial fused carry, and relocation between queues is paid at full price.
   This single model explains, post hoc, the failure pattern of essentially
   every overlap/sidecar/batching experiment in the registry.
3. **Closed the micro-tuning space.**  The config.txt "suggested" toggles, a
   fresh full `-tune ntt` on the current driver, and its one candidate delta
   were all verified at noise level or worse in matched long runs.  The
   production `-use` line is confirmed optimal.

Remaining routes to a real speedup, in order of credibility (all outside
what this machine/software state can express):
- Hardware with a higher power ceiling or more SMs (the workload is enforced
  power/issue limited at 300 W; 325 W is permission-blocked here).
- An algebraic reduction in total exact work per iteration — the reopen
  conditions in Sol's file stand (<=8-product M61 pair square, sub-15-us
  ownership-complete odd DFT, Tensor formulation with multiple recoverable
  modular products per accumulator).  All known families are measured out.
- Nsight Compute counters (ERR_NVGPUCTRPERM here) to find an invisible
  stall inside the M61 kernels, if a future host permits them.

Do not reopen without new evidence: width detachment in any direction or
queue arrangement; sidecar transforms of any length; queue-priority/PDL/
graph scheduling; persisting-L2 or full root tables; the tune option space.

## 600 W phase (2026-08-12): power-limit scaling on the full-power board

The campaign migrated (via `setup.sh`) to a Vast.ai box with an RTX PRO 6000
Blackwell **Workstation Edition**: the same GB202 silicon, SM count, and
97.9 GiB memory as the Max-Q board, but a **600 W power limit** (range
150-600 W, default 600).  Driver 595.84, CUDA 13.2.  This realizes the
300 W campaign's closing statement that a full-power GB202 board is the
single largest available speedup, and turns the previously "unavailable
control experiment" (power-limit scaling) into a measurable one.  Still
blocked in this container: `nvidia-smi -pl` and `-lgc` (Insufficient
Permissions; the limit is pinned at 600 W, so no direct sweep), and NCU
counters (ERR_NVGPUCTRPERM), so the invisible-stall route stays closed.

### New production baseline

A 5M-iteration production run (`.pl600-steady.ul7Cke`) passed the Gerbicz
check with exact residues at every checkpoint (1M = `52b03a7cc55e677d`,
matching the registry; new 5M reference `7eca65291732df02`):

- Cool (first ~100k): 137-138 us/iteration at ~2380-2470 MHz.
- Thermal steady state (reached ~1.2M, held through 5M): **148.9
  us/iteration at 2197-2247 MHz (mean 2210), 83-84 C, 600.0 W pinned**,
  active throttle reason = SW Power Cap only — never thermal, never the
  clock ceiling.
- Versus the 300 W steady state (209.5 us at 1552 MHz, 81 C): **1.407x
  sustained throughput from 2.00x power**, i.e. perf ∝ P^0.49 over this
  span.  Energy per iteration rose from 62.9 to 89.3 mJ (+42%).

### Scaling model (task: how does PRPLL scale with the power limit?)

With the limit unsweepable, the model uses the 300 W box as one anchor and
the 600 W thermal ramp as a natural clock sweep at constant power (the boost
governor walks 2467 -> 2200 MHz as the GPU heats at exactly 600 W).  A
single two-parameter model

```text
t(f_SM) = 6.0 us + 315,842 us*MHz / f_SM
```

fits the 300 W steady point, the 600 W steady point, and every 20k-block
(us/iter, mean clock) pair along the ramp — a 1552-2467 MHz span — within
about 1.5 us.  Consequences:

- The iteration is ~96% pure SM-clock-scaled; the clock-independent
  residue is only ~6 us.  Memory-controller activity is 18%; DRAM clock is
  immaterial.  Performance therefore scales with the power limit exactly as
  the power-capped clock does: measured perf ∝ P^0.49 between 300 and 600 W
  (strongly sublinear; V/f curve).
- At the 3090 MHz ceiling the model predicts ~108 us/iteration, but even
  600 W sustains only 2210 MHz: the workload is **still power-limited at
  the board maximum**, with ~29% clock headroom locked behind a power draw
  the board cannot deliver (extrapolating P ∝ f*V^2, roughly 900+ W).
- The `-time` profile at 600 W steady state shows every kernel scaling by
  ~the clock ratio (kCarryFused 71.5 -> 49.7, kfftMidOutGF31 50.8 -> 34.6,
  ktailSquareGF31 39.9 -> 28.4, kfftMidInGF31 17.9 -> 12.0; ratios
  1.40-1.49x vs clock ratio 1.42x).  No kernel became memory-bound, and the
  M31 dilation under concurrency is unchanged (~1.6x).  **The saturation
  model transfers to 600 W intact**: both queues remain work-saturated,
  only deleting work shortens the iteration, and every relocation/sidecar
  closure in the registry stands.

### Does 600 W open new optimization opportunities?

Three candidate reopenings were tested; none flips.

1. **Tune space (rerun per setup.sh).**  A fresh `-tune ntt` reconfirmed
   every structural production setting (TAIL_KERNELS=2, WMUL=2, MULTI_Q=1,
   L1CUDA=3, GRAPHS=0, NOREG=0, TABMUL_CHAIN32=1) and reconfirmed
   `1:512:8:512:202` as the only shape at its speed class covering p136
   (the faster table entries top out at exponent <= 133M; the next covering
   shape costs 168 us).  Its probe-level deltas (LOADS 22042->20042,
   STORES 21->24, MODM31 2->1; plus the commented
   TAIL_TRIGS32=0/UNROLL/ZEROHACK toggles) were put through a three-arm
   alternating 9x100k verification with exact residues at every 20k
   checkpoint: after thermal-drift correction both trial arms sit within
   +-0.5 us of control with inconsistent sign (round 3, near-steady:
   control 145.9, tune line 146.5, toggles 146.2).  **Do not adopt; the
   incumbent -use line remains optimal at 600 W.**  Notably GRAPHS=0 also
   survives the shorter-kernel regime.
2. **Two workers (reopened because the hardware condition materially
   changed).**  `-prps 136279841,136279879 -workers 2`, 200k/exponent,
   thermally matched, residues exact:

   | Config | per-exp us/iter | SM clock | aggregate it/s | vs single |
   |---|---:|---:|---:|---:|
   | 1 worker | 148.9 | 2210 | 6,716 | baseline |
   | 2 workers, TAIL_KERNELS=2 | 323.9 | 1848 | 6,175 | -8.1% |
   | 2 workers, TAIL_KERNELS=3 | 327.5 | ~1848 | 6,107 | -9.1% |

   The loss halves versus 300 W (-17.9%/-12.7%), and the decomposition is
   informative: doubling the resident work now costs only x0.840 in clock
   (x0.786 at 300 W — the V/f curve is steeper up here), while overlap
   efficiency *improves* 9.4% (the second exponent fills queue gaps and the
   serial fused carry the first cannot use).  But the clock penalty still
   dominates; break-even needs f_2w/f_1w >= 0.914.  The TAIL_KERNELS
   preference also flips (2 beats 3 with two workers at 600 W).  **Two
   workers remain rejected at 600 W**; on a still-higher-ceiling board the
   sign plausibly flips — first hardware condition to re-check on any
   future >600 W machine.
3. **Kernel balance / memory wall.**  Higher clocks with fixed DRAM
   bandwidth could have exposed a memory-bound kernel worth re-tuning or
   restructuring around; the `-time` profile above shows this did not
   happen (all kernels ~clock-scaled, memory activity 18%).  No reopening.

### Updated campaign boundary (600 W)

The full-power board delivers +40.8% sustained throughput (209.5 -> 148.9
us/iteration) purely from hardware; no software change is implicated, and
every software conclusion from both campaigns carries over unchanged.  The
GPU still runs at its power wall (SW Power Cap active at 600.0 W, 2210 of
3090 MHz).  Remaining routes, in order of credibility:

- Hardware again: more SMs or a higher ceiling (perf ∝ P^0.49 here; a
  hypothetical uncapped clock ceiling is worth a further ~27%, but needs
  ~900+ W).  On any >600 W or multi-GPU machine, re-check two workers
  first (-8.1% margin at 600 W, improving with ceiling).
- An algebraic reduction in total exact work per iteration.  Of Sol's three
  reopen conditions, two are now investigated to closure (see
  "Standing-gate investigations" below): the <=8-product pair square is
  unreachable (9 proven practically optimal by structured search) and the
  odd-DFT edge gates are unchanged by 600 W and remain 4-30x out of reach.
  The remaining condition (a Tensor formulation with multiple recoverable
  modular products per accumulator) has no supporting instruction on sm_120.
- NCU counters on a host that permits them (still ERR_NVGPUCTRPERM here).

## Standing-gate investigations (2026-08-12, 600 W phase)

The two credible remaining algebraic work-removal gates from Sol's registry
were investigated to completion at the user's direction.

### Gate: sub-15-us ownership-complete odd DFT — closed at 600 W by analysis

Every odd-radix candidate's decisive numbers from Sol's measured record
(300 W us):

| Candidate | Planes | Best-case saving | Odd-edge cost | Verdict |
|---|---|---:|---|---|
| 3M radix-3 (M31+M61+M19) | 3 | core -2.0 vs generic-q | free edge still fails complete budget | closed |
| 3.5M radix-7 GT | 3 | — | dependency-complete core 162.6 vs 129.6 (+25.4%) | closed before carry |
| 3.875M radix-31 | 2 | <=30.6 (generous stage bound) | measured 245.5-246.2; needs sub-8 | >=30x gap |
| 3.9375M radix-63 | 2 | <=38.5 (impossible bound) | measured 142.2-142.6 | 3.7x over the bound |
| 4.0625M radix-65 | 2 | 11.3-12.5 core win | zero-cost edge already insufficient | closed at zero |
| 33/32 near-M61 radix-33 | 2 | 4.6 total allowance | 48-60 estimated | >=10x gap |

600 W moves none of these: every odd-edge kernel is core-clock-domain work
(the radix-31 M61 edge is 502 IMAD + 440 SHFL + 208 warp syncs), and the
600 W profile showed uniform 1.40-1.49x clock scaling across all production
kernels including the shuffle/LDS-heavy middles, with no memory wall.  All
budget ratios are scale-invariant; the sub-8/sub-15-us edge allowances become
sub-5.7/sub-10.7 us at 600 W clocks while the measured edges shrink
identically.  The one named rescue (a GPU-mapped 31-point WFTA) has no
practical construction (flow-graph literature stops at 19; the addition
network grows super-linearly) and would not reduce the measured bottleneck,
which is ownership/shuffle traffic, not multiplication count.

Decision: **odd-radix routes stay closed on this hardware; reopen only per
Sol's original edge-demonstration conditions (scaled to 600 W clocks).**

### Gate: <=8-product M61 tail pair-square — closed as unreachable (9 is optimal)

Registry re-check before work: the "scaled two-multiply rotation audit"
leaves exactly this gate open; no prior audit had exploited the norm-one
relation (all tail trigs lie in the order-2^61 norm-one subgroup since
M61+1 = 2^61, so conj(t) = 1/t and t0^2+t1^2 = 1) as a search dimension.

Method ([`tools/pairsq-search/`](tools/pairsq-search/)): by Strassen
homogenization, any straight-line program for the four homogeneous quadratic
targets (adds/shifts/conjugations free, every wide product counted, arbitrary
precomputed T-constants free) reduces with no more products to k formations
(kappa x linear) + m bilinears + j scalings (kappa x quadratic), chained,
with T entering only through the kappas.  The continuous relaxation of every
shape (free real combo coefficients shared across 12 simultaneous norm-one
T samples, kappas free per T) was optimized with batched-restart Adam +
LBFGS polish.  A methodological trap worth recording: at 4 T samples the
model happily *interpolates* the samples (an apparent (2,6,0) "8-product
hit" failed the frozen-structure fresh-T test); 12 samples kill
interpolators.

Results:
- Flat bilinear rank of the system is **4** (with T-dependent operand
  coefficients free): -T is always a QR in GF(M61^2), so the pair-square
  splits into two independent complex squares of z+- = a +- i*t*b.  The
  honest cost of that route is 3 (form t*b) + 4 (squares) + 3 (unscale d
  by conj(t)) = 10 — the entire gap between 4 and 10 is T-scale
  manufacture, which is what the gate is really about.
- The c-subproblem alone (a^2 - T*b^2) costs exactly **6** (every
  5-product shape fails; (3,3,0) reaches an exact zero); d alone is the
  classical 3.  An 8-product joint scheme therefore requires sharing one
  product across the two subproblems.
- **No such sharing exists**: all 28 shapes at total 8 fail (best residual
  2.1e-4, reproducible across 3 independent T-seeds at heavy settings),
  while the known 9-product schemes are rediscovered to residual 1e-30 in
  every seed — including an independent rediscovery of the exact
  (a+b)^2 - a^2 - b^2 spelling already present in the kernel source.
  The norm-one relation does not rescue 8 (the whole search ran at
  norm-one T).

Decision: **close the <=8-product gate as unreachable: 9 wide products is
optimal for the canonical pair-square** (strong numerical evidence —
continuous infeasibility, not a formal proof).  Both 9-product spellings
are already implemented (`ENABLE_BETTER_ONEPAIRSQ`,
`ENABLE_PARALLEL_ONEPAIRSQ`) and measured tied-or-lost against the
production 10 at 300 W; the production tail pair-square is at its practical
optimum.  Do not reopen with basis changes, scaled rotations, or norm-one
identities — all are inside the searched model.  A future reopening needs
either a formal disproof of this search's negative (an exact 8-scheme it
somehow missed) or a representation change that removes the canonical
carry boundary (Sol's original alternative condition, unchanged).

### Literature check on the pair-square gate (2026-08-12, web)

The 8-vs-9 question was checked against published art after the search
closure.  The operation is exactly **squaring in a quadratic extension**
`F_q[j]/(j^2 - beta)` with `q = M61^2`, `beta = -T` (equivalently, doubling
in the Pell-conic/Brahmagupta group with `D = -T`) — a shape pairing-based
cryptography has optimized for two decades with beta-multiplications charged
as a first-class cost, i.e. the same cost model as the gate.

- [Devegili–O hEigeartaigh–Scott–Dahab, "Multiplication and Squaring on
  Pairing-Friendly Fields" (eprint 2006/471)](https://eprint.iacr.org/2006/471)
  catalogs the complete known inventory: schoolbook `M+2S+B` (= our
  canonical 10 after expanding over the GF(M61^2) tower), Karatsuba-adapted
  squaring `3S+B` (= our 9, literally the `(a+b)^2-a^2-b^2` spelling in
  `tailsquare.cl`), and complex squaring `2S` **only when beta is cheap**
  (beta = -1; unavailable for the dynamic tail trig).  Nothing below `3S+B`
  for generic charged beta appears there or in the follow-on literature
  (Chung–Hasan asymmetric squaring is cubic-extension material).
- Classical complexity theory —
  [Ja'Ja', "Optimal Evaluation of Pairs of Bilinear Forms"](https://epubs.siam.org/doi/10.1137/0208037),
  Winograd, de Groote,
  [Alder–Strassen algebra bounds](https://link.springer.com/chapter/10.1007/978-3-662-03338-8_17) —
  works in the scalar-multiplications-free model: it confirms our measured
  flat rank 4 (the split-algebra squares) but is structurally silent on the
  charged-constants question that separates 8 from 9.
  [Heideman–Burrus](https://link.springer.com/book/10.1007/978-1-4612-3912-3)
  is the charged-constants framework, but its results cover DFTs and
  convolutions, not this pointwise system.
- Mersenne-community sources (mersenneforum PRPLL/M61 threads) contain no
  deeper treatment of the pair-square multiplication count.

Conclusion: **no published formula reaches 8, and no published theorem
forbids it; the search's 9-optimality result is consistent with the state of
the art and appears to be a small novel result.**  The gate closure stands
unchanged; the community-standard best (`3S+B`) is what production already
ties against.

## Hardware-selection notes at the 600 W boundary (2026-08-12)

Assessment of candidate hardware for the next throughput step.  All
non-GB202 numbers are **estimates** from public specs scaled through the
measured workload profile (integer-IMAD issue-bound, ~96% SM-clock-scaled,
18% memory-controller activity, 48-MiB state); GB202 numbers are measured.

| Device | Path | us/iter | Power | Note |
|---|---|---:|---:|---|
| RTX PRO 6000 (GB202, 188 SM) | M31+M61 integer | **148.9 measured** | 600 W | this box |
| RTX 5090 (GB202, 170 SM, ~$2k) | M31+M61 integer | ~160-165 est. | 575 W | ~4x better it/s/$ than any alternative |
| A100 SXM | M31+M61 integer | ~700-900 est. | 400 W | 6,912 INT lanes x 1.41 GHz = 5.4x lane-clock deficit |
| A100 SXM | classic FP64 FFT | ~300-450 est. | 400 W | 9.7 TF FP64, 2 TB/s, 40 MiB L2 |
| H100 SXM | best of either | ~150-250 est. | 700 W | parity at best, ~10x the price |
| B200 | FP64/HBM path | ~60-100 est. | ~1 kW | only per-device candidate to beat GB202; hopeless per dollar/watt |
| MI300X (wildcard) | OpenCL FP64 path | ~120-180 est. | 750 W | 81.7 TF FP64, 5.3 TB/s; untuned path, high uncertainty |

- **A100 does not beat this box on either of its paths.**  Its strengths
  (FP64 tensor, NVLink, HBM) are all things this workload measured out or
  does not use.  M52's own discovery ran the FP64 backend on A100 — that
  fleet's economics were donated compute, not efficiency.
- **The per-dollar optimum for this codebase is a power-capped RTX 5090
  fleet**: same silicon, and the measured perf ∝ P^0.49 makes several
  capped cards strictly better than fewer full-power cards at equal wall
  power (two 300 W GB202 ≈ +42% over one 600 W).
- **The highest-value single box change is administrative, not silicon**: a
  host with working `nvidia-smi -pl`/`-lgc` and NCU counters unlocks the
  perf/W knee, the counters route, and the two-worker break-even map.
- On any >600 W or multi-GPU machine, re-check two workers first.

### FP64 tensor cores — rejected as an idea class for any current hardware

Raised as a candidate optimization; rejected on analysis without
implementation.  (1) sm_120 has **no FP64 MMA units at all** and 1/64-rate
scalar FP64 — nothing to optimize here (registry already rejected the FP64
backend and FP64-in-M61 offload at the rate gate).  (2) On A100/H100, which
do have DMMA: casting radix stages as dense matmuls inflates real flops
~10x for radix-8 (the DFT matrix is mostly +-1/+-i entries that butterflies
turn into additions), while the FP64 tensor/vector throughput ratio is only
**2x** (A100 19.5/9.7, H100 67/34 TF) — a ~5x net loss wherever
compute-bound, and no gain where memory-bound.  This is the same failure
mode the registry *measured* on dense INT8 tiles (42% slower despite a much
larger nominal TOPS ratio); published tensor-FFT successes (tcFFT etc.)
live in FP16/TF32 where the ratio is 8-16x.  (3) The exact-NTT variant
(integer limbs in FP64 mantissas via DMMA) starts below the GB202 integer
path's 5.3e13 IMAD/s before paying ~6x limb decomposition.  (4) Roadmap:
B200 cuts FP64 tensor to ~vector rate (~40 TF) and consumer Blackwell has
none — the unit is being removed, not grown.  FP64 MMA also fails the
standing tensor reopen condition (multiple recoverable modular products per
accumulator — it delivers one).  Do not revisit absent a part with a >=10x
FP64 tensor/vector ratio.

## Admin-host phase (2026-08-12): hardware controls unlocked

The campaign migrated (setup.sh) to a bare-metal-class host with an RTX PRO
6000 Blackwell **Server Edition**: the same GB202 silicon and 188 SMs, but
max SM clock 2430 MHz (vs the Workstation's 3090), memory 12481 MHz with
**ECC enabled**, power range 300-600 W (default 600), driver 580.126.09,
CUDA 13.0, passwordless sudo.  Baseline validation exact (2k/100k residues
match; ~150.8 us cool).  This host unlocks all three items the 600 W
boundary left open:

- `sudo nvidia-smi -pl` **works** (300-600 W, verified by setting).
- `sudo nvidia-smi -lgc/-rgc` **works** (~120 graphics steps to 2430).
- **NCU counters work** (`RmProfilingAdminOnly=1` + sudo; full metric sets
  captured on production kernels; residues stay exact under profiling).
- `-lmc` exists but the board exposes only two memory states (12481/405
  MHz): no usable memory-power knob.
- Thermal drift is absent (69 C max at 600 W): matched measurements no
  longer need drift correction; a 600 W point repeated after a 30-minute
  sweep reproduced within 0.2 us.

### Direct perf/W sweep (efficiency knee) — MEASURED

Fresh 1M-iteration production runs per power limit, all with exact registry
residues at 1M (`52b03a7cc55e677d`); tail = blocks >= 600k:

| W | us/iter | SM MHz | it/s | mJ/iter | it/s/W | local d ln perf/d ln P |
|---:|---:|---:|---:|---:|---:|---:|
| 300 | 215.5 | 1641 | 4641 | 64.7 | 15.5 | - |
| 350 | 198.9 | 1763 | 5027 | 69.7 | 14.4 | 0.52 |
| 400 | 182.6 | 1961 | 5476 | 72.8 | 13.7 | 0.64 |
| 450 | 170.7 | 2103 | 5858 | 76.6 | 13.1 | 0.57 |
| 500 | 163.4 | 2196 | 6120 | 81.7 | 12.3 | 0.42 |
| 550 | 157.5 | 2274 | 6351 | 86.5 | 11.6 | 0.39 |
| 600 | 153.3 | 2338 | 6522 | 91.6 | 10.9 | 0.31 |
| 600rep | 153.1 | 2339 | 6533 | 91.9 | 10.9 | (drift ctrl) |

- **The knee is at ~450-500 W**: the scaling exponent holds ~0.55-0.6 up to
  450 W then collapses to 0.31 by 600 W; the last 100 W buy +6% throughput.
  The Workstation two-point estimate (P^0.49) was the average of this curve.
- **Perf/W rises monotonically to the 300 W floor** (+42% it/s/W vs 600 W).
  Fleet rule measured, not estimated: two 300 W GB202 = 9,282 it/s vs one
  600 W = 6,522.  For it/s/$ on owned hardware, cap at the board minimum;
  for single-card latency, run at 600 W and accept the 0.31 exponent.
- **The Server Edition runs ~8.4% more cycles per iteration than the
  Workstation at every clock**: fitting t = c + k/f to the sweep gives
  c ~ 7 us, k ~ 342,000 us*MHz vs the Workstation's 6.0 + 315,842/f.  A
  multiplicative (per-cycle) penalty is the signature of the ECC-enabled
  memory path stretching latency-exposed cycles (see NCU findings), making
  ECC-off the top hardware experiment on this host.

### NCU counters route — the "invisible stall" is found and named

Full-set NCU capture (42 launches of the 7 steady-loop kernels, isolated,
~2.2 GHz, ECC on).  Per-kernel medians:

| kernel | us | SM% | DRAM% | issue% | CPI | occ theo/ach | regs | limiter |
|---|---:|---:|---:|---:|---:|---|---:|---|
| carryFused | 77.6 | 52.5 | 45.3 | 38.6 | 10.7 | 41.7/36.9 | 96 | regs (5 blk) |
| tailSquareGF61 | 46.9 | 53.1 | 50.6 | 38.4 | 10.9 | 41.7/37.6 | 96 | regs (5 blk) |
| fftMiddleInGF61 | 34.8 | 26.5 | 61.0 | 17.6 | 34.9 | 66.7/55.9 | 64 | regs (4 blk) |
| fftMiddleOutGF61 | 32.9 | 28.0 | 64.1 | 18.3 | 33.8 | 66.7/55.9 | 64 | regs (4 blk) |
| tailSquareGF31 | 27.2 | 58.6 | 43.6 | 41.9 | 11.2 | 50.0/42.7 | 60 | shmem (6 blk) |
| fftMiddleInGF31 | 16.5 | 29.8 | 64.0 | 20.9 | 33.1 | 83.3/69.6 | 48 | regs (5 blk) |
| fftMiddleOutGF31 | 16.7 | 29.9 | 63.8 | 21.0 | 33.8 | 100/73.1 | 40 | lg_throttle |

Findings (per-instruction stall sampling):

1. **No kernel is pipe-saturated in isolation.**  The hottest pipe anywhere
   is ALU-heavy at 52-59%; issue slots are 17-42% busy.
2. **The dominant stall everywhere is `long_scoreboard`** — unhidden
   global-memory latency landing on the first `IADD.64`/`IMAD.WIDE.U32`
   consumers of loaded residue words (33% of carryFused stall samples, 61%
   of fftMiddleOutGF61's).  Not an exotic pipeline hazard: plain exposed
   DRAM/L2 latency with too few resident warps to hide it.
3. **Occupancy is compile-time register-capped**: carryFused/tailSquareGF61
   at 96 regs -> 5 blocks of 128 threads -> 20 warps/SM; GF61 middles at 64
   regs -> 32 warps.  The defaults were tuned upstream on 4090/5070Ti
   (comments in `Gpu.cpp:numCudaRegisters`), without ECC, and are
   overridable per kernel via `-use REGCF3161/REGTS61/REGMI61/REGMO61/...`.
4. The GF31 middle-out is the exception: 100% theoretical occupancy with
   `lg_throttle` 10.4 (LSU queue full) — already latency-limited the other
   way; reg caps cannot help it.
5. If latency were fully hidden, the ALU-heavy pipe bound puts carryFused's
   floor at ~45-50 us at 2.2 GHz (vs 77.6 measured) — an upper bound on the
   occupancy prize, shaved in production by the power wall (higher
   issue/cycle -> more W/cycle -> lower f at the cap).
6. This mechanism also retro-explains the multiplicative Server/ECC
   penalty: ECC adds memory latency; latency-exposed cycles scale with
   core clock; hence k grows, not c.

### Directed occupancy experiment (register caps) — CLOSED, negative

The one-step register reductions the crossover table allows (96->80 on
carryFused/tailSquareGF61 for 20->24 warps/SM; 64->48 on the GF61 middles
for 32->40; a prefer-shared+64-reg arm for 8-block carry/tail) were run as
matched alternating fresh 100k runs at 600 W, 3 rounds, all residues exact,
controls 152.45 +- 0.5 us:

| arm | -use delta | d vs control |
|---|---|---:|
| R1 | REGCF3161=80 | +1.70 us |
| R2 | R1 + REGTS61=80 | +2.79 us |
| R3 | REGMI61=48,REGMO61=48 | +31.3 us (spills) |
| R4 | R2 + R3 | +33.1 us |
| R5 | L1CUDA=1,REGCF3161=64,REGTS61=64 | +20.6 us |

**Occupancy-via-register-caps is closed on GB202**: the extra
instructions/spills from tighter budgets always exceed the latency-hiding
gain, confirming upstream's 4090 tuning.  The intermediate steps do not
exist (88 regs buys no block at 128 threads; 56 buys none at 256), so the
space is exhausted, not sampled.  Corollary: the latency exposure is
structural at current occupancy — the remaining attacks on the
long_scoreboard stall are (i) reduce the latency itself (ECC-off; measured
next) and (ii) software prefetch/async-copy restructuring of the middles
(the registry's "TMA rewrite — no evidence" closure now has direct counter
evidence on its premise; see the Sol-registry audit).

### Sol-registry methodology audit (2026-08-12)

Six parallel auditors re-read all 9,578 lines of Sol's registry against the
post-Sol facts (F1 corrected concurrency model, F2 NCU latency findings,
F3 measured power scaling, F4 ~2-us noise floor, F5 scaffold-bound lesson),
hunting false negatives.  Verdict: **the registry's closures are largely
sound** — most rest on exact end-to-end measurements with 20-100+ us
margins or algebraic impossibilities, and Sol's late-era protocol
(alternated, order-reversed pairs) was good.  But the audit surfaced one
systematic pricing error, two class-level foreclosures made on premises now
overturned, and a handful of noise-level closures.  Cross-checked reopen
shortlist:

1. **TMA/cp.async prefetch of the middle kernels — REOPENED (HIGH).**
   Sol closed it analytically: "Without counters there is no evidence for
   a ... gain from merely replacing its existing cooperative loads"
   (counters were permission-blocked).  The counters now exist and show
   the middles are exactly what async prefetch targets: CPI 33-35,
   long_scoreboard 22-24 cycles/instruction, DRAM 61-64%, no pipe above
   30%.  The closure's stated reason is void.  Effort: days (kernel
   rewrite); the diagnostic half is already done.
2. **q24-replaces-M31 field swap — re-test the unmeasured charge (HIGH-MED).**
   The strongest measured M31-shrinking candidate (co-run advantage
   9.5-10.0 us at chain 8, 18.2 us at chain 16) was dismissed by an
   UNMEASURED analytical charge (generic roots/weights/CRT "eat the rest")
   against a 25.5-us denominator inflated by baseline drift.  Under F1,
   M31-side work reduction is paid back at the dilated co-run rate, not
   zero.  First gate: extend `q24_m61_overlap_bench.cu` with generic
   twiddle products + CF weights, matched alternation (hours-day).
3. **M31 co-run price: derived, probe deprioritized.**  The "hidden M31"
   decision rule priced all M31-side savings at zero; Sol's own data
   contradicts it, and the payback coefficient can be DERIVED from the
   registry without a new run: the detached-edge scaffold bought 18.8 us
   end-to-end for ~22-26 us of isolated M31 width work deleted, and the
   dilation data (65.6 us isolated -> ~108 us co-run) brackets the same
   quantity.  **Price M31-side deletions at ~0.7-0.9x payback.**  The
   correctness-off thinning probe adds little and is parked at lowest
   priority (user call, 2026-08-12).
4. **TAIL_TRIGS61 generate-vs-load (cheap, new).**  Cross-checking the
   audit's "twiddle generation no-go overgeneralization" flag against the
   code: the GF61/GF31 tails default to READING all trig values from
   memory (`TAIL_TRIGS61=0`) — the memory-heavy choice inside the
   latency-exposed tail — and no tune record toggles the GF61 knob.  The
   opposite direction (full middle-root TABLES) lost 0.7-0.9% in
   production, which under F2 pricing is evidence FOR generation.  Minutes
   to test.
5. **Resident 64-KiB M61 tile under real contention (MED).**  The isolated
   proxy that rejected it (+1.9-2.7%) was fully L2-resident, erasing
   exactly the global round trips the fusion deletes; under production
   L2/DRAM contention the sign could flip.  Re-run the existing bench
   beside an M31-shaped memory load (hours).
6. **Two-field radix-7 q7/M61 3.5M (MED).**  Rejected on a ~1.5-us wash
   between two DIFFERENT benchmark types, with the candidate's 12.5%
   smaller state/carry traffic priced at zero — the currency F2 says is
   binding.  Optimistic repricing brushes the old 25.5-us requirement.
   Re-run both existing benches alternated + NCU the edge kernel (hours).
7. **Two-61-bit-field 3M (MED-LOW).**  Rejected at ALU-chain depth 8-16
   (loses 10-13%) but ties/wins at chain 2 — the memory-side regime F2
   says production actually occupies; the candidate moves 25% fewer bytes.
   Needs a memory-realistic tile gate (1-2 days).
8. **Near-M61 radix-33 (LOW-MED)** — same invalid ALU-chain pricing, but
   corrected pricing must net under 4.6 us; likely still rejects (0.5-1 d).
9. **Batch-native M61 AoSoA (condition now resolved: stays closed here).**
   Robust for the 180-us gate, but its aggregate-throughput question was
   never closed by Sol.  The two-worker map (below) answers the condition:
   aggregate deltas are -11.3/-7.7/-5.3/-0.1% at 300/400/500/600 W —
   parity at this board's ceiling, never a gain.  Reopen batch-native (and
   two workers) only on a >600 W or multi-GPU host, where the monotonic
   trend implies a positive sign.

Flags that DISSOLVE on cross-check: the 9-vs-10-product pair-square
closures (noise-level as recorded, but re-closed decisively by this
campaign's structured search: 9 is optimal, both spellings tied-or-lost);
middle-root generation for the middles (production already generates —
the "no-go" never governed production); the middle-6 geometry scaffolds
(moot via the parent 3M closure); stream-priority and PDL (F1 supports the
closures); Sol's noise-level Harvey/lazy-q deltas (route independently
dead via the M19 bound).

Model corrections for the ledger: (i) F1's "deleting side-queue work pays
back ~fully" must NOT be generalized to mixed-pipe engines — the FP32
factor-four data shows ~76% absorption there; (ii) several robust closures
lean on one sub-noise number (the 1.28-us radix-9 core delta feeding the
PFA33 bound) — their margins survive 3x error, so no action; (iii) the 3M
hybrid scaffold (185.6 us "budget") is confirmed as the F5 exemplar: a
false POSITIVE that consumed the campaign's largest wasted effort.

### Locked-clock t(f) ladder (-lgc) — this box's model

100k production runs at locked SM clocks, 600 W limit, all residues exact:

| lock | meas MHz | us/iter | W drawn |
|---:|---:|---:|---:|
| 1500 | 1492 | 237.97 | 275 |
| 1700 | 1695 | 209.47 | 311 |
| 1900 | 1879 | 191.20 | 367 |
| 2100 | 2085 | 171.75 | 428 |
| 2300 | 2272 | 157.53 | 534 |

Fit: **t = 4.98 us + 347,700 us*MHz / f** (max residual 1.2 us; predicts
the unlocked power-sweep points within ~1.4 us).  Versus the Workstation's
6.0 + 315,842/f: **k is +10.1%** — the refined estimate of the Server
Edition's per-cycle handicap (two-point sweep estimate was +8.4%), now the
quantitative target for the ECC-off experiment (~15 us at 600 W if ECC
explains all of it; driver 580-vs-595 is the confound the on-box ECC
toggle isolates).  The locked P(f) points (275->534 W over 1.52x clock)
also give the V/f curve directly.

### ECC-off — measured; the ECC-tax hypothesis is FALSIFIED

Applied `-e 0` (in-place `--gpu-reset`, driver reload, and PCI
remove+rescan all failed to apply it; only a cold reboot works on this
board).  All residues exact (2k/…/100k and the 1M `52b03a7cc55e677d`).

| point | ECC on | ECC off | delta |
|---|---:|---:|---:|
| locked-2100 (2085 MHz), 100k | 171.75 us | 171.7 us | **0.0** |
| unlocked 600 W, 100k tail | 152.45 +- 0.5 (n=15) | 150.2-151.1 (n=3) | -1.8 us |
| unlocked 600 W, 1M tail | 153.3 us @ 2338 MHz | ~152 us @ 2370 MHz | ~-1.3 us |

- **At fixed clock ECC costs nothing**: the inline-ECC read path adds no
  visible latency to this workload.  The unlocked ~1.2% gain is purely the
  power channel — ECC logic/traffic draws board power, and removing it
  buys ~32 MHz at the 600 W wall.
- **The +10.1% k gap vs the Workstation box is therefore NOT ECC.**
  Remaining suspects: driver 580.126 vs 595.84, board firmware/memory
  timings, or bin.  Not actionable on this host; the next box with driver
  595+ should re-fit t(f) to isolate the driver term (setup.sh now
  snapshots `nvidia-smi -q`, so ECC state will be recorded).
- Decision: **keep ECC off on campaign boxes** (+1.2% free; Gerbicz + PRP
  double-checking cover integrity), with expectations calibrated: it is a
  power optimization, not a latency one.

### Driver 580 -> 595.84 experiment: setup and capture (in progress)

The +10.1% k gap's last suspects are driver and board.  This box's apt
offers exactly **595.84** — the Workstation box's driver — making a
controlled swap possible.  Mechanism confirmed before the swap: PRPLL
compiles OpenCL sources via NVRTC to PTX and loads PTX with
`cuModuleLoadData` (src/cuda/cudawrap.cpp:487,121), so **the driver's JIT
emits the final SASS of every kernel**; a driver swap changes the code
generator, the GSP firmware (clock/power arbitration), and power
management together.  Pre-swap capture under 580.126.09 in
`.driver580-capture/` (local to this box): JIT-cache cubins for all 19
kernels, disassembled SASS + instruction counts (carryFused 4664,
tailSquareGF61 3792, fftMiddleInGF61 1304, ...), NVRTC PTX (control —
toolkit-owned, driver-independent), GSP/VBIOS/package versions, supported
clocks.  Post-swap readout: SASS differs + perf moves = codegen; SASS
same + perf moves = firmware/power; nothing = board/bin (close the gap as
unactionable).

### TMA/cp.async middle-prefetch feasibility — GO (design ready)

Source-verified: the middles have NO latency hiding today (one tile per
thread, all MIDDLE=8 16-B loads batched up-front in `readMiddle*Line`
(middle.cl:845-847), first cmul eats the full round trip — the measured
61% long_scoreboard).  Kernels compile as CUDA C++ via NVRTC (sm_120)
with inline PTX already pervasive (base.cl:399-640), so cp.async staging
needs no new host plumbing; full TMA/cuTensorMap is NOT recommended (XOR-
swizzled 256-B chunks, no sm_120 clusters).  Upstream's PREFETCHL1 note
("tried in fftMiddleInGF61 on a 5080 with no benefit", base.cl:758) tested
scalar prefetch without restructuring — structurally unable to help,
doesn't price this design.  MVP: fftMiddleOutGF61 INPLACE branch, grid
1024->256 with a 4-tile block loop, 32-KiB cp.async stage (occupancy falls
4->2 blocks/SM — the main bet), plus a half-staging hedge arm (16 KiB,
keeps 4 blocks/SM).  DRAM-floor ceiling: 32.9 -> ~24 us isolated per
kernel; end-to-end estimate 3-6 us for the MVP, 9-17 us (6-11%) with all
four middles.  Risks and mitigations recorded (REGMO61 override hook for
pressure; spills are the known catastrophic mode; judge only by exact
100k end-to-end runs).  Sequence AFTER the driver decision — codegen
tuning must land on the final JIT.

### Driver 580 -> 595.84: RESOLVED — the gap was the driver's memory path

Swap executed (apt single-transaction 580-open -> 595-open, reboot; ECC
stayed off; revert path `nvidia-driver-580-open`).  SASS diff first: 18 of
19 kernels identical instruction counts (fftMiddleOutGF61 bit-identical;
tailSquareGF61 only cosmetic IADD3->IADD spellings; carryFused +16 of
4664, reordered prologue).  **Codegen is NOT the story.**  Perf battery
(all residues exact):

| point | 580/ECC-off | 595.84/ECC-off | delta |
|---|---:|---:|---:|
| unlocked 600 W, 100k tail | 150.2-151.1 | **141.8-143.2** | ~-7.5 us (-5%) |
| locked-2100 (2084 MHz) | 171.7 | **160.9** | -10.8 us (-6.3%) |
| ladder fit k | 347,700 | **325,116** (c 5.16) | -6.5% |

- Identical SASS at identical clock running 6.3% faster = the GSP
  firmware's memory-subsystem configuration (timings/latency), exactly
  where the NCU stall profile said the time goes.  The "+10.1% Server
  tax" was ~2/3 the 580 driver; residual vs the Workstation k is +2.9%
  (board/bin at most).
- **New campaign-best baseline: ~142-143 us/iteration at 600 W** — faster
  than the Workstation box's 148.9.  1M steady confirm + two-worker
  re-check (break-even was -0.1% under the old memory path) running.
- Fleet/ops rule: **driver major version is a first-class performance
  variable on this workload (~5-6%); pin and record it** (setup.sh
  snapshots it); prefer 595.84+ over 580.x everywhere.
- Confirmations: **1M steady 144.67 us** @ 2319 MHz / 598.5 W, 1M residue
  exact — the new production baseline.  Two-worker recheck under 595.84:
  6,842 aggregate it/s vs 6,912 single (-1.0% vs steady, ~-2% vs
  100k-scale) — the better memory path helps 1w slightly more than 2w;
  two workers remain rejected at <=600 W.
- Day summary: three stacked hardware-layer gains on one box (ECC-off
  +1.2%, driver +5-6%) took production 153.3 -> 144.7 us steady after
  ~200 software rejections across three campaigns.  Remaining software
  lead: the cp.async middle-prefetch MVP (re-anchor its NCU stall profile
  under 595.84 first — the firmware change shrank the very stalls it
  targets).

### cp.async middle-prefetch MVP — IMPLEMENTED, exact, and CLOSED negative

Implemented per the GO design (`-use ASYNC_MID61=1/2[,ASYNC_TILES=n]`,
default off): tile-looped fftMiddleOutGF61 with cp.async staging
(base.cl CP_ASYNC16/COMMIT/WAIT primitives; middle.cl
asyncReadMiddleOutLine; per-kernel 100% shared carveout via new
`cudaSetKernelSharedCarveout`).  Both arms residue-exact at every
checkpoint on the first build (the same-thread staging needs no barriers).
Alternating 3-round 100k A/B at 600 W on driver 595.84:

| arm | design | d vs control (145.2 warm) |
|---|---|---:|
| ASYNC_MID61=1 | full 32-KiB stage, 2 blocks/SM | **+4.55 us** |
| ASYNC_MID61=2 | half stage + direct prefix | +20.3 us |
| =1, ASYNC_TILES=8 | 128-block grid | +27 us (SM starvation) |

Mechanism: (i) per-tile direct loads serialize latency exposure (A2);
(ii) the full pipeline works structurally but halves resident blocks, and
in production the M31 co-run stream was ALREADY backfilling the middles'
latency stalls — the isolated NCU profile's idle-issue slack is consumed
by the other queue, so intra-kernel hiding buys nothing while the
occupancy cost is real.  This is the saturation model's "no idle
execution reservoir" expressed at warp granularity, now measured directly.
The one untested permutation (6-of-8 staging at 3 blocks/SM) is bounded
by the same trade and not worth its complexity.  **Code retained as
opt-in evidence; default path untouched.  Lesson for the registry: NCU
isolated-kernel latency slack is NOT exploitable for co-run kernels —
only the serial carryFused's slack is real, and its occupancy is
register-bound (closed).  The counters route is now fully walked.**

### q24-replaces-M31: charged gate reconstructed and MEASURED — rejection confirmed

Audit shortlist item 2 (HIGH-MED) resolved.  The lost overlap bench was
reconstructed as [`src/cuda/q24_m61_charged_bench.cu`](src/cuda/q24_m61_charged_bench.cu)
(committed this time) from the ledger's verbatim reducer proof, mirroring
the surviving near61 harness, with the q65-era alternation protocol and
full-array host oracles for both q24 kernels.  Fidelity anchors: the
short-quotient kernel compiles to 19 registers exactly as Sol recorded;
Sol's free-root advantage reproduces in sign and shape on this box
(q24+M61 beats generic-M31+M61 by ~2.0 us at chains 2-8, 10.7 at 16;
compressed from Sol's 300 W numbers by clocks/driver).  New charged arms:
per-round table-loaded generic qC twiddle product + entry/exit CF weight
products (candidate), and power-of-two-rotation M31 (production-realistic
control).  600 W, driver 595.84, all oracles PASS:

| chain 8 co-run | us |
|---|---:|
| M31cheap+M61 (realistic control) | 82.6 |
| M31gen+M61 (Sol's control) | 97.3 |
| q24+M61 (free roots) | 95.3 |
| **q24charged+M61** | **132.7** |

The root/weight charge MORE THAN DOUBLES the isolated q24 kernel (21.1 ->
47.7 us at chain 8); the charged candidate loses by +50.1 us against the
realistic control and +35.4 us even against Sol's pessimistic one, at
every chain length.  Production M31's cheap rotations are measured at
3.5x cheaper than generic products (8.0 vs 27.7 us isolated).
**Decision: the q24 field swap stays rejected; the audit flag is closed
by measurement, and Sol's unmeasured analytical charge turns out to have
UNDERSTATED the cost.**  Calibration for the remaining shortlist: this
was the strongest flagged candidate, and its measured charge exceeded the
assertion — the resident-tile and radix-7 reconstructions (days-class)
should be priced with that prior.

### Resident M61 tile under contention: sign FLIPS — reopened for production

Audit shortlist item 5 resolved by reconstruction
([`src/cuda/m61_resident_tile_contention_bench.cu`](src/cuda/m61_resident_tile_contention_bench.cu),
committed): 512 exact 8x512-tile 4096-point NTT squares (fused 64-KiB
resident kernel vs three-kernel control, shared device transforms),
gated by direct-convolution validation of the host NTT, full 2M-value
path agreement, and 8-tile host oracles.  q65 alternation, stream-A
makespan under a stream-B M31-shaped co-runner at two intensities:

| arm | fused - ctl3 |
|---|---:|
| isolated (Sol's regime) | +34.7 us (ratio 1.097; Sol's leaner design: 1.019-1.027) |
| moderate load (~0.7 TB/s, production-like DRAM pressure) | **-18.8 us — fused WINS (0.959)** |
| heavy load (~1.4 TB/s, saturating) | +2.9 us (parity; all queued) |

Sol's isolated gate was decided by the proxy's L2 residency: the control's
two extra 64-MiB round trips were free there and are NOT free under
production-like contention.  The ~53-us isolated->moderate swing is
traffic-driven (design-independent to first order); with a Sol-lean fused
kernel the moderate-load win would be larger still.  Notably the co-runner
also crowds SMs, so the fused kernel's 1-block/SM occupancy cost is
already priced in — unlike the cp.async case, deleting traffic (F1's
winning currency) beats the occupancy loss.

**Status: the M61 middle/height resident-fusion route REOPENS for
production consideration** — the first audit flag to survive measurement.
Production translation: fusing fftMiddleInGF61+tailSquareGF61+
fftMiddleOutGF61 deletes ~128 MiB/iteration of round trips during the
co-run window; naive scaling suggests a 10-25 us/iter potential at 600 W.
The hard part stands as Sol recorded: the real stripe/twiddle mapping
under INPLACE swizzling, the pair-square (not plain pointwise), and the
64-KiB/1-block/SM regime — a multi-day kernel project, now with a
measured motivation instead of a measured rejection.

### Radix-7 q7/M61: unified gate reconstructed — rejection CONFIRMED

Audit shortlist item 6 resolved
([`src/cuda/q7_m61_radix7_gate_bench.cu`](src/cuda/q7_m61_radix7_gate_bench.cu),
committed): the core saving and edge increment measured in ONE process
with q65 alternation (Sol subtracted them across two benchmarks), plus
the never-measured co-run hideability arm.  q7 Montgomery constants from
the ledger (the 65-bit REDC carry note reproduced as code); full-array
q7-chain oracle + 4096-group cyclic-7 convolution oracles both fields.
600 W / 595.84 medians: core saving S = +49.1 us (my instruments; Sol's
18.7-20.3 at 300 W), direct-DFT-7 edge E_seq = +107.9, co-run
E_conc = +95.4.  Decisive findings:

1. **Hideability measured at ~12%** (E_seq -> E_conc recovers 12.5 of
   107.9 us): the audit's F2 repricing hoped dense register-resident edge
   ALU would hide in the other stream's idle issue slots — it does not;
   both edge kernels are ALU-saturated and barely overlap.
2. Sensitivity bracket: my direct DFT-7 overcharges ~5-7x vs Sol's Rader
   graph; applying the measured 12% discount to SOL's edge numbers gives
   E ~ 17.9 vs S ~ 18.7-20.3 — a +1-2 us wash at best, before generic q7
   weights, wider CRT/carry, and the tight p150 capacity, all candidate-
   side charges.  In my internally consistent instruments E_conc >> S.
3. The ~12% co-run absorption figure retro-validates the whole odd-radix
   edge closure family (radix-31/63/65 edges are denser still).

**Decision: q7/M61 (and by extension the odd-radix family) stays
rejected; the audit flag is closed by measurement.**

### FUSED31: production M31 middle+tail fusion — BUILT, EXACT, and CLOSED negative

The resident-tile reopening was carried all the way into the production
kernel set (`-use FUSED31=1`, default off; src/cl/fusedmidtail31.cl plus
host plumbing in Gpu.cpp/replay).  Design: the Hermitian-pair closure of
tail lines over width-lines {w, WIDTH-w} fits one 64-KiB pair tile, so
(WIDTH/2+1) blocks fuse fftMiddleInGF31 + tailSquareGF31(+Zero),
running the production double-wide pair flow from the resident tile (the
fftbase LDS machinery auto-slices for concurrent pair-groups; tailutil
helpers cloned group-safe).  Correctness hazards found and solved on the
way, recorded for any future fusion attempt: (i) in-place aliasing — the
fused tail must write to the SCRATCH buffer, with fftMiddleOutGF31/
fftHinGF31 redirected to read scratch; (ii) replay sequences — mul bottom
halves and midIn-only replays must stay unfused (fuse only replays
containing BOTH KMIDIN and KTAILSQUARE); (iii) production TAIL_KERNELS=2
is the double-wide single-kernel flow with gt_line = group id.  Exact
2k/100k residues in all modes tested.

Measured (alternating 100k, 600 W, 595.84): control 144.7 us; fused
170.0 us at 256 threads (178.8 at 128).  `-time`: kFusedMidTail31 = 67-70
dilated us vs ~45 for the two kernels it replaces.  **Mechanism: the
INPLACE layout interleaves 16 widths per 256-B segment, so a
pair-resident tile (2 non-adjacent widths) reads its 64 KiB with ~16x
sector amplification and cannot coalesce — the contention-proxy's
CONTIGUOUS tiles do not transfer to this layout.**  Traffic deletion
(32 MiB) cannot pay for uncoalesced access + 1-block/SM occupancy + pair
serialization.  Decision: **the resident-fusion route is closed for the
production INPLACE layout**; reopening requires an inter-kernel layout
redesign (carryFused-side width de-interleaving) — architecture surgery,
not a kernel patch.  Code retained as opt-in evidence.

### Audit shortlist status update (post-inventory)

The three flagged benches (q24 overlap, resident tile, radix-7) were
**never committed** — they existed only in Sol's working tree and died
with that rental.  The ledger records their constants, algorithms, and
protocols verbatim, and `src/cuda/m61_near61_overlap_bench.cu` survives
as the harness template (compiles clean on CUDA 13.0/sm_120 with
`-ccbin g++`; nvcc's default host compiler is broken on this box).
Revised costs: q24 gate ~1-1.5 d reconstruction (LOW-MED risk), resident
tile 2-4 d (MED-HIGH fidelity risk; re-anchor the isolated 1.02-1.03x
ratio before adding contention), radix-7 3-4 d (MED-HIGH; resolve the
tension with the 600 W odd-radix closure first).  **Process rule added:
experiment benches must be committed, not just described in the ledger.**

### Quick probes: -lmc and TAIL_TRIGS61 — both closed

- **-lmc 405 (only alternative memory state)**: 404.6 us/iter — memory
  starvation also collapses SM clocks (877 MHz at 150 W).  2.6x slower,
  residue exact.  No usable memory-power knob exists on this board.
- **TAIL_TRIGS61=1/2 (generate M61 tail trigs instead of loading)**: +0.97
  and +1.12 us vs control (3 alternating rounds, exact).  The
  audit-derived candidate closes NEGATIVE: table loads beat regeneration
  in the real tail — consistent with the register-cap closure (the tail
  has no spare registers/ALU for trig chains; its trig loads evidently
  hit cache well enough).  Shortlist item 4 done.

### Two-worker break-even map across power points — MEASURED

`-prps 136279841,136279879 -workers 2` (TAIL_KERNELS=2 default), 200k per
exponent per power point, telemetry per run; residues exact and identical
across all points (e1's 200k residue matches the 1-worker sweep run):

| PL | 2w us/iter (e1/e2) | 2w SM MHz | 2w agg it/s | 1w it/s | delta |
|---:|---|---:|---:|---:|---:|
| 300 W | 486.2 / 485.2 | 1326 | 4,118 | 4,641 | -11.3% |
| 400 W | 395.5 / 395.7 | 1636 | 5,056 | 5,476 | -7.7% |
| 500 W | 345.1 / 344.7 | 1877 | 5,798 | 6,120 | -5.3% |
| 600 W | 307.1 / 306.9 | 2133 | 6,515 | 6,522 | **-0.1%** |

- **The break-even sits exactly at this board's 600 W ceiling**: the
  two-worker clock ratio f_2w/f_1w = 2133/2338 = 0.912 lands at the 0.914
  break-even threshold derived on the Workstation box.  The registry's
  prediction ("on a still-higher-ceiling board the sign plausibly flips")
  is confirmed and sharpened: parity at 600 W, monotonic improvement with
  power (-11.3 -> -0.1%), so any >600 W or multi-GPU host tips positive.
- At the per-watt optimum (300 W) two workers cost -11.3%: efficiency
  fleets run one worker per capped card.  This box's map also beats the
  Workstation's single 600 W point (-8.1% there vs -0.1% here at equal
  power) — better cooling retains more clock under doubled residency.

## Experiment log

### 2026-08-12: baseline reverification

Rebuilt `build-cuda/prpll` from the inherited working tree (all experimental
switches off by default).  100k-iteration production run reproduced Sol's
timings and residues exactly (201.5 us/iteration, `52775eea4730be87` at
100,000).  The machine and code state are unchanged; matched controls from
Sol's log remain valid references.

### 2026-08-12: exact detached M31 width edge (Lead 1) — implementation

Registry check: Sol's `DETACHED_M31_EDGE=1` was a correctness-disabled timing
scaffold only (186.3 us vs 205.1 us control); the exact external width kernels
were never implemented because the free-transform bound missed the 180-us
gate.  This experiment completes the exact version, targeting a real
end-to-end speedup rather than the 180-us gate.

Implementation (all opt-in via `DETACHED_M31_EDGE=1` on FFT3161):
- `bufDetach31`: one flat GF31 staging buffer (16 MiB at 4M).
- Producer A: existing `fftWGF31` (inverse width, `readCarryFusedLine` +
  `fft_WIDTH` + flat write) replayed as new bottom-half kernel `KDETACHA` at
  the end of the M31 cache group, writing `bufDetach31`.  Under MULTI_Q the
  M31 group runs on the main queue, so A overlaps the M61 tail/middle-out on
  the aux queue and finishes before the queues merge ahead of carryFused.
- Detached `carryFused` (FFT3161): skips both M31 width transforms and the
  `u31[NW]` live vector; reads `bufDetach31` flat per coefficient, writes the
  carried, forward-weighted pre-width M31 residues back to the same flat slot
  (same thread reads then writes its own slots — race-free; B(k) precedes
  A(k+1) on the in-order main queue, so one staging buffer suffices).
- Consumer B: new kernel `fftWOut31` (flat read + `fft_WIDTH2` +
  `writeCarryFusedLine`), replayed as `KDETACHB` at the START of the next M31
  cache group, i.e. after the split-queue event is recorded — so the M61 aux
  queue does NOT wait on it, and it precedes `fftMiddleInGF31` on main.
- Guards: requires in-place, short fused carry, no GOLD_PAIR, no L2_STRIPING,
  no GRAPHS.  The FFT323161 spelling remains the old timing-only scaffold.

Correctness: 4,000-iteration smoke run reproduced the exact production
iteration-2,000 residue `05d6515c416b83e2` with the Gerbicz check passing.

#### Mode 1 (both widths detached): exact, but a measured 10% loss

Matched alternating fresh 50k runs (control, trial, control, trial), all
residues identical to production at 2k/10k/20k/30k/40k/50k:

| Run | 10k block | 20k | 30k | 40k |
|---|---:|---:|---:|---:|
| control A | 197.3 | 202.1 | 202.1 | 202.4 |
| detached A | 217.4 | 222.7 | 222.5 | 222.8 |
| control B | 198.6 | 203.4 | 203.6 | 204.3 |
| detached B | 220.0 | 225.9 | 226.0 | 226.6 |

The exact mode loses ~20.5 us/iteration.  `-time` profiling explains it and
**corrects the inherited slack model**:

| Kernel (event-dilated us) | control | detached |
|---|---:|---:|
| kCarryFused | 71.5 | 59.7 |
| ktailSquareGF31 | 39.9 | 68.7 |
| kfftMidInGF31 | 17.9 | 29.6 |
| kfftMidOutGF31 | 50.8 | 22.9 |
| fftWGF31 (A, inverse width) | — | 15.0 |
| fftWOut31 (B, forward width) | — | 21.5 |

The serial fused-carry saving is real (-11.8 us, matching Sol's scaffold),
but the main queue is NOT 46 us idle: under concurrent execution the M31
kernels dilate 2-3x versus their isolated Nsight medians (main-queue busy is
~180 us of a 202-us iteration; only ~22 us idle).  Adding A+B's ~36 dilated
us to that queue overwhelms the 11.8-us serial gain.  **Lesson: PRPLL's
"hidden" M31 bottom half is only nominally hidden — both queues are
work-saturated, so relocating serial work into a parallel queue is paid at
the dilated rate, and any extra total work (64 MiB/iteration of staging
traffic here) is paid nearly in full.**  This also retroactively explains why
Sol's scaffold gain could never be realized exactly, and why batching and
sidecar transforms always lost: there is no idle execution reservoir.

#### Mode 2 (inverse width only): also a measured loss

`DETACHED_M31_EDGE=2` detaches only the inverse width (producer A, no B),
hoping A's ~15 dilated us would absorb into the ~22-us main-queue idle before
the merge.  Matched alternating 50k runs, all residues exact:

| Run | 10k | 20k | 30k | 40k |
|---|---:|---:|---:|---:|
| control 1 | 197.4 | 202.3 | 202.5 | 202.7 |
| inverse-only 1 | 207.3 | 212.1 | 212.3 | 212.9 |
| control 2 | 198.6 | 203.6 | 203.8 | 203.6 |
| inverse-only 2 | 210.1 | 215.3 | 215.5 | 216.3 |

Loss of ~10-12 us.  The idle window absorbs essentially nothing: A depends on
`fftMiddleOutGF31` (the last M31 kernel) so it lands at the end of the M31
chain where the M61 merge is imminent, and its 32 MiB of staging traffic plus
the flat reads inside the carry loop cost nearly full price.

Decision: **reject the detached M31 width edge in both exact forms.**  Keep
`DETACHED_M31_EDGE=1/2` as correct opt-in evidence (exact residues through
50k).  Do not retry by moving A/B to other queues: the aux queue has zero
idle outside carryFused's own execution window, which A and B cannot use
because both depend on or feed that very kernel.  Reopen only with a design
that reduces *total* M31 width work (not its location), or if hardware gains
a genuinely idle execution reservoir.  This closes the last unexploited
opening left by Sol's 186.3-us scaffold: the scaffold bound was real but
unreachable because it deleted work rather than relocating it.

### 2026-08-12: 600 W box — power-limit scaling investigation (see "600 W phase")

Machine change, not a code change.  Runs and protocol, all with the
unmodified production binary (commit 3071f48):

- `setup.sh` bootstrap validation: 100k exact (2k/100k residues match).
- `.pl600-steady.ul7Cke`: 5M iterations, production config, 1-s telemetry
  (power/SM clock/temp/throttle reasons).  Exact through 5M; steady state
  148.9 us at 2210 MHz / 600.0 W / 83-84 C, SW Power Cap the only active
  throttle reason for the entire run.
- `.time600.jCNBlT`: 100k with `-time` at matched thermal state; kernel
  table in the section above.
- `.w2tk2.qU8OAL`, `.w2tk3.pSUr5n`, `.w2clk.B5yGhE`: two-worker retest
  (`-prps 136279841,136279879 -workers 2`), 200k/exponent + clock capture;
  aggregate -8.1% (TK=2) / -9.1% (TK=3); 2-worker clock 1845-1852 MHz.
- `.tune600.uAnxPd`: full `-tune ntt`; `.v600-{C,A,B}{1,2,3}.*`: three-arm
  alternating 9x100k verification of its deltas (control / tune line /
  commented toggles), residues exact everywhere, no reproducible gain.

Conclusions recorded in the "600 W phase" section: perf ∝ P^0.49 (1.407x
sustained from 2.00x power); t(f) = 6.0 us + 315,842/f fits all points
1552-2467 MHz within ~1.5 us (~96% pure SM-clock scaling); still
power-limited at the 600 W board maximum; saturation model, tune optimum,
and the two-worker rejection all carry over; NCU still blocked; `-pl`/`-lgc`
still permission-blocked, so the limit itself cannot be swept from this
container.

### Toolkit NVRTC 13.2 and the clock-offset endgame (post-closure addendum)

Two further routes toward 135 us after the fusion closure:

1. **NVRTC 13.2 PTX codegen (MEASURED, adoptable): 143.2-143.3 us, exact**
   (`LD_LIBRARY_PATH=/usr/local/cuda-13.2/lib64`, package cuda-nvrtc-13-2).
   The 13.0-vs-13.2 toolkit explains ~1.4 us of the +2.9% k-residual vs
   the Workstation box (same driver JIT, different PTX).  Adopt.
2. **NVML SM clock offset (V/f shift, undervolt-equivalent): prepared,
   blocked by permissions** — scratchpad/clkoff.c
   (nvmlDeviceSetClockOffsets, SM, P0).  Arithmetic: at k~322,000 (with
   NVRTC 13.2) sustained 2450 MHz gives ~136.6; 2520 gives ~133.0.
   Offsets +60..+200 with Gerbicz/residue gating per step are the last
   route to 135 on this box; offsets reset on reboot.

Final permitted increment: `-lgc 2430,2430` biases the governor +20 MHz
(2392 sustained) -> **142.8 us exact** with NVRTC 13.2.  This is the
box's floor within the permission envelope; the residue-gated NVML
offset ladder (user-gated) is the sole path below it.
