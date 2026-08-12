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
