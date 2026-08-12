---
name: prpll-campaign-state
description: "PRPLL speedup campaign status, baseline numbers, saturation model, and ledger conventions"
metadata: 
  node_type: memory
  type: project
  originSessionId: 29ee1069-0039-4be4-8ffe-bebd253e2d4d
  modified: 2026-08-12T15:55:36.049Z
---

The task is speeding up PRPLL for exponent 136279841 (4M `1:512:8:512:202`
M31+M61 NTT).  Ledgers: `sol-prpll-speedup-attempts.md` (~95 rejected
experiments; the registry that must be checked before any new attempt) and
`fable-prpll-speedup-attempts.md` (continuation).  Both must be read at
session start; nothing recorded as rejected may be repeated without meeting
its reopen condition.

**Machine history:** (1) Max-Q GB202, 300 W cap: 201.5 us/100k, 209.5 steady.
(2) Workstation GB202, 600 W pinned: 148.9 us steady @ 2210 MHz; perf ∝ P^0.49;
t(f) = 6.0 us + 315,842/f.  (3) Current (2026-08-12): RTX PRO 6000 Blackwell
**Server Edition** (GB202, 188 SM, max SM clock 2430 MHz, 300–600 W range,
driver 580.126.09, CUDA 13.0) — **sudo nvidia-smi -pl AND -lgc work, NCU
counters work (RmProfilingAdminOnly=1 + passwordless sudo)**.  Validation
run exact (residues `05d6515c416b83e2`@2k, `52775eea4730be87`@100k).
Admin-host results (all in the fable ledger): perf/W knee ~450-500 W;
ECC OFF (worth +1.2%, power channel only — zero latency cost);
**driver 595.84 replaced 580.126 (+5-6%: GSP memory-path firmware, SASS
near-identical — driver major version is a first-class perf variable)**;
new baseline **144.67 us 1M-steady @600 W** (~142-143 at 100k scale),
t(f) = 5.16 + 325,116/f (+2.9% residual vs Workstation = board/bin);
two workers still -1-2% at 600 W under 595; occupancy reg-caps closed
negative; **cp.async middle-prefetch implemented (exact, opt-in
ASYNC_MID61) and closed negative (+4.6 us best arm) — co-run backfill
consumes isolated-kernel latency slack; the counters route is fully
walked.**  Audit shortlist fully walked (all benches committed): q24 charged gate —
rejection CONFIRMED (+50 us); radix-7 unified gate — rejection CONFIRMED
(edge co-run hideability measured ~12%, retro-validating the odd-radix
closures); resident-tile proxy WON under contention (-18.8 us) BUT the
production fusion (BUILT: `-use FUSED31=1`, exact residues, code in tree)
CLOSED negative +25 us: the INPLACE layout interleaves 16 widths per
256-B segment, so pair-resident tiles read uncoalesced (~16x sector
amplification) — contiguous-tile proxies do not transfer.  Reopen only
via a carryFused-side layout redesign.  **No software leads remain.
Production stands at 144.7 us steady @600 W; the 135-us goal needs
~6.7% more with no identified software path on this layout/hardware —
remaining routes are hardware (>600 W board, two workers there).**

**Key model (corrects Sol's):** both MULTI_Q queues are work-saturated — M31
kernels dilate 2-3x under concurrency; only deleting work (not relocating it)
shortens the iteration.  Do not reopen relocation/sidecar/tune ideas.
Registry reference residues: 1M `52b03a7cc55e677d`, 5M `7eca65291732df02`.

**Conventions:** matched alternating fresh-directory runs with residue checks
at every checkpoint; commit each experiment with a "Record ..." message;
scratch run dirs `.name.XXXX/` stay untracked; `setup.sh pack` migrates
artifacts to the next box.

**Session end (2026-08-12):** best exact config = production `-use` line +
`LD_LIBRARY_PATH=/usr/local/cuda-13.2/lib64` (NVRTC 13.2): **143.2 us**.
PDL: no effect.  140M probe exact (shape 4:1K:8:256:101, gains transfer).
135-us floor needs sustained ~2450+ MHz: ONLY remaining lever = NVML SM
clock offset (tool: scratchpad clkoff.c; classifier-gated to the USER —
they must run `sudo /tmp/clkoff 60` first; then ladder +60/+120/+150/+180
with 100k residue gate each, 1M confirm at best stable).

**Architectural endgame:** approach A (carry decoupling: shuttle = 6-us
algorithmic price, upstream design vindicated) and approach B (layout
co-design: cluster-DSM 17x worse, amplification L2-absorbed, inversion
analytically closed) both measured out.  The architectural space is
exhausted end to end.  Terminal: 142.8 us exact.

**Two-stage factorization: BUILT end-to-end and CLOSED negative
(2026-08-12).**  Final shape 1:512:1:4K (MIDDLE=1; 2048 is impossible —
the generic NTT ladder only does pure radix powers; 4K:1:512 loses to a
102.5-us 1-block/SM carryFused).  `-use TWO_STAGE=1,TAIL_KERNELS=0`
elides fftMiddleIn/Out on square passes; tail does direct carryFused-
layout IO + folded W_ND^(w*y) twiddles (TS_TAIL kernels; note -use keys
leak into the global define string — device gate must be a separate
macro).  EXACT at 2k/100k/1M, but 261.1 us vs 152.0 production: the
M=1 transpose amortization loss and the direct-IO sector scatter
(+55 us/field) dwarf the 96-MiB traffic saving.  The three-stage
512x8x512 wins on hardware REGIME (8-way middle amortization,
double-wide 33-KiB tail, 20-warp carry), not butterfly count.  Reopen
only on 2x-register-file or 96-KiB-static-shared hardware.  Full record
in ledger "Two-stage backend BUILT end-to-end".  Shapes 1:4K:1:512 and
1:512:1:4K remain runnable (exact) through the normal pipeline.

**Power-aware codegen (E1+E2): DONE and CLOSED (2026-08-12).**  E1 =
first per-class power table on GB202 (`src/cuda/clockcost.cu`,
`e1-clockcost/`): IADD.64 0.78 pJ/cyc and IMAD.WIDE 0.75 are the hogs
(the M61 field), LOP3 0.53 cheapest; instruction mix alone spans 510 MHz
of sustained clock at a 300 W cap; production binds at 600 W / 2313 MHz
with **1 W = 0.083 us (12 W per us)**.  E2 = NVRTC 13.0/13.2/13.3
SASS-diff + A/B: 13.3 (via new `PRPLL_PTX_VERSION` clamp; 595 JIT rejects
its PTX 9.3) and `PRPLL_JIT_OPT` levels all FLAT; the banked 13.0->13.2
-1.4 us = carryFused 4680->4504 instr; GF61 reg-reg MOVs (187/200) are
pinned by the JIT register allocator — not pullable.  Channel yield
<=0.5-1 us, captured none; reopen only on driver/ptxas major bump
(re-run `e2ab.sh`+`e2stats.py`).

**Clock offset: CLOSED by hardware (2026-08-12).**  User ran the clkoff
bootstrap; Server Edition vBIOS reports GRAPHICS P0 offset range [0,0]
(SM type unsupported) — OC is fused off on this SKU.  **CAMPAIGN
TERMINAL: 142-143 us; every software and quasi-software channel closed.
135 needs different hardware: >600 W board, offset-unlocked
Workstation-edition GB202, or next-gen part.**
