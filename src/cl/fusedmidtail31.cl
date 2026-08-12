// Copyright (C) Mihai Preda and George Woltman
// FUSED31 experiment (2026-08-12): fuse fftMiddleInGF31 + tailSquareGF31
// (+tailSquareZeroGF31) into one pair-resident kernel for the M31 plane.
//
// Motivation (fable ledger, "Resident M61 tile under contention"): the
// middle->tail global round trip (16 MiB write + 16 MiB read per iteration
// on the M31 plane) is free in isolated benchmarks but real under the
// production co-run's DRAM contention.  This kernel deletes it: the
// Hermitian-pair closure of tail lines over width-lines {w, WIDTH-w} fits
// in a 64-KiB resident tile (2 x MIDDLE x SMALL_HEIGHT GF31), so one block
// holds width-lines w and WIDTH-w, performs middle-in on both, then runs
// the production double-wide tail pair flow over the 8 resident pairs,
// two pair-groups (of G_H*2 lanes) at a time.
//
// The full 3-kernel fusion is layout-impossible (fftMiddleOut's columns
// gather transPos'd lines from 8 width-lines across the plane), so
// middle-out stays a separate kernel.  IN-PLACE HAZARD: without the
// middle/tail kernel boundary, tail writes would clobber other blocks'
// unread carryFused-layout data; therefore the tail output goes to the
// SCRATCH buffer (host passes buf3 as 'out'), and fftMiddleOutGF31 (and
// fftHinGF31) read from scratch in FUSED31 mode.  Same traffic, no race.
//
// Grid: (WIDTH/2 + 1) blocks of G_H*8 = 512 threads.  Block b owns
// width-lines b and (WIDTH-b)%WIDTH (b=0 and b=WIDTH/2 are self-paired,
// half work).  The fftbase machinery slices its LDS by get_local_id/WG,
// so an 8-slice lds array serves the four concurrent pair-groups without
// modification; only the tailutil double-wide helpers (which assume a
// 128-thread block) are cloned below with an explicit pair-lane argument.

#include "base.cl"
#include "fft-middle.cl"
#include "fftheight.cl"
#include "tailutil.cl"
#include "middle.cl"

#if NTT_GF31 && INPLACE == 1 && !L2_STRIPING

// Dynamic shared memory (matches carrymiddle.cl's convention)
#ifndef CM_SHARED
#define CM_SHARED(name) extern __shared__ unsigned char name[]
#endif

// --- group-safe clones of the tailutil double-wide helpers ---
// pairMe is the lane within one pair-group [0, 2*G_H); ldsPair is the base
// of the pair-group's two LDS slices.

void revCrossLineG(local GF31* ldsPair, GF31 *u, u32 pairMe) {
  u32 lowMe = pairMe % WG;
  u32 revLowMe = WG - 1 - lowMe;
  local GF31 *ldsOut = ldsPair;
  local GF31 *ldsIn = ldsPair;
  if (pairMe < WG) ldsOut += LDS_BYTES / sizeof(GF31);
  else ldsIn += LDS_BYTES / sizeof(GF31);
  bar();
  for (u32 i = 0; i < NH/2; ++i) { ldsOut[WG * (NH/2 - 1 - i) + revLowMe] = u[i + NH/2]; }
  bar();
  for (u32 i = 0; i < NH/2; ++i) { u[i + NH/2] = ldsIn[WG * i + lowMe]; }
}

void reverse2G(local GF31 *lds, GF31 *u, u32 pairMe) {
  u32 lowMe = pairMe % WG;
  if (pairMe >= WG) lds += LDS_BYTES / sizeof(GF31);
  bar();
  for (u32 i = 0; i < NH/2; ++i) { lds[((NH/2 - i) * WG - (pairMe >= WG ? 1 : 0) - lowMe) % (NH/2 * WG)] = u[NH/2 + i]; }
  bar();
  for (u32 i = 0; i < NH/2; ++i) { u[NH/2 + i] = lds[i * WG + lowMe]; }
}

// --- copied from tailsquare.cl (GF31 section) with an explicit pair lane ---
void OVERLOAD onePairSq(GF31* pa, GF31* pb, GF31 t_squared, const u32 t_squared_type) {
  GF31 a = *pa, b = *pb;
  GF31 b2t2, c, d;

  X2conjb(a, b);
  b2t2 = cmul(csq(b), t_squared);     // b2t2 = b^2 * t_squared
  if (t_squared_type == 0)            // mul t_squared by 1
    c = csq_sub(a, b2t2);             // a^2 - (b^2 * t_squared)
  if (t_squared_type == 1)            // mul t_squared by i
    c = csq_subi(a, b2t2);            // a^2 - i*(b^2 * t_squared)
  if (t_squared_type == 2)            // mul t_squared by -1
    c = csq_add(a, b2t2);             // a^2 - -1*(b^2 * t_squared)
  if (t_squared_type == 3)            // mul t_squared by -i
    c = csq_addi(a, b2t2);            // a^2 - -i*(b^2 * t_squared)
  d = mul2(cmul(a, b));
  X2_conjb(c, d);
  *pa = SWAP_XY(c), *pb = SWAP_XY(d);
}

void pairSqG(u32 N, GF31 *u, GF31 *v, GF31 base_squared, bool special, u32 pairMe) {
  for (i32 i = 0; i < NH / 4; ++i, base_squared = mul_t8(base_squared)) {
    if (special && i == 0 && pairMe == 0) {
      u[i] = SWAP_XY(mul2(foo(u[i])));
      v[i] = SWAP_XY(shl(csq(v[i]), 2));
    } else {
      onePairSq(&u[i], &v[i], base_squared, 0);
    }
    if (N == NH) {
      onePairSq(&u[i+NH/2], &v[i+NH/2], base_squared, 2);
    }
    onePairSq(&u[i+NH/4], &v[i+NH/4], base_squared, 1);
    if (N == NH) {
      onePairSq(&u[i+3*NH/4], &v[i+3*NH/4], base_squared, 3);
    }
  }
}

// Special pairSq for double-wide line 0
void pairSq2G_special(GF31 *u, GF31 base_squared, u32 pairMe) {
  for (i32 i = 0; i < NH / 4; ++i, base_squared = mul_t8(base_squared)) {
    if (i == 0 && pairMe == 0) {
      u[0] = SWAP_XY(mul2(foo(u[0])));
      u[NH/2] = SWAP_XY(shl(csq(u[NH/2]), 2));
    } else {
      onePairSq(&u[i], &u[NH/2+i], base_squared, 0);
    }
    onePairSq(&u[i+NH/4], &u[NH/2+i+NH/4], base_squared, 1);
  }
}
// --- end copies ---

KERNEL(G_H * 4) fusedMidTail31(P(T2) out, P(T2) in, u32 base, Trig trigM, Trig trigH) {
  CM_SHARED(rawShared);
  local GF31 * const tile = (local GF31 *) rawShared;      // 2*MIDDLE*SMALL_HEIGHT GF31 = 64 KiB dynamic
  local GF31 lds[4 * LDS_BYTES / sizeof(GF31)];            // 4 auto-selected slices (fftbase offsets by me/WG)
  const u32 H = ND / SMALL_HEIGHT;

  CP(GF31) in31 = (CP(GF31)) (in + DISTGF31);              // carryFused output (read-only here)
  P(GF31) out31 = (P(GF31)) (out + DISTGF31);              // tail output -> scratch buffer
  TrigGF31 trig31 = (TrigGF31) (trigM + DISTMTRIGGF31);
  TrigGF31 smallTrig31 = (TrigGF31) (trigH + DISTHTRIGGF31);

  u32 b = get_group_id(0);                                 // 0 .. WIDTH/2 inclusive
  u32 me = get_local_id(0);
  u32 w1 = b, w2 = (WIDTH - b) % WIDTH;
  u32 nW = (w1 == w2) ? 1 : 2;

  // P1: middle-in for the resident width line(s), written straight into the
  // tile (production's middleShuffle exists only to coalesce global writes).
  for (u32 slot = 0; slot < nW; ++slot) {
    u32 w = slot ? w2 : w1;
    for (u32 h = me; h < SMALL_HEIGHT; h += G_H * 4) {
      GF31 uu[MIDDLE];
      readMiddleInLine(uu, in31, h, w);
      middleMul2(uu, w, h, trig31);
      fft_MIDDLE_IN(uu);
      middleMul(uu, h, trig31);
      for (u32 m = 0; m < MIDDLE; ++m) tile[(slot * MIDDLE + m) * SMALL_HEIGHT + h] = uu[m];
    }
  }
  bar();

  // P2: the resident lines pair among themselves (line m*WIDTH+w1 pairs with
  // H - that).  Four pair-groups run the production double-wide flow
  // concurrently; nW==2 blocks have 8 pairs (2 rounds), self-paired blocks 4.
  u32 pairMe = me % (G_H * 2);
  u32 group = me / (G_H * 2);
  u32 lowMe = pairMe % G_H;
  bool isSecondHalf = pairMe >= G_H;
  local GF31 * const ldsPair = lds + group * 2 * (LDS_BYTES / sizeof(GF31));
  u32 pairs = nW == 2 ? MIDDLE : MIDDLE / 2;
  for (u32 r = 0; r < pairs / 2; ++r) {
    u32 p = r * 2 + group;
    u32 lineA = p * WIDTH + w1;
    u32 gt_line = lineA <= H - lineA ? lineA : H - lineA;  // production's line_u ordinal
    u32 line_v = gt_line ? H - gt_line : H / 2;
    u32 line = !isSecondHalf ? gt_line : line_v;

    u32 lm = line / WIDTH, lw = line % WIDTH;
    u32 slot = (lw == w1) ? 0 : 1;
    GF31 u[NH];
    for (u32 i = 0; i < NH; ++i) u[i] = tile[(slot * MIDDLE + lm) * SMALL_HEIGHT + i * G_H + lowMe];

    fft_HEIGHT1(lds, u, smallTrig31, 2, lowMe);

#if TAIL_TRIGS31 >= 1
    u32 height_trigs = SMALL_HEIGHT*1;
    GF31 trig = TFLOAD(&smallTrig31[height_trigs + lowMe]);
    GF31 mult = TSLOAD(&smallTrig31[height_trigs + G_H + gt_line*2 + isSecondHalf]);
    trig = cmul(trig, mult);
#else
    u32 height_trigs = SMALL_HEIGHT*1;
    GF31 trig = TOLOAD(&smallTrig31[height_trigs + gt_line*G_H*2 + pairMe]);
#endif

    if (gt_line == 0) {
      reverse2G(ldsPair, u, pairMe);
      pairSq2G_special(u, trig, pairMe);
      reverse2G(ldsPair, u, pairMe);
    } else {
      revCrossLineG(ldsPair, u, pairMe);
      pairSqG(NH/2, u, u + NH/2, trig, false, pairMe);
      revCrossLineG(ldsPair, u, pairMe);
    }

    fft_HEIGHT2(lds, u, smallTrig31, 2, lowMe);
    writeTailFusedLine(u, out31, transPos(line, MIDDLE, WIDTH), lowMe);
    bar();
  }
}

#endif
