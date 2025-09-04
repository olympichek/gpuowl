// Copyright (C) Mihai Preda and George Woltman

#include "base.cl"
#include "tailutil.cl"
#include "trig.cl"
#include "fftheight.cl"

#if FFT_FP64

// Handle the final multiplication step on a pair of complex numbers.  Swap real and imaginary results for the inverse FFT.
// We used to conjugate the results, but swapping real and imaginary can save some negations in carry propagation.

void OVERLOAD onePairMul(T2* pa, T2* pb, T2* pc, T2* pd, T2 t_squared) {
  T2 a = *pa, b = *pb, c = *pc, d = *pd;

  X2conjb(a, b);
  X2conjb(c, d);

  *pa = cfma(a, c, cmul(cmul(b, d), -t_squared));
  *pb = cfma(b, c, cmul(a, d));

  X2_conjb(*pa, *pb);

  *pa = SWAP_XY(*pa), *pb = SWAP_XY(*pb);
}

void OVERLOAD pairMul(u32 N, T2 *u, T2 *v, T2 *p, T2 *q, T2 base_squared, bool special) {
  u32 me = get_local_id(0);

  for (i32 i = 0; i < NH / 4; ++i, base_squared = mul_t8(base_squared)) {
    if (special && i == 0 && me == 0) {
      u[i] = SWAP_XY(2 * foo2(u[i], p[i]));
      v[i] = SWAP_XY(4 * cmul(v[i], q[i]));
    } else {
      onePairMul(&u[i], &v[i], &p[i], &q[i], base_squared);
    }

    if (N == NH) {
      onePairMul(&u[i+NH/2], &v[i+NH/2], &p[i+NH/2], &q[i+NH/2], -base_squared);
    }

    T2 new_base_squared = mul_t4(base_squared);
    onePairMul(&u[i+NH/4], &v[i+NH/4], &p[i+NH/4], &q[i+NH/4], new_base_squared);

    if (N == NH) {
      onePairMul(&u[i+3*NH/4], &v[i+3*NH/4], &p[i+3*NH/4], &q[i+3*NH/4], -new_base_squared);
    }
  }
}

KERNEL(G_H) tailMul(P(T2) out, CP(T2) in, CP(T2) a, Trig smallTrig) {
  local T2 lds[SMALL_HEIGHT];

  T2 u[NH], v[NH];
  T2 p[NH], q[NH];

  u32 H = ND / SMALL_HEIGHT;

  u32 line1 = get_group_id(0);
  u32 line2 = line1 ? H - line1 : (H / 2);
  u32 memline1 = transPos(line1, MIDDLE, WIDTH);
  u32 memline2 = transPos(line2, MIDDLE, WIDTH);

  u32 me = get_local_id(0);
  readTailFusedLine(in, u, line1, me);
  readTailFusedLine(in, v, line2, me);

#if NH == 8
  T2 w = fancyTrig_N(ND / SMALL_HEIGHT * me);
#else
  T2 w = slowTrig_N(ND / SMALL_HEIGHT * me, ND / NH);
#endif

#if MUL_LOW
  read(G_H, NH, p, a, memline1 * SMALL_HEIGHT);
  read(G_H, NH, q, a, memline2 * SMALL_HEIGHT);
  fft_HEIGHT(lds, u, smallTrig, w);
  bar();
  fft_HEIGHT(lds, v, smallTrig, w);
#else
  readTailFusedLine(a, p, line1, me);
  readTailFusedLine(a, q, line2, me);
  fft_HEIGHT(lds, u, smallTrig, w);
  bar();
  fft_HEIGHT(lds, v, smallTrig, w);
  bar();
  fft_HEIGHT(lds, p, smallTrig, w);
  bar();
  fft_HEIGHT(lds, q, smallTrig, w);
#endif

  T2 trig = slowTrig_N(line1 + me * H, ND / NH);

  if (line1 == 0) {
    reverse(G_H, lds, u + NH/2, true);
    reverse(G_H, lds, p + NH/2, true);
    pairMul(NH/2, u,  u + NH/2, p, p + NH/2, trig, true);
    reverse(G_H, lds, u + NH/2, true);

    T2 trig2 = cmulFancy(trig, TAILT);
    reverse(G_H, lds, v + NH/2, false);
    reverse(G_H, lds, q + NH/2, false);
    pairMul(NH/2, v,  v + NH/2, q, q + NH/2, trig2, false);
    reverse(G_H, lds, v + NH/2, false);
  } else {
    reverseLine(G_H, lds, v);
    reverseLine(G_H, lds, q);
    pairMul(NH, u, v, p, q, trig, false);
    reverseLine(G_H, lds, v);
  }

  bar();
  fft_HEIGHT(lds, v, smallTrig, w);
  bar();
  fft_HEIGHT(lds, u, smallTrig, w);
  writeTailFusedLine(v, out, memline2, me);
  writeTailFusedLine(u, out, memline1, me);
}

#endif



/**************************************************************************/
/*          Similar to above, but for an NTT based on GF(M31^2)           */
/**************************************************************************/

#if NTT_GF31

void OVERLOAD onePairMul(GF31* pa, GF31* pb, GF31* pc, GF31* pd, GF31 t_squared) {
  GF31 a = *pa, b = *pb, c = *pc, d = *pd;

  X2conjb(a, b);
  X2conjb(c, d);

  *pa = sub(cmul(a, c), cmul(cmul(b, d), t_squared));
  *pb = add(cmul(b, c), cmul(a, d));

  X2_conjb(*pa, *pb);
  *pa = SWAP_XY(*pa), *pb = SWAP_XY(*pb);
}

void OVERLOAD pairMul(u32 N, GF31 *u, GF31 *v, GF31 *p, GF31 *q, GF31 base_squared, bool special) {
  u32 me = get_local_id(0);

  for (i32 i = 0; i < NH / 4; ++i, base_squared = mul_t8(base_squared)) {
    if (special && i == 0 && me == 0) {
      u[i] = SWAP_XY(mul2(foo2(u[i], p[i])));
      v[i] = SWAP_XY(shl(cmul(v[i], q[i]), 2));
   } else {
      onePairMul(&u[i], &v[i], &p[i], &q[i], base_squared);
    }

    if (N == NH) {
      onePairMul(&u[i+NH/2], &v[i+NH/2], &p[i+NH/2], &q[i+NH/2], neg(base_squared));
    }

    GF31 new_base_squared = mul_t4(base_squared);
    onePairMul(&u[i+NH/4], &v[i+NH/4], &p[i+NH/4], &q[i+NH/4], new_base_squared);

    if (N == NH) {
      onePairMul(&u[i+3*NH/4], &v[i+3*NH/4], &p[i+3*NH/4], &q[i+3*NH/4], neg(new_base_squared));
    }
  }
}

KERNEL(G_H) tailMulGF31(P(GF31) out, CP(GF31) in, CP(GF31) a, TrigGF31 smallTrig) {
  local GF31 lds[SMALL_HEIGHT];

  in += DISTGF31, out += DISTGF31, a += DISTGF31;

  GF31 u[NH], v[NH];
  GF31 p[NH], q[NH];

  u32 H = ND / SMALL_HEIGHT;

  u32 line1 = get_group_id(0);
  u32 line2 = line1 ? H - line1 : (H / 2);
  u32 memline1 = transPos(line1, MIDDLE, WIDTH);
  u32 memline2 = transPos(line2, MIDDLE, WIDTH);

  u32 me = get_local_id(0);
  readTailFusedLine(in, u, line1, me);
  readTailFusedLine(in, v, line2, me);

#if MUL_LOW
  read(G_H, NH, p, a, memline1 * SMALL_HEIGHT);
  read(G_H, NH, q, a, memline2 * SMALL_HEIGHT);
  fft_HEIGHT(lds, u, smallTrig);
  bar();
  fft_HEIGHT(lds, v, smallTrig);
#else
  readTailFusedLine(a, p, line1, me);
  readTailFusedLine(a, q, line2, me);
  fft_HEIGHT(lds, u, smallTrig);
  bar();
  fft_HEIGHT(lds, v, smallTrig);
  bar();
  fft_HEIGHT(lds, p, smallTrig);
  bar();
  fft_HEIGHT(lds, q, smallTrig);
#endif

  // Calculate number of trig values used by fft_HEIGHT (see genSmallTrigCombo in trigBufCache.cpp)
  // The trig values used here are pre-computed and stored after the fft_HEIGHT trig values.
  u32 height_trigs = SMALL_HEIGHT*1;
#if TAIL_TRIGS >= 1
  GF31 trig = smallTrig[height_trigs + me];                    // Trig values for line zero, should be cached
#if SINGLE_WIDE
  GF31 mult = smallTrig[height_trigs + G_H + line1];
#else
  GF31 mult = smallTrig[height_trigs + G_H + line1 * 2];
#endif
  trig = cmul(trig, mult);
#else
#if SINGLE_WIDE
  GF31 trig = NTLOAD(smallTrig[height_trigs + line1*G_H + me]);
#else
  GF31 trig = NTLOAD(smallTrig[height_trigs + line1*2*G_H + me]);
#endif
#endif

  if (line1 == 0) {
    reverse(G_H, lds, u + NH/2, true);
    reverse(G_H, lds, p + NH/2, true);
    pairMul(NH/2, u,  u + NH/2, p, p + NH/2, trig, true);
    reverse(G_H, lds, u + NH/2, true);

    GF31 trig2 = cmul(trig, TAILTGF31);
    reverse(G_H, lds, v + NH/2, false);
    reverse(G_H, lds, q + NH/2, false);
    pairMul(NH/2, v,  v + NH/2, q, q + NH/2, trig2, false);
    reverse(G_H, lds, v + NH/2, false);
  } else {
    reverseLine(G_H, lds, v);
    reverseLine(G_H, lds, q);
    pairMul(NH, u, v, p, q, trig, false);
    reverseLine(G_H, lds, v);
  }

  bar();
  fft_HEIGHT(lds, v, smallTrig);
  bar();
  fft_HEIGHT(lds, u, smallTrig);
  writeTailFusedLine(v, out, memline2, me);
  writeTailFusedLine(u, out, memline1, me);
}

#endif



/**************************************************************************/
/*          Similar to above, but for an NTT based on GF(M61^2)           */
/**************************************************************************/

#if NTT_GF61

void OVERLOAD onePairMul(GF61* pa, GF61* pb, GF61* pc, GF61* pd, GF61 t_squared) {
  GF61 a = *pa, b = *pb, c = *pc, d = *pd;

  X2conjb(a, b);
  X2conjb(c, d);

  *pa = sub(cmul(a, c), cmul(cmul(b, d), t_squared));
  *pb = add(cmul(b, c), cmul(a, d));

  X2_conjb(*pa, *pb);
  *pa = SWAP_XY(*pa), *pb = SWAP_XY(*pb);
}

void OVERLOAD pairMul(u32 N, GF61 *u, GF61 *v, GF61 *p, GF61 *q, GF61 base_squared, bool special) {
  u32 me = get_local_id(0);

  for (i32 i = 0; i < NH / 4; ++i, base_squared = mul_t8(base_squared)) {
    if (special && i == 0 && me == 0) {
      u[i] = SWAP_XY(mul2(foo2(u[i], p[i])));
      v[i] = SWAP_XY(shl(cmul(v[i], q[i]), 2));
   } else {
      onePairMul(&u[i], &v[i], &p[i], &q[i], base_squared);
    }

    if (N == NH) {
      onePairMul(&u[i+NH/2], &v[i+NH/2], &p[i+NH/2], &q[i+NH/2], neg(base_squared));
    }

    GF61 new_base_squared = mul_t4(base_squared);
    onePairMul(&u[i+NH/4], &v[i+NH/4], &p[i+NH/4], &q[i+NH/4], new_base_squared);

    if (N == NH) {
      onePairMul(&u[i+3*NH/4], &v[i+3*NH/4], &p[i+3*NH/4], &q[i+3*NH/4], neg(new_base_squared));
    }
  }
}

KERNEL(G_H) tailMulGF61(P(GF61) out, CP(GF61) in, CP(GF61) a, TrigGF61 smallTrig) {
  local GF61 lds[SMALL_HEIGHT];

  in += DISTGF61, out += DISTGF61, a += DISTGF61;

  GF61 u[NH], v[NH];
  GF61 p[NH], q[NH];

  u32 H = ND / SMALL_HEIGHT;

  u32 line1 = get_group_id(0);
  u32 line2 = line1 ? H - line1 : (H / 2);
  u32 memline1 = transPos(line1, MIDDLE, WIDTH);
  u32 memline2 = transPos(line2, MIDDLE, WIDTH);

  u32 me = get_local_id(0);
  readTailFusedLine(in, u, line1, me);
  readTailFusedLine(in, v, line2, me);

#if MUL_LOW
  read(G_H, NH, p, a, memline1 * SMALL_HEIGHT);
  read(G_H, NH, q, a, memline2 * SMALL_HEIGHT);
  fft_HEIGHT(lds, u, smallTrig);
  bar();
  fft_HEIGHT(lds, v, smallTrig);
#else
  readTailFusedLine(a, p, line1, me);
  readTailFusedLine(a, q, line2, me);
  fft_HEIGHT(lds, u, smallTrig);
  bar();
  fft_HEIGHT(lds, v, smallTrig);
  bar();
  fft_HEIGHT(lds, p, smallTrig);
  bar();
  fft_HEIGHT(lds, q, smallTrig);
#endif

  // Calculate number of trig values used by fft_HEIGHT (see genSmallTrigCombo in trigBufCache.cpp)
  // The trig values used here are pre-computed and stored after the fft_HEIGHT trig values.
  u32 height_trigs = SMALL_HEIGHT*1;
#if TAIL_TRIGS >= 1
  GF61 trig = smallTrig[height_trigs + me];                    // Trig values for line zero, should be cached
#if SINGLE_WIDE
  GF61 mult = smallTrig[height_trigs + G_H + line1];
#else
  GF61 mult = smallTrig[height_trigs + G_H + line1 * 2];
#endif
  trig = cmul(trig, mult);
#else
#if SINGLE_WIDE
  GF61 trig = NTLOAD(smallTrig[height_trigs + line1*G_H + me]);
#else
  GF61 trig = NTLOAD(smallTrig[height_trigs + line1*2*G_H + me]);
#endif
#endif

  if (line1 == 0) {
    reverse(G_H, lds, u + NH/2, true);
    reverse(G_H, lds, p + NH/2, true);
    pairMul(NH/2, u,  u + NH/2, p, p + NH/2, trig, true);
    reverse(G_H, lds, u + NH/2, true);

    GF61 trig2 = cmul(trig, U2((Z61)TAILTGF61.x, (Z61)TAILTGF61.y));
    reverse(G_H, lds, v + NH/2, false);
    reverse(G_H, lds, q + NH/2, false);
    pairMul(NH/2, v,  v + NH/2, q, q + NH/2, trig2, false);
    reverse(G_H, lds, v + NH/2, false);
  } else {
    reverseLine(G_H, lds, v);
    reverseLine(G_H, lds, q);
    pairMul(NH, u, v, p, q, trig, false);
    reverseLine(G_H, lds, v);
  }

  bar();
  fft_HEIGHT(lds, v, smallTrig);
  bar();
  fft_HEIGHT(lds, u, smallTrig);
  writeTailFusedLine(v, out, memline2, me);
  writeTailFusedLine(u, out, memline1, me);
}

#endif
