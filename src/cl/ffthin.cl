// Copyright (C) Mihai Preda

#include "base.cl"
#include "math.cl"
#include "fftheight.cl"

#if FFT_FP64

// Do an FFT Height after an fftMiddleIn (which may not have fully transposed data, leading to non-sequential input)
KERNEL(G_H) fftHin(P(T2) out, CP(T2) in, Trig smallTrig) {
  local T2 lds[SMALL_HEIGHT / 2];
  
  T2 u[NH];
  u32 g = get_group_id(0);

  u32 me = get_local_id(0);

  readTailFusedLine(in, u, g, me);

#if NH == 8
  T2 w = fancyTrig_N(ND / SMALL_HEIGHT * me);
#else
  T2 w = slowTrig_N(ND / SMALL_HEIGHT * me, ND / NH);
#endif

  fft_HEIGHT(lds, u, smallTrig, w);

  write(G_H, NH, u, out, SMALL_HEIGHT * transPos(g, MIDDLE, WIDTH));
}

#endif



/**************************************************************************/
/*          Similar to above, but for an NTT based on GF(M31^2)           */
/**************************************************************************/

#if NTT_GF31

// Do an FFT Height after an fftMiddleIn (which may not have fully transposed data, leading to non-sequential input)
KERNEL(G_H) fftHinGF31(P(GF31) out, CP(GF31) in, TrigGF31 smallTrig) {
  local GF31 lds[SMALL_HEIGHT / 2];

  in += DISTGF31, out += DISTGF31;

  GF31 u[NH];
  u32 g = get_group_id(0);

  u32 me = get_local_id(0);

  readTailFusedLine(in, u, g, me);

  fft_HEIGHT(lds, u, smallTrig);

  write(G_H, NH, u, out, SMALL_HEIGHT * transPos(g, MIDDLE, WIDTH));
}

#endif



/**************************************************************************/
/*          Similar to above, but for an NTT based on GF(M61^2)           */
/**************************************************************************/

#if NTT_GF61

// Do an FFT Height after an fftMiddleIn (which may not have fully transposed data, leading to non-sequential input)
KERNEL(G_H) fftHinGF61(P(GF61) out, CP(GF61) in, TrigGF61 smallTrig) {
  local GF61 lds[SMALL_HEIGHT / 2];

  in += DISTGF61, out += DISTGF61;

  GF61 u[NH];
  u32 g = get_group_id(0);

  u32 me = get_local_id(0);

  readTailFusedLine(in, u, g, me);

  fft_HEIGHT(lds, u, smallTrig);

  write(G_H, NH, u, out, SMALL_HEIGHT * transPos(g, MIDDLE, WIDTH));
}

#endif
