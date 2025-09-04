// Copyright (C) Mihai Preda

#include "base.cl"
#include "math.cl"
#include "fftwidth.cl"
#include "middle.cl"

#if FFT_FP64

// Do the ending fft_WIDTH after an fftMiddleOut.  This is the same as the first half of carryFused.
KERNEL(G_W) fftW(P(T2) out, CP(T2) in, Trig smallTrig) {
  local T2 lds[WIDTH / 2];

  T2 u[NW];
  u32 g = get_group_id(0);

  readCarryFusedLine(in, u, g);
  fft_WIDTH(lds, u, smallTrig);  
  out += WIDTH * g;
  write(G_W, NW, u, out, 0);
}

#endif


/**************************************************************************/
/*          Similar to above, but for an NTT based on GF(M31^2)           */
/**************************************************************************/

#if NTT_GF31

KERNEL(G_W) fftWGF31(P(GF31) out, CP(GF31) in, TrigGF31 smallTrig) {
  local GF31 lds[WIDTH / 2];

  in += DISTGF31, out += DISTGF31;

  GF31 u[NW];
  u32 g = get_group_id(0);

  readCarryFusedLine(in, u, g);
  fft_WIDTH(lds, u, smallTrig);  
  out += WIDTH * g;
  write(G_W, NW, u, out, 0);
}

#endif


/**************************************************************************/
/*          Similar to above, but for an NTT based on GF(M61^2)           */
/**************************************************************************/

#if NTT_GF61

KERNEL(G_W) fftWGF61(P(GF61) out, CP(GF61) in, TrigGF61 smallTrig) {
  local GF61 lds[WIDTH / 2];

  in += DISTGF61, out += DISTGF61;

  GF61 u[NW];
  u32 g = get_group_id(0);

  readCarryFusedLine(in, u, g);
  fft_WIDTH(lds, u, smallTrig);  
  out += WIDTH * g;
  write(G_W, NW, u, out, 0);
}

#endif
