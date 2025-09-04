// Copyright (C) Mihai Preda and George Woltman

#include "base.cl"
#include "math.cl"
#include "fft-middle.cl"
#include "middle.cl"

#if FFT_FP64

KERNEL(IN_WG) fftMiddleIn(P(T2) out, CP(T2) in, Trig trig) {
  T2 u[MIDDLE];
  
  u32 SIZEY = IN_WG / IN_SIZEX;

  u32 N = WIDTH / IN_SIZEX;
  
  u32 g = get_group_id(0);
  u32 gx = g % N;
  u32 gy = g / N;

  u32 me = get_local_id(0);
  u32 mx = me % IN_SIZEX;
  u32 my = me / IN_SIZEX;

  u32 startx = gx * IN_SIZEX;
  u32 starty = gy * SIZEY;

  u32 x = startx + mx;
  u32 y = starty + my;

  readMiddleInLine(u, in, y, x);

  middleMul2(u, x, y, 1, trig);

  fft_MIDDLE(u);

  middleMul(u, y, trig);

#if MIDDLE_IN_LDS_TRANSPOSE
  // Transpose the x and y values
  local T lds[IN_WG / 2 * (MIDDLE <= 8 ? 2 * MIDDLE : MIDDLE)];
  middleShuffle(lds, u, IN_WG, IN_SIZEX);
  out += me;  // Threads write sequentially to memory since x and y values are already transposed
#else
  // Adjust out pointer to effect a transpose of x and y values
  out += mx * SIZEY + my;
#endif

  writeMiddleInLine(out, u, gy, gx);
}

#endif



/**************************************************************************/
/*          Similar to above, but for an NTT based on GF(M31^2)           */
/**************************************************************************/

#if NTT_GF31

KERNEL(IN_WG) fftMiddleInGF31(P(GF31) out, CP(GF31) in, TrigGF31 trig) {
  GF31 u[MIDDLE];

  in += DISTGF31, out += DISTGF31;

  u32 SIZEY = IN_WG / IN_SIZEX;

  u32 N = WIDTH / IN_SIZEX;
  
  u32 g = get_group_id(0);
  u32 gx = g % N;
  u32 gy = g / N;

  u32 me = get_local_id(0);
  u32 mx = me % IN_SIZEX;
  u32 my = me / IN_SIZEX;

  u32 startx = gx * IN_SIZEX;
  u32 starty = gy * SIZEY;

  u32 x = startx + mx;
  u32 y = starty + my;

  readMiddleInLine(u, in, y, x);

  middleMul2(u, x, y, trig);

  fft_MIDDLE(u);

  middleMul(u, y, trig);

#if MIDDLE_IN_LDS_TRANSPOSE
  // Transpose the x and y values
  local Z31 lds[IN_WG / 2 * (MIDDLE <= 8 ? 2 * MIDDLE : MIDDLE)];
  middleShuffle(lds, u, IN_WG, IN_SIZEX);
  out += me;  // Threads write sequentially to memory since x and y values are already transposed
#else
  // Adjust out pointer to effect a transpose of x and y values
  out += mx * SIZEY + my;
#endif

  writeMiddleInLine(out, u, gy, gx);
}

#endif



/**************************************************************************/
/*          Similar to above, but for an NTT based on GF(M61^2)           */
/**************************************************************************/

#if NTT_GF61

KERNEL(IN_WG) fftMiddleInGF61(P(GF61) out, CP(GF61) in, TrigGF61 trig) {
  GF61 u[MIDDLE];

  in += DISTGF61, out += DISTGF61;

  u32 SIZEY = IN_WG / IN_SIZEX;

  u32 N = WIDTH / IN_SIZEX;
  
  u32 g = get_group_id(0);
  u32 gx = g % N;
  u32 gy = g / N;

  u32 me = get_local_id(0);
  u32 mx = me % IN_SIZEX;
  u32 my = me / IN_SIZEX;

  u32 startx = gx * IN_SIZEX;
  u32 starty = gy * SIZEY;

  u32 x = startx + mx;
  u32 y = starty + my;

  readMiddleInLine(u, in, y, x);

  middleMul2(u, x, y, trig);

  fft_MIDDLE(u);

  middleMul(u, y, trig);

#if MIDDLE_IN_LDS_TRANSPOSE
  // Transpose the x and y values
  local Z61 lds[IN_WG / 2 * (MIDDLE <= 8 ? 2 * MIDDLE : MIDDLE)];
  middleShuffle(lds, u, IN_WG, IN_SIZEX);
  out += me;  // Threads write sequentially to memory since x and y values are already transposed
#else
  // Adjust out pointer to effect a transpose of x and y values
  out += mx * SIZEY + my;
#endif

  writeMiddleInLine(out, u, gy, gx);
}

#endif
