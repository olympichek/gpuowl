// Forward width edge for an already M31-weighted folded syndrome input.

#include "base.cl"
#include "fftwidth.cl"
#include "middle.cl"

#if FOLD_TRANSFORM && FFT_TYPE == FFT31
#if FOLD_FACTOR == 16
// carryFused naturally produces eight aliases per value.  Combine pairs of
// those partial bins into the sixteen-alias fold, and transpose directly into
// this shortened transform's line-major input layout.
KERNEL(256) foldCompact(P(GF31) out, CP(GF31) in) {
  u32 const logical = get_global_id(0);
  if (logical >= ND) return;
  GF31 const value = add(in[logical], in[logical + ND]);
  u32 const x = logical / BIG_HEIGHT;
  u32 const line = logical - x * BIG_HEIGHT;
  out[line * WIDTH + x] = value;
}
#else
KERNEL(256) foldCompact(P(GF31) out, CP(GF31) in) {
  u32 const logical = get_global_id(0);
  if (logical < ND) out[logical] = in[logical];
}
#endif

KERNEL(G_W) foldP(P(GF31) out, CP(GF31) in, TrigGF31 smallTrig) {
  local GF31 lds[LDS_BYTES / sizeof(GF31)];
  GF31 u[NW];

  u32 const g = get_group_id(0);
  u32 const me = get_local_id(0);
  in += g * WIDTH;
  for (u32 i = 0; i < NW; ++i) {
    u[i] = in[G_W * i + me];
  }

  fft_WIDTH(lds, u, smallTrig, 1, me);
  writeCarryFusedLine(u, out, g, me);
}
#endif
