// Copyright (C) Mihai Preda

#include "base.cl"

#if PARITY_PACKED
u32 parityBallot(bool predicate) {
#if CUDA_BACKEND
  return __ballot_sync(0xffffffffu, predicate);
#else
  return predicate ? 1u : 0u;
#endif
}
#endif

#if PARITY_INIT
#if PARITY_DIRECT
void parityInitStoreExpected(P(u32) parity, u32 sourceWord, u32 bit) {
  u32 const outputPair = sourceWord & (NWORDS / 2 - 1);
  u32 const frac = ((u32) EXP * outputPair) & (NWORDS - 1);
  bool const selectsUpper = frac != 0 && frac <= NWORDS / 2;
  bool const isUpper = sourceWord >= NWORDS / 2;
  if (selectsUpper == isUpper) {
    u32 const x = outputPair / BIG_HEIGHT;
    u32 const line = outputPair - x * BIG_HEIGHT;
    parity[line * WIDTH + x] = bit & 1u;
  }
}
#endif

// Build a sidecar containing the parity of each logical pair of input words.
// The integer buffer is physically transposed; the sidecar is in logical
// Crandall--Fagin coefficient order so carry kernels can index it directly.
KERNEL(256) parityInit(P(u32) parity, CP(Word2) in) {
  u32 pair = get_global_id(0);
  if (pair >= ND) return;
#if PARITY_PACKED
  Word2 const words = in[pair];
  u32 const xMask = parityBallot(((u32) words.x & 1u) != 0);
  u32 const yMask = parityBallot(((u32) words.y & 1u) != 0);
  if ((pair & 31u) == 0) {
    parity[pair >> 5] = xMask;
    parity[ND / 32 + (pair >> 5)] = yMask;
  }
#elif PARITY_DIRECT
  u32 const line = pair / WIDTH;
  u32 const x = pair - line * WIDTH;
  u32 const logicalPair = x * BIG_HEIGHT + line;
  Word2 const words = in[pair];
  parityInitStoreExpected(parity, 2 * logicalPair, (u32) words.x);
  parityInitStoreExpected(parity, 2 * logicalPair + 1, (u32) words.y);
#elif PARITY_PREPARED || PARITY_PHYSICAL
  // Physical carry layout: x changes fastest, so both this initialization
  // and carry's later reads/writes are coalesced.
  parity[pair] = ((u32) in[pair].x & 1u) | (((u32) in[pair].y & 1u) << 1);
#else
  u32 x = pair / BIG_HEIGHT;
  u32 y = pair - x * BIG_HEIGHT;
  Word2 words = in[WIDTH * y + x];
  parity[pair] = ((u32) words.x & 1u) | (((u32) words.y & 1u) << 1);
#endif
}
#endif

#if PARITY_PREPARE
// For a weighted cyclic square, off-diagonal terms vanish modulo two.  Map
// the packed input-word parity into the exact even-coefficient parity needed
// by carry.  Input and output use physical (line-major) carry layout.
KERNEL(256) parityPrepare(CP(u32) parity, P(u32) expected) {
  u32 physical = get_global_id(0);
  if (physical >= ND) return;
  u32 line = physical / WIDTH;
  u32 x = physical - line * WIDTH;
  u32 outputPair = x * BIG_HEIGHT + line;

  u32 sourceWord = outputPair;
  u32 frac = ((u32) EXP * sourceWord) & (NWORDS - 1);
  if (frac != 0 && frac <= NWORDS / 2) sourceWord += NWORDS / 2;

  u32 sourcePair = sourceWord >> 1;
  u32 sourceX = sourcePair / BIG_HEIGHT;
  u32 sourceLine = sourcePair - sourceX * BIG_HEIGHT;
  u32 const sourcePhysical = sourceLine * WIDTH + sourceX;
#if PARITY_PACKED
  u32 const plane = (sourceWord & 1u) * (ND / 32);
  u32 const bit = (parity[plane + (sourcePhysical >> 5)] >> (sourcePhysical & 31u)) & 1u;
  u32 const mask = parityBallot(bit != 0);
  if ((physical & 31u) == 0) expected[physical >> 5] = mask;
#else
  u32 packed = parity[sourcePhysical];
  expected[physical] = (packed >> (sourceWord & 1)) & 1u;
#endif
}
#endif

#if READRESIDUE

// Because the data "in" is stored transposed, and we want to read
// a number of logically successive values, we have a very bad read access pattern
KERNEL(32) readResidue(P(Word2) out, CP(Word2) in) {
  u32 me = get_local_id(0);
  u32 k = (ND - 16 + me) % ND;
  u32 y = k % BIG_HEIGHT;
  u32 x = k / BIG_HEIGHT;
  out[me] = in[WIDTH * y + x];
}
#endif

#if SUM64
KERNEL(64) sum64(global ulong* out, u32 count, CP(Word) in) {
  ulong sum = 0;
  for (i32 p = get_global_id(0); p < count; p += get_global_size(0)) {
    sum += in[p];
  }
  u32 prev = atomic_add((global u32*)out, (u32) sum);
  u32 high = (sum + prev) >> 32;
  atomic_add(((global u32*)out) + 1, high);
}
#endif

#if ISEQUAL
// outEqual must be "true" on entry.
KERNEL(256) isEqual(global i64 *in1, global i64 *in2, P(int) outEqual) {
  for (i32 p = get_global_id(0); p < NWORDS * sizeof(Word) / sizeof(i64); p += get_global_size(0)) {
    if (in1[p] != in2[p]) {
      *outEqual = 0;
      return;
    }
  }
}
#endif

#if TEST_KERNEL
// Generate a small unused kernel so developers can look at how well individual macros assemble and optimize
kernel void testKernel(global int* in, global double* out) {
  const double TAB[8] = {M_PI/13, M_PI/17, M_PI, M_SQRT2, M_SQRT1_2, M_PI/7, M_PI*7, M_PI/15};

  int me = get_local_id(0);
  int p = me * in[me] % 8; // % 15;
  out[me] = TAB[p];
}
#endif
