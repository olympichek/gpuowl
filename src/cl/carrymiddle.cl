// Experimental production-shaped fusion for the steady-state M31/M61 path.
//
// The middle transform groups eight carry lines that are 512 lines apart,
// whereas carries advance through consecutive lines.  Keeping the complete
// boundary in one kernel would therefore require an unsafe grid-wide wait.
// Split it at the integer representation instead:
//
//   carryMiddleOut: middle-out -> inverse width -> CRT/sloppy carry
//   carryMiddleIn:  incoming carry -> forward width -> middle-in
//
// The bridge stores two low 32-bit words per quadratic value.  Their two
// 33rd bits are packed into warp-ballot planes, reducing the intermediate
// from 32 MiB to 16.5 MiB for a 4M transform.

#include "base.cl"
#include "fftwidth.cl"
#include "fft-middle.cl"
#include "carryutil.cl"
#include "weight.cl"
#include "middle.cl"

#if FFT_TYPE == FFT3161

#if WIDTH != 512 || SMALL_HEIGHT != 512 || MIDDLE != 8 || NW != 8 || !INPLACE
#error carryMiddle currently supports only in-place 512:8:512 FFT3161
#endif

#if EXP / NWORDS != 32
#error carryMiddle packed bridge requires 32 <= bits/word < 33
#endif

#define CM_THREADS 512
#define CM_SHARED_BYTES (MIDDLE * WIDTH * (sizeof(GF31) + sizeof(GF61)))

#ifndef CM_SKIP_MIDDLE
#define CM_SKIP_MIDDLE 0
#endif

#if CUDA_BACKEND
#define CM_SHARED(name) extern __shared__ unsigned char name[]
#else
// This path is deliberately CUDA-only.  The declaration keeps preprocessing
// intelligible on other backends, while startup rejects the mode there.
#define CM_SHARED(name) local unsigned char name[CM_SHARED_BYTES]
#endif

u32 cmBallot(bool predicate) {
#if CUDA_BACKEND
  return __ballot_sync(0xffffffffu, predicate);
#else
  return predicate ? 0xffffffffu : 0u;
#endif
}

void cmInitCounters(u32 word_index,
                    u64 *m31_counter, u64 *m61_counter,
                    u64 *m31_start, u64 *m61_start) {
  const u32 m31_log2_root_two = (u32)(((1ULL << 30) / NWORDS) % 31);
  const u32 m31_bigword_weight_shift =
    (NWORDS - EXP % NWORDS) * m31_log2_root_two % 31;
  const u32 m31_shift_step = (m31_bigword_weight_shift + 30) % 31;
  const u32 m61_log2_root_two = (u32)(((1ULL << 60) / NWORDS) % 61);
  const u32 m61_bigword_weight_shift =
    (NWORDS - EXP % NWORDS) * m61_log2_root_two % 61;
  const u32 m61_shift_step = (m61_bigword_weight_shift + 60) % 61;

  *m31_counter = comboFracBits(word_index) +
                 make_u64(word_index * m31_shift_step, 0xFFFFFFFF);
  *m61_counter = comboFracBits(word_index) +
                 make_u64(word_index * m61_shift_step, 0xFFFFFFFF);

  union { uint2 a; u64 b; } c31, c61;
  c31.b = *m31_counter;
  c61.b = *m61_counter;
  c31.a[1] %= 31;
  c61.a[1] %= 61;
  *m31_counter = c31.b;
  *m61_counter = c61.b;
  *m31_start = c31.b;
  *m61_start = c61.b;

  const u32 log2_NWORDS = 22;
  c31.a[1] = adjust_m31_weight_shift(c31.a[1] + log2_NWORDS + 1);
  c61.a[1] = adjust_m61_weight_shift(c61.a[1] + log2_NWORDS + 1);
  *m31_counter = c31.b;
  *m61_counter = c61.b;
}

// Perform the inverse middle and width edges, exact two-prime reconstruction,
// and the independent/sloppy portion of carry propagation.
KERNEL(CM_THREADS) carryMiddleOut(CP(T2) in, P(u32) packed,
                                  P(i64) carryOut, u32 posROE,
                                  Trig middleTrig, Trig widthTrig,
                                  P(uint) bufROE) {
  CM_SHARED(sharedBytes);
  local GF31 *tile31 = (local GF31 *)sharedBytes;
  local GF61 *tile61 = (local GF61 *)(sharedBytes +
                                      MIDDLE * WIDTH * sizeof(GF31));

  CP(GF31) in31 = (CP(GF31))(in + DISTGF31);
  CP(GF61) in61 = (CP(GF61))(in + DISTGF61);
  TrigGF31 middle31 = (TrigGF31)(middleTrig + DISTMTRIGGF31);
  TrigGF61 middle61 = (TrigGF61)(middleTrig + DISTMTRIGGF61);
  TrigGF31 width31 = (TrigGF31)(widthTrig + DISTWTRIGGF31);
  TrigGF61 width61 = (TrigGF61)(widthTrig + DISTWTRIGGF61);

  u32 const x = get_group_id(0);
  u32 const tid = get_local_id(0);
  u32 const middle = tid / G_W;
  u32 const lowMe = tid % G_W;
  u32 const line = middle * SMALL_HEIGHT + x;
  GF31 u31[NW];
  GF61 u61[NW];

#if CM_SKIP_MIDDLE
  // Diagnostic form: consume the normal fftMiddleOut layout.  This isolates
  // the two-kernel carry bridge from the experimental middle-stage fusion.
  readCarryFusedLine(in31, u31, line, lowMe);
  readCarryFusedLine(in61, u61, line, lowMe);
#else
  // A fixed SMALL_HEIGHT coordinate owns all WIDTH columns and all eight
  // middle outputs.  The direct logical load avoids both 16x16 transpose
  // passes used by the separate middle kernels.
  {
    GF31 u[MIDDLE];
    readMiddleOutLine(u, in31, tid, x);
    middleMul(u, x, middle31);
    fft_MIDDLE(u);
    middleMul2(u, tid, x, middle31);
    for (u32 m = 0; m != MIDDLE; ++m) tile31[m * WIDTH + tid] = u[m];
  }
  {
    GF61 u[MIDDLE];
    readMiddleOutLine(u, in61, tid, x);
    middleMul(u, x, middle61);
    fft_MIDDLE(u);
    middleMul2(u, tid, x, middle61);
    for (u32 m = 0; m != MIDDLE; ++m) tile61[m * WIDTH + tid] = u[m];
  }
  bar();

  for (u32 i = 0; i != NW; ++i) {
    u32 const y = lowMe + i * G_W;
    u31[i] = tile31[middle * WIDTH + y];
  }
#endif
  fft_WIDTH1(tile31, u31, width31, MIDDLE, lowMe);
  bar();
  for (u32 i = 0; i != NW; ++i) {
    tile31[middle * WIDTH + lowMe + i * G_W] = u31[i];
  }
  bar();
#if !CM_SKIP_MIDDLE
  for (u32 i = 0; i != NW; ++i) {
    u32 const y = lowMe + i * G_W;
    u61[i] = tile61[middle * WIDTH + y];
  }
#endif
  fft_WIDTH1(tile61, u61, width61, MIDDLE, lowMe);
  bar();
  for (u32 i = 0; i != NW; ++i) {
    tile61[middle * WIDTH + lowMe + i * G_W] = u61[i];
  }
  bar();

  const u32 m31_log2_root_two = (u32)(((1ULL << 30) / NWORDS) % 31);
  const u32 m31_bigword_weight_shift =
    (NWORDS - EXP % NWORDS) * m31_log2_root_two % 31;
  const u32 m31_shift_step = (m31_bigword_weight_shift + 30) % 31;
  const u32 m61_log2_root_two = (u32)(((1ULL << 60) / NWORDS) % 61);
  const u32 m61_bigword_weight_shift =
    (NWORDS - EXP % NWORDS) * m61_log2_root_two % 61;
  const u32 m61_shift_step = (m61_bigword_weight_shift + 60) % 61;
  const u64 m31_combo_step = make_u64(m31_shift_step, FRAC_BPW_HI);
  const u64 m61_combo_step = make_u64(m61_shift_step, FRAC_BPW_HI);
  const u64 m31_combo_bigstep =
    (comboFracBits(G_W * BIG_HEIGHT * 2 - 1) +
     make_u64((G_W * BIG_HEIGHT * 2 - 1) * m31_shift_step, 0)) %
    (31ULL << 32);
  const u64 m61_combo_bigstep =
    (comboFracBits(G_W * BIG_HEIGHT * 2 - 1) +
     make_u64((G_W * BIG_HEIGHT * 2 - 1) * m61_shift_step, 0)) %
    (61ULL << 32);

  union { uint2 a; u64 b; } c31, c61;
  u64 start31, start61;
  cmInitCounters((lowMe * BIG_HEIGHT + line) * 2,
                 &c31.b, &c61.b, &start31, &start61);

  u32 roundMax = 0;
  float carryMax = 0;
  for (u32 i = 0; i != NW; ++i) {
    u32 const shift310 = c31.a[1];
    c31.b += m31_combo_step;
    c31.a[1] = adjust_m31_weight_shift(c31.a[1]);
    u32 const shift311 = c31.a[1];
    u32 const shift610 = c61.a[1];
    c61.b += m61_combo_step;
    c61.a[1] = adjust_m61_weight_shift(c61.a[1]);
    u32 const shift611 = c61.a[1];
    bool const biglit0 = c31.a[0] <= FRAC_BPW_HI;
    bool const biglit1 = c31.a[0] >= -FRAC_BPW_HI;

    i64 carry;
    Word2 const wu = weightAndCarryPairSloppy(
      SWAP_XY(u31[i]), SWAP_XY(u61[i]),
      shift310, shift311, shift610, shift611,
      LL != 0,
      (LL & (i == 0) & (line == 0) & (lowMe == 0)) ? -2 : 0,
      biglit0, biglit1, &carry, &roundMax, &carryMax);

    u32 const y = lowMe + i * G_W;
    u32 const pair = line * WIDTH + y;
    packed[2 * pair] = (u32)wu.x;
    packed[2 * pair + 1] = (u32)wu.y;
    carryOut[pair] = carry;

    // At this bpw the unsigned first word and signed second word each need
    // exactly one bit beyond their low u32.  Store those bits as two ballot
    // planes: 8 bytes/value plus 2 bits/value in total.
    u32 const xMask = cmBallot(((u64)wu.x >> 32) != 0);
    u32 const yMask = cmBallot((((u64)wu.y >> 32) & 1) != 0);
    if ((tid & 31) == 0) {
      packed[NWORDS + pair / 32] = xMask;
      packed[NWORDS + ND / 32 + pair / 32] = yMask;
    }

    c31.b += m31_combo_bigstep;
    c31.a[1] = adjust_m31_weight_shift(c31.a[1]);
    c61.b += m61_combo_bigstep;
    c61.a[1] = adjust_m61_weight_shift(c61.a[1]);
  }

#if ROE
  updateStats((local u32 *)sharedBytes, CM_THREADS, SMALL_HEIGHT,
              bufROE, posROE,
              (float)roundMax / (float)0x1FFFFFFF);
#elif STATS & (1 << MUL3)
  updateStats((local u32 *)sharedBytes, CM_THREADS, SMALL_HEIGHT,
              bufROE, posROE, carryMax);
#endif
}

// Complete carry propagation after the first kernel's grid-wide barrier,
// reconstruct both weighted fields, and perform the forward width/middle
// edges directly into the tail-ready layout.
KERNEL(CM_THREADS) carryMiddleIn(P(T2) out, CP(u32) packed,
                                 CP(i64) carryOut,
                                 Trig middleTrig, Trig widthTrig) {
  CM_SHARED(sharedBytes);
  local GF31 *tile31 = (local GF31 *)sharedBytes;
  local GF61 *tile61 = (local GF61 *)(sharedBytes +
                                      MIDDLE * WIDTH * sizeof(GF31));

  P(GF31) out31 = (P(GF31))(out + DISTGF31);
  P(GF61) out61 = (P(GF61))(out + DISTGF61);
  TrigGF31 middle31 = (TrigGF31)(middleTrig + DISTMTRIGGF31);
  TrigGF61 middle61 = (TrigGF61)(middleTrig + DISTMTRIGGF61);
  TrigGF31 width31 = (TrigGF31)(widthTrig + DISTWTRIGGF31);
  TrigGF61 width61 = (TrigGF61)(widthTrig + DISTWTRIGGF61);

  u32 const x = get_group_id(0);
  u32 const tid = get_local_id(0);
  u32 const middle = tid / G_W;
  u32 const lowMe = tid % G_W;
  u32 const line = middle * SMALL_HEIGHT + x;

  const u32 m31_log2_root_two = (u32)(((1ULL << 30) / NWORDS) % 31);
  const u32 m31_bigword_weight_shift =
    (NWORDS - EXP % NWORDS) * m31_log2_root_two % 31;
  const u32 m31_shift_step = (m31_bigword_weight_shift + 30) % 31;
  const u32 m61_log2_root_two = (u32)(((1ULL << 60) / NWORDS) % 61);
  const u32 m61_bigword_weight_shift =
    (NWORDS - EXP % NWORDS) * m61_log2_root_two % 61;
  const u32 m61_shift_step = (m61_bigword_weight_shift + 60) % 61;
  const u64 m31_combo_step = make_u64(m31_shift_step, FRAC_BPW_HI);
  const u64 m61_combo_step = make_u64(m61_shift_step, FRAC_BPW_HI);
  const u64 m31_combo_bigstep =
    (comboFracBits(G_W * BIG_HEIGHT * 2 - 1) +
     make_u64((G_W * BIG_HEIGHT * 2 - 1) * m31_shift_step, 0)) %
    (31ULL << 32);
  const u64 m61_combo_bigstep =
    (comboFracBits(G_W * BIG_HEIGHT * 2 - 1) +
     make_u64((G_W * BIG_HEIGHT * 2 - 1) * m61_shift_step, 0)) %
    (61ULL << 32);

  union { uint2 a; u64 b; } c31, c61;
  u64 start31, start61;
  cmInitCounters((lowMe * BIG_HEIGHT + line) * 2,
                 &c31.b, &c61.b, &start31, &start61);
  // The inverse CRT uses the shifts adjusted for the 2*NWORDS transform
  // scale.  Forward weighting starts from the original, unadjusted counters,
  // exactly as carryFused restores its saved starting_combo_counter values.
  c31.b = start31;
  c61.b = start61;

  GF31 u31[NW];
  GF61 u61[NW];
  for (u32 i = 0; i != NW; ++i) {
    u32 const shift310 = c31.a[1];
    c31.b += m31_combo_step;
    c31.a[1] = adjust_m31_weight_shift(c31.a[1]);
    u32 const shift311 = c31.a[1];
    u32 const shift610 = c61.a[1];
    c61.b += m61_combo_step;
    c61.a[1] = adjust_m61_weight_shift(c61.a[1]);
    u32 const shift611 = c61.a[1];

    u32 const y = lowMe + i * G_W;
    u32 const pair = line * WIDTH + y;
    u32 const bit = pair & 31;
    u32 const xMask = packed[NWORDS + pair / 32];
    u32 const yMask = packed[NWORDS + ND / 32 + pair / 32];
    u64 const rawX = (u64)packed[2 * pair] |
                     ((u64)((xMask >> bit) & 1) << 32);
    u64 const rawY = (u64)packed[2 * pair + 1] |
                     ((u64)((yMask >> bit) & 1) << 32);
    Word2 wu = U2((i64)rawX, (i64)(rawY << 31) >> 31);

    i64 const inCarry = line != 0 ?
      carryOut[(line - 1) * WIDTH + y] :
      carryOut[(BIG_HEIGHT - 1) * WIDTH + (y + WIDTH - 1) % WIDTH];
    wu = carryFinal(wu, inCarry, c31.a[0] <= FRAC_BPW_HI);
    u31[i] = U2(shl(make_Z31(wu.x), shift310),
                shl(make_Z31(wu.y), shift311));
    u61[i] = U2(shl(make_Z61(wu.x), shift610),
                shl(make_Z61(wu.y), shift611));

    c31.b += m31_combo_bigstep;
    c31.a[1] = adjust_m31_weight_shift(c31.a[1]);
    c61.b += m61_combo_bigstep;
    c61.a[1] = adjust_m61_weight_shift(c61.a[1]);
  }

  fft_WIDTH2(tile31, u31, width31, MIDDLE, lowMe);
#if CM_SKIP_MIDDLE
  writeCarryFusedLine(u31, out31, line, lowMe);
#else
  bar();
  for (u32 i = 0; i != NW; ++i) {
    tile31[middle * WIDTH + lowMe + i * G_W] = u31[i];
  }
#endif
  bar();
  fft_WIDTH2(tile61, u61, width61, MIDDLE, lowMe);
#if CM_SKIP_MIDDLE
  writeCarryFusedLine(u61, out61, line, lowMe);
#else
  bar();
  for (u32 i = 0; i != NW; ++i) {
    tile61[middle * WIDTH + lowMe + i * G_W] = u61[i];
  }
  bar();

  // Reinterpret the same 512 threads as WIDTH coordinates and gather the
  // eight middle values.  Direct logical stores make the transpose-only LDS
  // stages of fftMiddleIn unnecessary.
  {
    GF31 u[MIDDLE];
    for (u32 m = 0; m != MIDDLE; ++m) u[m] = tile31[m * WIDTH + tid];
    middleMul2(u, tid, x, middle31);
    fft_MIDDLE(u);
    middleMul(u, x, middle31);
    // writeMiddleInLine's physical layout expects the value transpose that
    // the standalone 16x16 kernel performs immediately before its store.
    u32 const storeY = (x / 16) * 16 + tid % 16;
    u32 const storeX = (tid / 16) * 16 + x % 16;
    writeMiddleInLine(out31, u, storeY, storeX);
  }
  {
    GF61 u[MIDDLE];
    for (u32 m = 0; m != MIDDLE; ++m) u[m] = tile61[m * WIDTH + tid];
    middleMul2(u, tid, x, middle61);
    fft_MIDDLE(u);
    middleMul(u, x, middle61);
    u32 const storeY = (x / 16) * 16 + tid % 16;
    u32 const storeX = (tid / 16) * 16 + x % 16;
    writeMiddleInLine(out61, u, storeY, storeX);
  }
#endif
}

#undef CM_SHARED
#undef CM_SHARED_BYTES
#undef CM_THREADS

#endif
