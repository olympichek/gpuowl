// Copyright (C) Mihai Preda

#include "base.cl"
#include "math.cl"
#include "trig.cl"
#include "carryutil.cl"
#include "weight.cl"

#if FFT_TYPE == FFT64

// Carry propagation with optional MUL-3, over CARRY_LEN words.
// Input arrives with real and imaginary values swapped and weighted.

KERNEL(G_W) carry(P(Word2) out, CP(T2) in, u32 posROE, P(CarryABM) carryOut, BigTab THREAD_WEIGHTS, P(uint) bufROE) {
  u32 g  = get_group_id(0);
  u32 me = get_local_id(0);
  u32 gx = g % NW;
  u32 gy = g / NW;
  u32 H = BIG_HEIGHT;

  // & vs. && to workaround spurious warning
  CarryABM carry = (LL & (me == 0) & (g == 0)) ? -2 : 0;
  float roundMax = 0;
  float carryMax = 0;

  // Calculate the most significant 32-bits of FRAC_BPW * the index of the FFT word.  Also add FRAC_BPW_HI to test first biglit flag.
  u32 line = gy * CARRY_LEN;
  u32 word_index = (gx * G_W * H + me * H + line) * 2;
  u32 frac_bits = fracBits(word_index) + FRAC_BPW_HI;

  T base = optionalDouble(fancyMul(THREAD_WEIGHTS[me].x, iweightStep(gx)));

  for (i32 i = 0; i < CARRY_LEN; ++i) {
    u32 p = G_W * gx + WIDTH * (CARRY_LEN * gy + i) + me;
    T w1 = optionalDouble(fancyMul(base, THREAD_WEIGHTS[G_W + gy * CARRY_LEN + i].x));
    T w2 = optionalDouble(fancyMul(w1, IWEIGHT_STEP));
    bool biglit0 = frac_bits + (2*i) * FRAC_BPW_HI <= FRAC_BPW_HI;
    bool biglit1 = frac_bits + (2*i) * FRAC_BPW_HI >= -FRAC_BPW_HI;   // Same as frac_bits + (2*i) * FRAC_BPW_HI + FRAC_BPW_HI <= FRAC_BPW_HI;
    out[p] = weightAndCarryPair(SWAP_XY(in[p]), U2(w1, w2), carry, biglit0, biglit1, &carry, &roundMax, &carryMax);
  }
  carryOut[G_W * g + me] = carry;

#if ROE
  local u32 lds[G_W];
  updateStats(lds, G_W, H, bufROE, posROE, roundMax);
#elif (STATS & (1 << (2 + MUL3)))
  local u32 lds[G_W];
  updateStats(bufROE, posROE, carryMax);
#endif
}


/**************************************************************************/
/*            Similar to above, but for an FFT based on FP32              */
/**************************************************************************/

#elif FFT_TYPE == FFT32

// Carry propagation with optional MUL-3, over CARRY_LEN words.
// Input arrives with real and imaginary values swapped and weighted.

KERNEL(G_W) carry(P(Word2) out, CP(F2) in, u32 posROE, P(CarryABM) carryOut, BigTabFP32 THREAD_WEIGHTS, P(uint) bufROE) {
  u32 g  = get_group_id(0);
  u32 me = get_local_id(0);
  u32 gx = g % NW;
  u32 gy = g / NW;
  u32 H = BIG_HEIGHT;

  CarryABM carry = (LL & (me == 0) & (g == 0)) ? -2 : 0;
  float roundMax = 0;
  float carryMax = 0;

  // Calculate the most significant 32-bits of FRAC_BPW * the index of the FFT word.
  u32 line = gy * CARRY_LEN;
  u32 word_index = (gx * G_W * H + me * H + line) * 2;

  F base = fancyMul(THREAD_WEIGHTS[me].x, iweightStep(gx));
  u32 me_frac_bits = fracBits(me * H * 2);
  u32 step_frac_bits = weightStepFracBits(gx);
  u32 base_frac_bits = me_frac_bits + step_frac_bits;
  base = optionalDouble(base, base_frac_bits > step_frac_bits);

  u32 frac_bits = fracBits(word_index);

  // Base_frac_bits and frac_bits are inexact values.  We only want to trigger an optional double when it is clear to do so.
  // Fudge base_frac_bits to make it harder to trigger a double when the two inexact values are equal.
  base_frac_bits++;

  for (i32 i = 0; i < CARRY_LEN; ++i) {
    u32 p = G_W * gx + WIDTH * (line + i) + me;
    F w1 = optionalDouble(fancyMul(base, THREAD_WEIGHTS[G_W + line + i].x), frac_bits > base_frac_bits);
    F w2 = optionalDouble(fancyMul(w1, IWEIGHT_STEP), frac_bits + FRAC_BPW_HI > FRAC_BPW_HI);
    frac_bits += FRAC_BPW_HI;
    bool biglit0 = frac_bits <= FRAC_BPW_HI;
    bool biglit1 = frac_bits >= -FRAC_BPW_HI;   // Same as frac_bits + FRAC_BPW_HI <= FRAC_BPW_HI;
    out[p] = weightAndCarryPair(SWAP_XY(in[p]), U2(w1, w2), carry, biglit0, biglit1, &carry, &roundMax, &carryMax);
    // Generate frac_bits for next pair
    frac_bits += FRAC_BPW_HI;
  }
  carryOut[G_W * g + me] = carry;

#if ROE
  local u32 lds[G_W];
  updateStats(lds, G_W, H, bufROE, posROE, roundMax);
#elif (STATS & (1 << (2 + MUL3)))
  local u32 lds[G_W];
  updateStats(lds, G_W, H, bufROE, posROE, carryMax);
#endif
}


/**************************************************************************/
/*          Similar to above, but for an NTT based on GF(M31^2)           */
/**************************************************************************/

#elif FFT_TYPE == FFT31

KERNEL(G_W) carry(P(Word2) out, CP(GF31) in, u32 posROE, P(CarryABM) carryOut, P(uint) bufROE) {
  u32 g  = get_group_id(0);
  u32 me = get_local_id(0);
  u32 gx = g % NW;
  u32 gy = g / NW;
  u32 H = BIG_HEIGHT;
  u32 line = gy * CARRY_LEN;

  // & vs. && to workaround spurious warning
  CarryABM carry = (LL & (me == 0) & (g == 0)) ? -2 : 0;
  u32 roundMax = 0;
  float carryMax = 0;

  u32 word_index = (gx * G_W * H + me * H + line) * 2;

  // Weight is 2^[ceil(qj / n) - qj/n] where j is the word index, q is the Mersenne exponent, and n is the number of words.
  // Weights can be applied with shifts because 2 is the 30th root GF31.
  // Let s be the shift amount for word 1.  The shift amount for word x is ceil(x * (s - 1) + num_big_words_less_than_x) % 31.
  const u32 log2_root_two = (u32) (((1ULL << 30) / NWORDS) % 31);
  const u32 bigword_weight_shift = (NWORDS - EXP % NWORDS) * log2_root_two % 31;
  const u32 bigword_weight_shift_minus1 = (bigword_weight_shift + 30) % 31;

  // Derive the big vs. little flags from the fractional number of bits in each word.
  // Create a 64-bit counter that tracks both weight shifts and frac_bits (adding 0xFFFFFFFF to effect the ceil operation required for weight shift).
  union { uint2 a; u64 b; } combo;
#define frac_bits       combo.a[0]
#define weight_shift    combo.a[1]
#define combo_counter   combo.b

  const u64 combo_step = make_u64(bigword_weight_shift_minus1, FRAC_BPW_HI);
  combo_counter = comboFracBits(word_index) + make_u64(word_index * bigword_weight_shift_minus1, 0xFFFFFFFF);

  // We also adjust shift amount for the fact that NTT returns results multiplied by 2*NWORDS.
  const u32 log2_NWORDS = (WIDTH == 256 ? 8 : WIDTH == 512 ? 9 : WIDTH == 1024 ? 10 : 12) +
                          (MIDDLE == 1 ? 0 : MIDDLE == 2 ? 1 : MIDDLE == 4 ? 2 : MIDDLE == 8 ? 3 : 4) +
                          (SMALL_HEIGHT == 256 ? 8 : SMALL_HEIGHT == 512 ? 9 : SMALL_HEIGHT == 1024 ? 10 : 12) + 1;
  weight_shift = (weight_shift + log2_NWORDS + 1) % 31;

  for (i32 i = 0; i < CARRY_LEN; ++i) {
    u32 p = G_W * gx + WIDTH * (CARRY_LEN * gy + i) + me;

    // Generate the second weight shift
    u32 weight_shift0 = weight_shift;
    combo_counter += combo_step;
    if (weight_shift > 31) weight_shift -= 31;
    u32 weight_shift1 = weight_shift;

    // Generate big-word/little-word flags
    bool biglit0 = frac_bits <= FRAC_BPW_HI;
    bool biglit1 = frac_bits >= -FRAC_BPW_HI;   // Same as frac_bits + FRAC_BPW_HI <= FRAC_BPW_HI;

    // Compute result
    out[p] = weightAndCarryPair(SWAP_XY(in[p]), weight_shift0, weight_shift1, carry, biglit0, biglit1, &carry, &roundMax, &carryMax);
    // Generate weight shifts and frac_bits for next pair
    combo_counter += combo_step;
    if (weight_shift > 31) weight_shift -= 31;
  }
  carryOut[G_W * g + me] = carry;

#if ROE
  local u32 lds[G_W];
  float fltRoundMax = (float) roundMax / (float) M31;      // For speed, roundoff was computed as 32-bit integer.  Convert to float.
  updateStats(lds, G_W, H, bufROE, posROE, fltRoundMax);
#elif (STATS & (1 << (2 + MUL3)))
  local u32 lds[G_W];
  updateStats(lds, G_W, H, bufROE, posROE, carryMax);
#endif
}


/**************************************************************************/
/*          Similar to above, but for an NTT based on GF(M61^2)           */
/**************************************************************************/

#elif FFT_TYPE == FFT61

KERNEL(G_W) carry(P(Word2) out, CP(GF61) in, u32 posROE, P(CarryABM) carryOut, P(uint) bufROE) {
  u32 g  = get_group_id(0);
  u32 me = get_local_id(0);
  u32 gx = g % NW;
  u32 gy = g / NW;
  u32 H = BIG_HEIGHT;
  u32 line = gy * CARRY_LEN;

  // & vs. && to workaround spurious warning
  CarryABM carry = (LL & (me == 0) & (g == 0)) ? -2 : 0;
  u32 roundMax = 0;
  float carryMax = 0;

  u32 word_index = (gx * G_W * H + me * H + line) * 2;

  // Weight is 2^[ceil(qj / n) - qj/n] where j is the word index, q is the Mersenne exponent, and n is the number of words.
  // Weights can be applied with shifts because 2 is the 60th root GF61.
  // Let s be the shift amount for word 1.  The shift amount for word x is ceil(x * (s - 1) + num_big_words_less_than_x) % 61.
  const u32 log2_root_two = (u32) (((1ULL << 60) / NWORDS) % 61);
  const u32 bigword_weight_shift = (NWORDS - EXP % NWORDS) * log2_root_two % 61;
  const u32 bigword_weight_shift_minus1 = (bigword_weight_shift + 60) % 61;

  // Derive the big vs. little flags from the fractional number of bits in each word.
  // Create a 64-bit counter that tracks both weight shifts and frac_bits (adding 0xFFFFFFFF to effect the ceil operation required for weight shift).
  union { uint2 a; u64 b; } combo;
#define frac_bits       combo.a[0]
#define weight_shift    combo.a[1]
#define combo_counter   combo.b

  const u64 combo_step = make_u64(bigword_weight_shift_minus1, FRAC_BPW_HI);
  combo_counter = comboFracBits(word_index) + make_u64(word_index * bigword_weight_shift_minus1, 0xFFFFFFFF);

  // We also adjust shift amount for the fact that NTT returns results multiplied by 2*NWORDS.
  const u32 log2_NWORDS = (WIDTH == 256 ? 8 : WIDTH == 512 ? 9 : WIDTH == 1024 ? 10 : 12) +
                          (MIDDLE == 1 ? 0 : MIDDLE == 2 ? 1 : MIDDLE == 4 ? 2 : MIDDLE == 8 ? 3 : 4) +
                          (SMALL_HEIGHT == 256 ? 8 : SMALL_HEIGHT == 512 ? 9 : SMALL_HEIGHT == 1024 ? 10 : 12) + 1;
  weight_shift = (weight_shift + log2_NWORDS + 1) % 61;

  for (i32 i = 0; i < CARRY_LEN; ++i) {
    u32 p = G_W * gx + WIDTH * (CARRY_LEN * gy + i) + me;

    // Generate the second weight shift
    u32 weight_shift0 = weight_shift;
    combo_counter += combo_step;
    if (weight_shift > 61) weight_shift -= 61;
    u32 weight_shift1 = weight_shift;

    // Generate big-word/little-word flags
    bool biglit0 = frac_bits <= FRAC_BPW_HI;
    bool biglit1 = frac_bits >= -FRAC_BPW_HI;   // Same as frac_bits + FRAC_BPW_HI <= FRAC_BPW_HI;

    // Compute result
    out[p] = weightAndCarryPair(SWAP_XY(in[p]), weight_shift0, weight_shift1, carry, biglit0, biglit1, &carry, &roundMax, &carryMax);
    // Generate weight shifts and frac_bits for next pair
    combo_counter += combo_step;
    if (weight_shift > 61) weight_shift -= 61;
  }
  carryOut[G_W * g + me] = carry;

#if ROE
  local u32 lds[G_W];
  float fltRoundMax = (float) roundMax / (float) (M61 >> 32);      // For speed, roundoff was computed as 32-bit integer.  Convert to float.
  updateStats(lds, G_W, H, bufROE, posROE, fltRoundMax);
#elif (STATS & (1 << (2 + MUL3)))
  local u32 lds[G_W];
  updateStats(lds, G_W, H, bufROE, posROE, carryMax);
#endif
}


/**************************************************************************/
/*    Similar to above, but for a hybrid FFT based on FP64 & GF(M31^2)    */
/**************************************************************************/

#elif FFT_TYPE == FFT6431

KERNEL(G_W) carry(P(Word2) out, CP(T2) in, u32 posROE, P(CarryABM) carryOut, BigTab THREAD_WEIGHTS, P(uint) bufROE) {
  u32 g  = get_group_id(0);
  u32 me = get_local_id(0);
  u32 gx = g % NW;
  u32 gy = g / NW;
  u32 H = BIG_HEIGHT;
  u32 line = gy * CARRY_LEN;

  CP(GF31) in31 = (CP(GF31)) (in + DISTGF31);

  CarryABM carry = (LL & (me == 0) & (g == 0)) ? -2 : 0;
  float roundMax = 0;
  float carryMax = 0;

  u32 word_index = (gx * G_W * H + me * H + line) * 2;

  T base = optionalDouble(fancyMul(THREAD_WEIGHTS[me].x, iweightStep(gx)));

  // Weight is 2^[ceil(qj / n) - qj/n] where j is the word index, q is the Mersenne exponent, and n is the number of words.
  const u32 log2_root_two = (u32) (((1ULL << 30) / NWORDS) % 31);
  const u32 bigword_weight_shift = (NWORDS - EXP % NWORDS) * log2_root_two % 31;
  const u32 bigword_weight_shift_minus1 = (bigword_weight_shift + 30) % 31;

  // Derive the big vs. little flags from the fractional number of bits in each word.
  // Create a 64-bit counter that tracks both weight shifts and frac_bits (adding 0xFFFFFFFF to effect the ceil operation required for weight shift).
  union { uint2 a; u64 b; } combo;
#define frac_bits       combo.a[0]
#define weight_shift    combo.a[1]
#define combo_counter   combo.b

  const u64 combo_step = make_u64(bigword_weight_shift_minus1, FRAC_BPW_HI);
  combo_counter = comboFracBits(word_index) + make_u64(word_index * bigword_weight_shift_minus1, 0xFFFFFFFF);

  // We also adjust shift amount for the fact that NTT returns results multiplied by 2*NWORDS.
  const u32 log2_NWORDS = (WIDTH == 256 ? 8 : WIDTH == 512 ? 9 : WIDTH == 1024 ? 10 : 12) +
                          (MIDDLE == 1 ? 0 : MIDDLE == 2 ? 1 : MIDDLE == 4 ? 2 : MIDDLE == 8 ? 3 : 4) +
                          (SMALL_HEIGHT == 256 ? 8 : SMALL_HEIGHT == 512 ? 9 : SMALL_HEIGHT == 1024 ? 10 : 12) + 1;
  weight_shift = (weight_shift + log2_NWORDS + 1) % 31;

  for (i32 i = 0; i < CARRY_LEN; ++i) {
    u32 p = G_W * gx + WIDTH * (CARRY_LEN * gy + i) + me;

    // Generate the FP64 and second GF31 weight shift
    T w1 = optionalDouble(fancyMul(base, THREAD_WEIGHTS[G_W + gy * CARRY_LEN + i].x));
    T w2 = optionalDouble(fancyMul(w1, IWEIGHT_STEP));
    u32 weight_shift0 = weight_shift;
    combo_counter += combo_step;
    if (weight_shift > 31) weight_shift -= 31;
    u32 weight_shift1 = weight_shift;

    // Generate big-word/little-word flags
    bool biglit0 = frac_bits <= FRAC_BPW_HI;
    bool biglit1 = frac_bits >= -FRAC_BPW_HI;   // Same as frac_bits + FRAC_BPW_HI <= FRAC_BPW_HI;

    // Compute result
    out[p] = weightAndCarryPair(SWAP_XY(in[p]), SWAP_XY(in31[p]), w1, w2, weight_shift0, weight_shift1,
                                LL != 0 || i != 0, carry, biglit0, biglit1, &carry, &roundMax, &carryMax);

    // Generate weight shifts and frac_bits for next pair
    combo_counter += combo_step;
    if (weight_shift > 31) weight_shift -= 31;
  }
  carryOut[G_W * g + me] = carry;

#if ROE
  local u32 lds[G_W];
  updateStats(lds, G_W, H, bufROE, posROE, roundMax);
#elif (STATS & (1 << (2 + MUL3)))
  local u32 lds[G_W];
  updateStats(lds, G_W, H, bufROE, posROE, carryMax);
#endif
}


/**************************************************************************/
/*    Similar to above, but for a hybrid FFT based on FP32 & GF(M31^2)    */
/**************************************************************************/

#elif FFT_TYPE == FFT3231

KERNEL(G_W) carry(P(Word2) out, CP(T2) in, u32 posROE, P(CarryABM) carryOut, BigTabFP32 THREAD_WEIGHTS, P(uint) bufROE) {
  u32 g  = get_group_id(0);
  u32 me = get_local_id(0);
  u32 gx = g % NW;
  u32 gy = g / NW;
  u32 H = BIG_HEIGHT;
  u32 line = gy * CARRY_LEN;

  CP(F2) inF2 = (CP(F2)) in;
  CP(GF31) in31 = (CP(GF31)) (in + DISTGF31);

  CarryABM carry = (LL & (me == 0) & (g == 0)) ? -2 : 0;
  float roundMax = 0;
  float carryMax = 0;

  u32 word_index = (gx * G_W * H + me * H + line) * 2;

  F base = fancyMul(THREAD_WEIGHTS[me].x, iweightStep(gx));
  u32 me_frac_bits = fracBits(me * H * 2);
  u32 step_frac_bits = weightStepFracBits(gx);
  u32 base_frac_bits = me_frac_bits + step_frac_bits;
  base = optionalDouble(base, base_frac_bits > step_frac_bits);

  // Weight is 2^[ceil(qj / n) - qj/n] where j is the word index, q is the Mersenne exponent, and n is the number of words.
  const u32 log2_root_two = (u32) (((1ULL << 30) / NWORDS) % 31);
  const u32 bigword_weight_shift = (NWORDS - EXP % NWORDS) * log2_root_two % 31;
  const u32 bigword_weight_shift_minus1 = (bigword_weight_shift + 30) % 31;

  // Derive the big vs. little flags from the fractional number of bits in each word.
  // Create a 64-bit counter that tracks both weight shifts and frac_bits (adding 0xFFFFFFFF to effect the ceil operation required for weight shift).
  union { uint2 a; u64 b; } combo;
#define frac_bits       combo.a[0]
#define weight_shift    combo.a[1]
#define combo_counter   combo.b

  const u64 combo_step = make_u64(bigword_weight_shift_minus1, FRAC_BPW_HI);
  combo_counter = comboFracBits(word_index) + make_u64(word_index * bigword_weight_shift_minus1, 0xFFFFFFFF);

  // We also adjust shift amount for the fact that NTT returns results multiplied by 2*NWORDS.
  const u32 log2_NWORDS = (WIDTH == 256 ? 8 : WIDTH == 512 ? 9 : WIDTH == 1024 ? 10 : 12) +
                          (MIDDLE == 1 ? 0 : MIDDLE == 2 ? 1 : MIDDLE == 4 ? 2 : MIDDLE == 8 ? 3 : 4) +
                          (SMALL_HEIGHT == 256 ? 8 : SMALL_HEIGHT == 512 ? 9 : SMALL_HEIGHT == 1024 ? 10 : 12) + 1;
  weight_shift = (weight_shift + log2_NWORDS + 1) % 31;

  for (i32 i = 0; i < CARRY_LEN; ++i) {
    u32 p = G_W * gx + WIDTH * (line + i) + me;

    // Generate the FP32 and second GF31 weight shift
    F w1 = optionalDouble(fancyMul(base, THREAD_WEIGHTS[G_W + line + i].x), frac_bits > base_frac_bits);
    F w2 = optionalDouble(fancyMul(w1, IWEIGHT_STEP), frac_bits + FRAC_BPW_HI > FRAC_BPW_HI);
    u32 weight_shift0 = weight_shift;
    combo_counter += combo_step;
    if (weight_shift > 31) weight_shift -= 31;
    u32 weight_shift1 = weight_shift;

    // Generate big-word/little-word flags
    bool biglit0 = frac_bits <= FRAC_BPW_HI;
    bool biglit1 = frac_bits >= -FRAC_BPW_HI;   // Same as frac_bits + FRAC_BPW_HI <= FRAC_BPW_HI;

    // Compute result
    out[p] = weightAndCarryPair(SWAP_XY(inF2[p]), SWAP_XY(in31[p]), w1, w2, weight_shift0, weight_shift1,
                                carry, biglit0, biglit1, &carry, &roundMax, &carryMax);

    // Generate weight shifts and frac_bits for next pair
    combo_counter += combo_step;
    if (weight_shift > 31) weight_shift -= 31;
  }
  carryOut[G_W * g + me] = carry;

#if ROE
  local u32 lds[G_W];
  updateStats(lds, G_W, H, bufROE, posROE, roundMax);
#elif (STATS & (1 << (2 + MUL3)))
  local u32 lds[G_W];
  updateStats(lds, G_W, H, bufROE, posROE, carryMax);
#endif
}


/**************************************************************************/
/*    Similar to above, but for a hybrid FFT based on FP32 & GF(M61^2)    */
/**************************************************************************/

#elif FFT_TYPE == FFT3261

KERNEL(G_W) carry(P(Word2) out, CP(T2) in, u32 posROE, P(CarryABM) carryOut, BigTabFP32 THREAD_WEIGHTS, P(uint) bufROE
#if PARITY_SQUARE
                    , CP(u32) parityIn
#endif
                    ) {
  u32 g  = get_group_id(0);
  u32 me = get_local_id(0);
  u32 gx = g % NW;
  u32 gy = g / NW;
  u32 H = BIG_HEIGHT;
  u32 line = gy * CARRY_LEN;

  CP(F2) inF2 = (CP(F2)) in;
  CP(GF61) in61 = (CP(GF61)) (in + DISTGF61);

  CarryABM carry = (LL & (me == 0) & (g == 0)) ? -2 : 0;
  float roundMax = 0;
  float carryMax = 0;

  u32 word_index = (gx * G_W * H + me * H + line) * 2;

  F base = fancyMul(THREAD_WEIGHTS[me].x, iweightStep(gx));
  u32 me_frac_bits = fracBits(me * H * 2);
  u32 step_frac_bits = weightStepFracBits(gx);
  u32 base_frac_bits = me_frac_bits + step_frac_bits;
  base = optionalDouble(base, base_frac_bits > step_frac_bits);

  // Weight is 2^[ceil(qj / n) - qj/n] where j is the word index, q is the Mersenne exponent, and n is the number of words.
#if GOOD_THOMAS9
  const u32 log2_root_two = 30;
#else
  const u32 log2_root_two = (u32) (((1ULL << 60) / NWORDS) % 61);
#endif
  const u32 bigword_weight_shift = (NWORDS - EXP % NWORDS) * log2_root_two % 61;
  const u32 bigword_weight_shift_minus1 = (bigword_weight_shift + 60) % 61;

  // Derive the big vs. little flags from the fractional number of bits in each word.
  // Create a 64-bit counter that tracks both weight shifts and frac_bits (adding 0xFFFFFFFF to effect the ceil operation required for weight shift).
  union { uint2 a; u64 b; } combo;
#define frac_bits       combo.a[0]
#define weight_shift    combo.a[1]
#define combo_counter   combo.b

  const u64 combo_step = make_u64(bigword_weight_shift_minus1, FRAC_BPW_HI);
  combo_counter = comboFracBits(word_index) + make_u64(word_index * bigword_weight_shift_minus1, 0xFFFFFFFF);

  // We also adjust shift amount for the fact that NTT returns results multiplied by 2*NWORDS.
#if GOOD_THOMAS9
  // The inverse radix-nine normalization is applied in weightAndCarryOne;
  // only the 2^19-word binary channel scale remains here.
  const u32 log2_NWORDS = 19;
#else
  const u32 log2_NWORDS = (WIDTH == 256 ? 8 : WIDTH == 512 ? 9 : WIDTH == 1024 ? 10 : 12) +
                          (MIDDLE == 1 ? 0 : MIDDLE == 2 ? 1 : MIDDLE == 4 ? 2 : MIDDLE == 8 ? 3 : 4) +
                          (SMALL_HEIGHT == 256 ? 8 : SMALL_HEIGHT == 512 ? 9 : SMALL_HEIGHT == 1024 ? 10 : 12) + 1;
#endif
  weight_shift = (weight_shift + log2_NWORDS + 1) % 61;

  for (i32 i = 0; i < CARRY_LEN; ++i) {
    u32 p = G_W * gx + WIDTH * (line + i) + me;

    // Generate the FP32 and second GF61 weight shift
    F w1 = optionalDouble(fancyMul(base, THREAD_WEIGHTS[G_W + line + i].x), frac_bits > base_frac_bits);
    F w2 = optionalDouble(fancyMul(w1, IWEIGHT_STEP), frac_bits + FRAC_BPW_HI > FRAC_BPW_HI);
    u32 weight_shift0 = weight_shift;
    combo_counter += combo_step;
    if (weight_shift > 61) weight_shift -= 61;
    u32 weight_shift1 = weight_shift;

    // Generate big-word/little-word flags
    bool biglit0 = frac_bits <= FRAC_BPW_HI;
    bool biglit1 = frac_bits >= -FRAC_BPW_HI;   // Same as frac_bits + FRAC_BPW_HI <= FRAC_BPW_HI;

    // Compute result
#if GOOD_THOMAS9
    // Convert the natural output pair to its two scalar Good--Thomas planes.
    // fftW has already completed independently in every channel plane.
    const u32 gt_line = CARRY_LEN * gy + i;
    const u32 gt_m = gt_line / SMALL_HEIGHT;
    const u32 gt_y = gt_line % SMALL_HEIGHT;
    const u32 gt_x = G_W * gx + me;
    const u32 gt_x_channel = (9 * gt_x + gt_m) & (WIDTH - 1);
    const u32 gt_channel_even = (2 * (gt_y + 9 - gt_m)) % 9;
    const u32 gt_channel_odd = (gt_channel_even + 1) % 9;
    const u32 gt_even_storage =
      (gt_channel_even * SMALL_HEIGHT + gt_y) * WIDTH + gt_x_channel;
    const u32 gt_odd_storage =
      (gt_channel_odd * SMALL_HEIGHT + gt_y) * WIDTH + gt_x_channel;
    const GF61 gt_u61 = U2(in61[gt_odd_storage].x,
                           in61[gt_even_storage].y);
#else
    const GF61 gt_u61 = in61[p];
#endif
    out[p] = weightAndCarryPair(SWAP_XY(inF2[p]), SWAP_XY(gt_u61), w1, w2, weight_shift0, weight_shift1,
#if PARITY_SQUARE
#if PARITY_LAZY
#if PARITY_PHYSICAL
                                parityIn, gx * G_W * H + me * H + line + i,
#else
                                parityIn, p,
#endif
#else
#if PARITY_PREPARED
                                squareCoefficientParity(parityIn, p), 0,
#else
                                squareCoefficientParity(parityIn, gx * G_W * H + me * H + line + i), 0,
#endif
#endif
#endif
                                LL != 0 || i != 0, carry, biglit0, biglit1, &carry, &roundMax, &carryMax);

    // Generate weight shifts and frac_bits for next pair
    combo_counter += combo_step;
    if (weight_shift > 61) weight_shift -= 61;
  }
  carryOut[G_W * g + me] = carry;

#if ROE
  local u32 lds[G_W];
  updateStats(lds, G_W, H, bufROE, posROE, roundMax);
#elif (STATS & (1 << (2 + MUL3)))
  local u32 lds[G_W];
  updateStats(lds, G_W, H, bufROE, posROE, carryMax);
#endif
}


/**************************************************************************/
/*    Similar to above, but for an NTT based on GF(M31^2)*GF(M61^2)       */
/**************************************************************************/

#elif FFT_TYPE == FFT3161

KERNEL(G_W) carry(P(Word2) out, CP(T2) in, u32 posROE, P(CarryABM) carryOut
#if GOLD_PAIR
                   , BigTab THREAD_WEIGHTS
#endif
                   , P(uint) bufROE) {
  u32 g  = get_group_id(0);
  u32 me = get_local_id(0);
  u32 gx = g % NW;
  u32 gy = g / NW;
  u32 H = BIG_HEIGHT;
  u32 line = gy * CARRY_LEN;

  CP(GF31) in31 = (CP(GF31)) (in + DISTGF31);
  CP(GF61) in61 = (CP(GF61)) (in + DISTGF61);

  // & vs. && to workaround spurious warning
  CarryABM carry = (LL & (me == 0) & (g == 0)) ? -2 : 0;
  u32 roundMax = 0;
  float carryMax = 0;

  u32 word_index = (gx * G_W * H + me * H + line) * 2;

  // Weight is 2^[ceil(qj / n) - qj/n] where j is the word index, q is the Mersenne exponent, and n is the number of words.
  const u32 m31_log2_root_two = (u32) (((1ULL << 30) / NWORDS) % 31);
  const u32 m31_bigword_weight_shift = (NWORDS - EXP % NWORDS) * m31_log2_root_two % 31;
  const u32 m31_bigword_weight_shift_minus1 = (m31_bigword_weight_shift + 30) % 31;
#if !GOLD_PAIR
  const u32 m61_log2_root_two = (u32) (((1ULL << 60) / NWORDS) % 61);
  const u32 m61_bigword_weight_shift = (NWORDS - EXP % NWORDS) * m61_log2_root_two % 61;
  const u32 m61_bigword_weight_shift_minus1 = (m61_bigword_weight_shift + 60) % 61;
#endif

  // Derive the big vs. little flags from the fractional number of bits in each word.
  // Create a 64-bit counter that tracks both weight shifts and frac_bits (adding 0xFFFFFFFF to effect the ceil operation required for weight shift).
  union { uint2 a; u64 b; } m31_combo;
#define frac_bits           m31_combo.a[0]
#define m31_weight_shift    m31_combo.a[1]
#define m31_combo_counter   m31_combo.b
#if !GOLD_PAIR
  union { uint2 a; u64 b; } m61_combo;
#define m61_weight_shift    m61_combo.a[1]
#define m61_combo_counter   m61_combo.b
#endif

  const u64 m31_combo_step = make_u64(m31_bigword_weight_shift_minus1, FRAC_BPW_HI);
  m31_combo_counter = comboFracBits(word_index) + make_u64(word_index * m31_bigword_weight_shift_minus1, 0xFFFFFFFF);
#if !GOLD_PAIR
  const u64 m61_combo_step = make_u64(m61_bigword_weight_shift_minus1, FRAC_BPW_HI);
  m61_combo_counter = comboFracBits(word_index) + make_u64(word_index * m61_bigword_weight_shift_minus1, 0xFFFFFFFF);
#endif

  // We also adjust shift amount for the fact that NTT returns results multiplied by 2*NWORDS.
  const u32 log2_NWORDS = (WIDTH == 256 ? 8 : WIDTH == 512 ? 9 : WIDTH == 1024 ? 10 : 12) +
                          (MIDDLE == 1 ? 0 : MIDDLE == 2 ? 1 : MIDDLE == 4 ? 2 : MIDDLE == 8 ? 3 : 4) +
                          (SMALL_HEIGHT == 256 ? 8 : SMALL_HEIGHT == 512 ? 9 : SMALL_HEIGHT == 1024 ? 10 : 12) + 1;
  m31_weight_shift = (m31_weight_shift + log2_NWORDS + 1) % 31;
#if !GOLD_PAIR
  m61_weight_shift = (m61_weight_shift + log2_NWORDS + 1) % 61;
#endif

  for (i32 i = 0; i < CARRY_LEN; ++i) {
    u32 p = G_W * gx + WIDTH * (CARRY_LEN * gy + i) + me;

    // Generate the second weight shifts
    u32 m31_weight_shift0 = m31_weight_shift;
    m31_combo_counter += m31_combo_step;
    m31_weight_shift = adjust_m31_weight_shift(m31_weight_shift);
    u32 m31_weight_shift1 = m31_weight_shift;
#if GOLD_PAIR
    u32 goldExponent;
    Z61 goldInvWeight0 = goldStartingWeight(
      THREAD_WEIGHTS, gx * G_W + me, line + i, true, &goldExponent);
    goldInvWeight0 = mul(goldInvWeight0, (Z61)GOLD_INV_ND);
    Z61 const goldInvWeight1 = advanceGoldInverse(
      goldInvWeight0, GOLD_DELTA_ONE, (Z61)GOLD_INV_ONE,
      &goldExponent);
#else
    u32 m61_weight_shift0 = m61_weight_shift;
    m61_combo_counter += m61_combo_step;
    m61_weight_shift = adjust_m61_weight_shift(m61_weight_shift);
    u32 m61_weight_shift1 = m61_weight_shift;
#endif

    // Generate big-word/little-word flags
    bool biglit0 = frac_bits <= FRAC_BPW_HI;
    bool biglit1 = frac_bits >= -FRAC_BPW_HI;   // Same as frac_bits + FRAC_BPW_HI <= FRAC_BPW_HI;

    // Compute result
    out[p] = weightAndCarryPair(SWAP_XY(in31[p]),
#if GOLD_PAIR
                                in61[p],
#else
                                SWAP_XY(in61[p]),
#endif
                                m31_weight_shift0, m31_weight_shift1,
#if GOLD_PAIR
                                goldInvWeight0, goldInvWeight1,
#else
                                m61_weight_shift0, m61_weight_shift1,
#endif
                                LL != 0 || i != 0, carry, biglit0, biglit1, &carry, &roundMax, &carryMax);

    // Generate weight shifts and frac_bits for next pair
    m31_combo_counter += m31_combo_step;
    m31_weight_shift = adjust_m31_weight_shift(m31_weight_shift);
#if !GOLD_PAIR
    m61_combo_counter += m61_combo_step;
    m61_weight_shift = adjust_m61_weight_shift(m61_weight_shift);
#endif
  }
  carryOut[G_W * g + me] = carry;

#if ROE
  local u32 lds[G_W];
  float fltRoundMax = (float) roundMax / (float) 0x1FFFFFFF;      // For speed, roundoff was computed as 32-bit integer.  Convert to float.
  updateStats(lds, G_W, H, bufROE, posROE, fltRoundMax);
#elif (STATS & (1 << (2 + MUL3)))
  local u32 lds[G_W];
  updateStats(lds, G_W, H, bufROE, posROE, carryMax);
#endif
}


/**************************************************************************/
/*        Three independent 31-bit residue planes (FFT31R2 / FFT54)      */
/**************************************************************************/

#elif FFT_TYPE == FFT31R2

KERNEL(G_W) carry(P(Word2) out, CP(T2) in, u32 posROE, P(CarryABM) carryOut,
                  BigTab THREAD_WEIGHTS, P(uint) bufROE) {
  u32 const g = get_group_id(0);
  u32 const me = get_local_id(0);
  u32 const gx = g % NW;
  u32 const gy = g / NW;
  u32 const H = BIG_HEIGHT;
  u32 const line = gy * CARRY_LEN;

  CP(GF31) in31 = (CP(GF31))(in + DISTGF31);
#if RIESEL_PAIR
  CP(GF61) inPair = (CP(GF61))(in + DISTGF61);
#else
  CP(GF31) in0 = (CP(GF31))(in + DISTR0);
  CP(GF31) in1 = (CP(GF31))(in + DISTR1);
#endif
  CarryABM carry = (LL & (me == 0) & (g == 0)) ? -2 : 0;
  u32 roundMax = 0;
  float carryMax = 0;

  u32 const word_index = (gx * G_W * H + me * H + line) * 2;
  const u32 log2_root_two = (u32)(((1ULL << 30) / NWORDS) % 31);
  const u32 bigword_weight_shift = (NWORDS - EXP % NWORDS) * log2_root_two % 31;
  const u32 bigword_weight_shift_minus1 = (bigword_weight_shift + 30) % 31;
  union { uint2 a; u64 b; } combo;
#define frac_bits       combo.a[0]
#define weight_shift    combo.a[1]
#define combo_counter   combo.b
  const u64 combo_step = make_u64(bigword_weight_shift_minus1, FRAC_BPW_HI);
  combo_counter = comboFracBits(word_index) +
                  make_u64(word_index * bigword_weight_shift_minus1, 0xFFFFFFFF);
  const u32 log2_NWORDS = (WIDTH == 256 ? 8 : WIDTH == 512 ? 9 : WIDTH == 1024 ? 10 : 12) +
                          (MIDDLE == 1 ? 0 : MIDDLE == 2 ? 1 : MIDDLE == 4 ? 2 : MIDDLE == 8 ? 3 : 4) +
                          (SMALL_HEIGHT == 256 ? 8 : SMALL_HEIGHT == 512 ? 9 : SMALL_HEIGHT == 1024 ? 10 : 12) + 1;
  weight_shift = (weight_shift + log2_NWORDS + 1) % 31;

  #if RIESEL_PAIR
  u32 exponent0;
  Z61 invPair = rieselStartingWeight(THREAD_WEIGHTS, gx * G_W + me,
                                     line, true, &exponent0);
  invPair = mul(invPair, (Z61)RIESEL_INV_SCALE);
  #else
  u32 exponent0, exponent1;
  Z31 inv0 = rieselCarryStartingInverse(THREAD_WEIGHTS, 0, gx * G_W + me,
                                         line, &exponent0);
  Z31 inv1 = rieselCarryStartingInverse(THREAD_WEIGHTS, 1, gx * G_W + me,
                                         line, &exponent1);
  inv0 = riesel0Mul(inv0, (Z31)RIESEL0_INV_SCALE);
  inv1 = riesel1Mul(inv1, (Z31)RIESEL1_INV_SCALE);
  #endif

  for (i32 i = 0; i < CARRY_LEN; ++i) {
    u32 const p = G_W * gx + WIDTH * (line + i) + me;
    u32 const shift0 = weight_shift;
    combo_counter += combo_step;
    weight_shift = adjust_m31_weight_shift(weight_shift);
    u32 const shift1 = weight_shift;

    #if RIESEL_PAIR
    u32 oddExponent0 = exponent0;
    Z61 const oddInvPair = advanceRieselInverse(
      invPair, RIESEL_DELTA_ONE, (Z61)RIESEL_INV_ONE, &oddExponent0);
    Z31 const inv0 = rq0(invPair), inv1 = rq1(invPair);
    Z31 const oddInv0 = rq0(oddInvPair), oddInv1 = rq1(oddInvPair);
    GF61 const packedValue = inPair[p];
    GF31 const value0 = U2(rq0(packedValue.x), rq0(packedValue.y));
    GF31 const value1 = U2(rq1(packedValue.x), rq1(packedValue.y));
    #else
    u32 oddExponent0 = exponent0;
    u32 oddExponent1 = exponent1;
    Z31 const oddInv0 = advanceRieselCarryInverse(
      inv0, 0, RIESEL_DELTA_ONE, (Z31)RIESEL0_INV_ONE, &oddExponent0);
    Z31 const oddInv1 = advanceRieselCarryInverse(
      inv1, 1, RIESEL_DELTA_ONE, (Z31)RIESEL1_INV_ONE, &oddExponent1);
    GF31 const value0 = in0[p], value1 = in1[p];
    #endif
    bool const biglit0 = frac_bits <= FRAC_BPW_HI;
    bool const biglit1 = frac_bits >= -FRAC_BPW_HI;
    out[p] = weightAndCarryPair(SWAP_XY(in31[p]), SWAP_XY(value0), SWAP_XY(value1),
                                shift0, shift1, inv0, oddInv0, inv1, oddInv1,
                                LL != 0 || i != 0, carry,
                                biglit0, biglit1, &carry, &roundMax, &carryMax);

    combo_counter += combo_step;
    weight_shift = adjust_m31_weight_shift(weight_shift);
    #if RIESEL_PAIR
    invPair = advanceRieselInverse(oddInvPair, RIESEL_DELTA_ONE,
                                   (Z61)RIESEL_INV_ONE, &oddExponent0);
    exponent0 = oddExponent0;
    #else
    inv0 = advanceRieselCarryInverse(oddInv0, 0, RIESEL_DELTA_ONE,
                                     (Z31)RIESEL0_INV_ONE, &oddExponent0);
    inv1 = advanceRieselCarryInverse(oddInv1, 1, RIESEL_DELTA_ONE,
                                     (Z31)RIESEL1_INV_ONE, &oddExponent1);
    exponent0 = oddExponent0;
    exponent1 = oddExponent1;
    #endif
  }
#undef frac_bits
#undef weight_shift
#undef combo_counter
  carryOut[G_W * g + me] = carry;

#if ROE
  local u32 lds[G_W];
  updateStats(lds, G_W, H, bufROE, posROE, 0.0f);
#elif (STATS & (1 << (2 + MUL3)))
  local u32 lds[G_W];
  updateStats(lds, G_W, H, bufROE, posROE, carryMax);
#endif
}

/******************************************************************************/
/*  Similar to above, but for a hybrid FFT based on FP32*GF(M31^2)*GF(M61^2)  */
/******************************************************************************/

#elif FFT_TYPE == FFT323161

KERNEL(G_W) carry(P(Word2) out, CP(T2) in, u32 posROE, P(CarryABM) carryOut, BigTabFP32 THREAD_WEIGHTS, P(uint) bufROE
#if PARITY_SQUARE
                    , CP(u32) parityIn
#endif
                    ) {
  u32 g  = get_group_id(0);
  u32 me = get_local_id(0);
  u32 gx = g % NW;
  u32 gy = g / NW;
  u32 H = BIG_HEIGHT;
  u32 line = gy * CARRY_LEN;

  CP(F2) inF2 = (CP(F2)) in;
  CP(GF31) in31 = (CP(GF31)) (in + DISTGF31);
  CP(GF61) in61 = (CP(GF61)) (in + DISTGF61);

  // & vs. && to workaround spurious warning
  CarryABM carry = (LL & (me == 0) & (g == 0)) ? -2 : 0;
  float roundMax = 0;
  float carryMax = 0;

  u32 word_index = (gx * G_W * H + me * H + line) * 2;

  F base = fancyMul(THREAD_WEIGHTS[me].x, iweightStep(gx));
  u32 me_frac_bits = fracBits(me * H * 2);
  u32 step_frac_bits = weightStepFracBits(gx);
  u32 base_frac_bits = me_frac_bits + step_frac_bits;
  base = optionalDouble(base, base_frac_bits > step_frac_bits);

  // Weight is 2^[ceil(qj / n) - qj/n] where j is the word index, q is the Mersenne exponent, and n is the number of words.
#if GOOD_THOMAS3
  const u32 m31_log2_root_two = 21;
#else
  const u32 m31_log2_root_two = (u32) (((1ULL << 30) / NWORDS) % 31);
#endif
  const u32 m31_bigword_weight_shift = (NWORDS - EXP % NWORDS) * m31_log2_root_two % 31;
  const u32 m31_bigword_weight_shift_minus1 = (m31_bigword_weight_shift + 30) % 31;
#if GOOD_THOMAS3
  const u32 m61_log2_root_two = 45;
#else
  const u32 m61_log2_root_two = (u32) (((1ULL << 60) / NWORDS) % 61);
#endif
  const u32 m61_bigword_weight_shift = (NWORDS - EXP % NWORDS) * m61_log2_root_two % 61;
  const u32 m61_bigword_weight_shift_minus1 = (m61_bigword_weight_shift + 60) % 61;

  // Derive the big vs. little flags from the fractional number of bits in each word.
  // Create a 64-bit counter that tracks both weight shifts and frac_bits (adding 0xFFFFFFFF to effect the ceil operation required for weight shift).
  union { uint2 a; u64 b; } m31_combo, m61_combo;
#define frac_bits           m31_combo.a[0]
#define m31_weight_shift    m31_combo.a[1]
#define m31_combo_counter   m31_combo.b
#define m61_weight_shift    m61_combo.a[1]
#define m61_combo_counter   m61_combo.b

  const u64 m31_combo_step = ((u64) m31_bigword_weight_shift_minus1 << 32) + FRAC_BPW_HI;
  m31_combo_counter = comboFracBits(word_index) + make_u64(word_index * m31_bigword_weight_shift_minus1, 0xFFFFFFFF);
  const u64 m61_combo_step = ((u64) m61_bigword_weight_shift_minus1 << 32) + FRAC_BPW_HI;
  m61_combo_counter = comboFracBits(word_index) + make_u64(word_index * m61_bigword_weight_shift_minus1, 0xFFFFFFFF);

  // We also adjust shift amount for the fact that NTT returns results multiplied by 2*NWORDS.
#if GOOD_THOMAS3
  // Each Good--Thomas exact channel transforms 2^20 scalar words.  The
  // common +1 below accounts for the packed DGT's additional factor of two.
  const u32 log2_NWORDS = 20;
#else
  const u32 log2_NWORDS = (WIDTH == 256 ? 8 : WIDTH == 512 ? 9 : WIDTH == 1024 ? 10 : 12) +
                          (MIDDLE == 1 ? 0 : MIDDLE == 2 ? 1 : MIDDLE == 4 ? 2 : MIDDLE == 8 ? 3 : 4) +
                          (SMALL_HEIGHT == 256 ? 8 : SMALL_HEIGHT == 512 ? 9 : SMALL_HEIGHT == 1024 ? 10 : 12) + 1;
#endif
  m31_weight_shift = (m31_weight_shift + log2_NWORDS + 1) % 31;
  m61_weight_shift = (m61_weight_shift + log2_NWORDS + 1) % 61;

  for (i32 i = 0; i < CARRY_LEN; ++i) {
    u32 p = G_W * gx + WIDTH * (CARRY_LEN * gy + i) + me;

    // Generate the FP32 and second GF31 and GF61 weight shift
    F w1 = optionalDouble(fancyMul(base, THREAD_WEIGHTS[G_W + line + i].x), frac_bits > base_frac_bits);
    F w2 = optionalDouble(fancyMul(w1, IWEIGHT_STEP), frac_bits + FRAC_BPW_HI > FRAC_BPW_HI);
    u32 m31_weight_shift0 = m31_weight_shift;
    m31_combo_counter += m31_combo_step;
    m31_weight_shift = adjust_m31_weight_shift(m31_weight_shift);
    u32 m31_weight_shift1 = m31_weight_shift;
    u32 m61_weight_shift0 = m61_weight_shift;
    m61_combo_counter += m61_combo_step;
    m61_weight_shift = adjust_m61_weight_shift(m61_weight_shift);
    u32 m61_weight_shift1 = m61_weight_shift;

    // Generate big-word/little-word flags
    bool biglit0 = frac_bits <= FRAC_BPW_HI;
    bool biglit1 = frac_bits >= -FRAC_BPW_HI;   // Same as frac_bits + FRAC_BPW_HI <= FRAC_BPW_HI;

    // Compute result.  fftW leaves FP32 in natural order but the exact fields
    // in three Good--Thomas channel planes.  Gather the two scalar residues
    // of a natural pair from their respective planes.  The inverse transform
    // convention swaps packed components, hence raw odd/even below is
    // followed by the existing SWAP_XY at reconstruction.
#if GOOD_THOMAS3
    const u32 gt_line = CARRY_LEN * gy + i;
    const u32 gt_m = gt_line / SMALL_HEIGHT;
    const u32 gt_y = gt_line % SMALL_HEIGHT;
    const u32 gt_x = G_W * gx + me;
    const u32 gt_x_channel = (3 * gt_x + gt_m) & (WIDTH - 1);
    const u32 gt_channel_even = (gt_m + 2 * gt_y) % 3;
    const u32 gt_channel_odd = (gt_channel_even + 1) % 3;
    const u32 gt_even_storage =
      (gt_channel_even * SMALL_HEIGHT + gt_y) * WIDTH + gt_x_channel;
    const u32 gt_odd_storage =
      (gt_channel_odd * SMALL_HEIGHT + gt_y) * WIDTH + gt_x_channel;
    const GF31 gt_u31 = U2(in31[gt_odd_storage].x,
                           in31[gt_even_storage].y);
    const GF61 gt_u61 = U2(in61[gt_odd_storage].x,
                           in61[gt_even_storage].y);
#else
    const GF31 gt_u31 = in31[p];
    const GF61 gt_u61 = in61[p];
#endif
    out[p] = weightAndCarryPair(SWAP_XY(inF2[p]), SWAP_XY(gt_u31), SWAP_XY(gt_u61), w1, w2, m31_weight_shift0, m31_weight_shift1, m61_weight_shift0, m61_weight_shift1,
#if PARITY_SQUARE
#if PARITY_LAZY
#if PARITY_PHYSICAL
                                parityIn, gx * G_W * H + me * H + line + i,
#else
                                parityIn, p,
#endif
#else
#if PARITY_PREPARED
                                squareCoefficientParity(parityIn, p), 0,
#else
                                squareCoefficientParity(parityIn, gx * G_W * H + me * H + line + i), 0,
#endif
#endif
#endif
                                LL != 0 || i != 0, carry, biglit0, biglit1, &carry, &roundMax, &carryMax);

    // Generate weight shifts and frac_bits for next pair
    m31_combo_counter += m31_combo_step;
    m31_weight_shift = adjust_m31_weight_shift(m31_weight_shift);
    m61_combo_counter += m61_combo_step;
    m61_weight_shift = adjust_m61_weight_shift(m61_weight_shift);
  }
  carryOut[G_W * g + me] = carry;

#if ROE
  local u32 lds[G_W];
  updateStats(lds, G_W, H, bufROE, posROE, roundMax);
#elif (STATS & (1 << (2 + MUL3)))
  local u32 lds[G_W];
  updateStats(lds, G_W, H, bufROE, posROE, carryMax);
#endif
}


#else
error - missing Carry kernel implementation
#endif
