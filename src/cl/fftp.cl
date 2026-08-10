// Copyright (C) Mihai Preda

#include "base.cl"
#include "fftwidth.cl"
#include "weight.cl"
#include "middle.cl"

#if FFT_TYPE == FFT64

// fftPremul: weight words with IBDWT weights followed by FFT-width.
KERNEL(G_W) fftP(P(T2) out, CP(Word2) in, Trig smallTrig, BigTab THREAD_WEIGHTS) {
  local T2 lds[LDS_BYTES / sizeof(T2)];
  T2 u[NW];

  u32 g = get_group_id(0);
  u32 me = get_local_id(0);

  in += g * WIDTH;

  T base = optionalHalve(fancyMul(THREAD_WEIGHTS[me].y, THREAD_WEIGHTS[G_W + g].y));

  for (u32 i = 0; i < NW; ++i) {
    T w1 = i == 0 ? base : optionalHalve(fancyMul(base, fweightStep(i)));
    T w2 = optionalHalve(fancyMul(w1, WEIGHT_STEP));
    u32 p = G_W * i + me;
    u[i] = U2(in[p].x * w1, in[p].y * w2);
  }

  fft_WIDTH(lds, u, smallTrig, 1, me);

  writeCarryFusedLine(u, out, g, me);
}


/**************************************************************************/
/*            Similar to above, but for an FFT based on FP32              */
/**************************************************************************/

#elif FFT_TYPE == FFT32

// fftPremul: weight words with IBDWT weights followed by FFT-width.
KERNEL(G_W) fftP(P(F2) out, CP(Word2) in, TrigFP32 smallTrig, BigTabFP32 THREAD_WEIGHTS) {
  local F2 lds[LDS_BYTES / sizeof(F2)];
  F2 u[NW];

  u32 g = get_group_id(0);
  u32 me = get_local_id(0);

  in += g * WIDTH;

  u32 word_index = (me * BIG_HEIGHT + g) * 2;

  u32 me_frac_bits = fracBits(me * BIG_HEIGHT * 2);
  u32 line_frac_bits = fracBits(g * 2);
  u32 base_frac_bits = me_frac_bits + line_frac_bits;
  F base = optionalHalve(fancyMul(THREAD_WEIGHTS[me].y, THREAD_WEIGHTS[G_W + g].y), base_frac_bits > line_frac_bits);

  const u32 frac_bits_bigstep = fracBits(G_W * BIG_HEIGHT * 2);

  u32 frac_bits = base_frac_bits;
  for (u32 i = 0; i < NW; ++i) {
    F w1 = i == 0 ? base : optionalHalve(fancyMul(base, fweightStep(i)), frac_bits > base_frac_bits);
    F w2 = optionalHalve(fancyMul(w1, WEIGHT_STEP), frac_bits + FRAC_BPW_HI > FRAC_BPW_HI);
    u32 p = G_W * i + me;
    u[i] = U2(in[p].x * w1, in[p].y * w2);
    // Generate frac_bits for next pair
    frac_bits += frac_bits_bigstep;
  }

  fft_WIDTH(lds, u, smallTrig, 1, me);

  writeCarryFusedLine(u, out, g, me);
}


/**************************************************************************/
/*     Two independently reduced Riesel fields vectorized in one u64      */
/**************************************************************************/

#elif FFT_TYPE == FFT31R2 && RIESEL_PAIR && NTT_GF61

KERNEL(G_W) fftP(P(T2) out, CP(Word2) in, Trig smallTrig,
                  BigTab THREAD_WEIGHTS) {
  local GF61 lds[LDS_BYTES / sizeof(GF61)];
  GF61 u[NW];

  u32 const g = get_group_id(0);
  u32 const me = get_local_id(0);
  P(GF61) fieldOut = (P(GF61))(out + DISTGF61);
  TrigGF61 fieldTrig = (TrigGF61)(smallTrig + DISTWTRIGGF61);
  in += g * WIDTH;

  u32 exponent;
  Z61 weight = rieselStartingWeight(THREAD_WEIGHTS, me, g, false,
                                    &exponent);
  for (u32 i = 0; i < NW; ++i) {
    u32 const p = G_W * i + me;
    u32 oddExponent = exponent;
    Z61 const oddWeight = advanceRieselForward(
      weight, RIESEL_DELTA_ONE, (Z61)RIESEL_FWD_ONE, &oddExponent);
    u[i] = U2(mul(make_Z61_word(in[p].x), weight),
              mul(make_Z61_word(in[p].y), oddWeight));
    weight = advanceRieselForward(weight, RIESEL_DELTA_X,
                                  (Z61)RIESEL_FWD_X, &exponent);
  }

  fft_WIDTH(lds, u, fieldTrig, 1, me);
  writeCarryFusedLine(u, fieldOut, g, me);
}


/**************************************************************************/
/*          Similar to above, but for an NTT based on GF(M31^2)           */
/**************************************************************************/

#elif FFT_TYPE == FFT31 || (FFT_TYPE == FFT31R2 && !RIESEL_FIELD && (!RIESEL_FOUR || NTT_GF31))

// fftPremul: weight words with IBDWT weights followed by FFT-width.
KERNEL(G_W) fftP(P(GF31) out, CP(Word2) in, TrigGF31 smallTrig) {
  local GF31 lds[LDS_BYTES / sizeof(GF31)];
  GF31 u[NW];

  u32 g = get_group_id(0);
  u32 me = get_local_id(0);

  in += g * WIDTH;

  u32 word_index = (me * BIG_HEIGHT + g) * 2;

  // Weight is 2^[ceil(qj / n) - qj/n] where j is the word index, q is the Mersenne exponent, and n is the number of words.
  const u32 log2_root_two = (u32) (((1ULL << 30) / NWORDS) % 31);
  const u32 bigword_weight_shift = (NWORDS - EXP % NWORDS) * log2_root_two % 31;
  const u32 bigword_weight_shift_minus1 = (bigword_weight_shift + 30) % 31;

  // Derive the big vs. little flags from the fractional number of bits in each word.
  // Create a 64-bit counter that tracks both weight shifts and frac_bits (adding 0xFFFFFFFF to effect the ceil operation required for weight shift).
  union { uint2 a; u64 b; } combo;
#define frac_bits           combo.a[0]
#define weight_shift        combo.a[1]
#define combo_counter       combo.b

  const u64 combo_step = make_u64(bigword_weight_shift_minus1, FRAC_BPW_HI);
  const u64 combo_bigstep = (comboFracBits(G_W * BIG_HEIGHT * 2 - 1) + make_u64((G_W * BIG_HEIGHT * 2 - 1) * bigword_weight_shift_minus1, 0)) % (31ULL << 32);
  combo_counter = comboFracBits(word_index) + make_u64(word_index * bigword_weight_shift_minus1, 0xFFFFFFFF);
  weight_shift = weight_shift % 31;

  for (u32 i = 0; i < NW; ++i) {
    u32 p = G_W * i + me;
    // Generate the second weight shift
    u32 weight_shift0 = weight_shift;
    combo_counter += combo_step;
    if (weight_shift > 31) weight_shift -= 31;
    u32 weight_shift1 = weight_shift;
    // Convert and weight inputs
    u[i] = U2(shl(make_Z31(in[p].x), weight_shift0), shl(make_Z31(in[p].y), weight_shift1));      // Form a GF31 from each pair of input words
    // Generate weight shifts and frac_bits for next pair
    combo_counter += combo_bigstep;
    if (weight_shift > 31) weight_shift -= 31;
  }

  fft_WIDTH(lds, u, smallTrig, 1, me);

  writeCarryFusedLine(u, out, g, me);
}


/**************************************************************************/
/*       One independent 32-bit quadratic Riesel-prime transform         */
/**************************************************************************/

#elif FFT_TYPE == FFT31R2 && RIESEL_FIELD

KERNEL(G_W) fftP(P(T2) out, CP(Word2) in, Trig smallTrig,
                  BigTab THREAD_WEIGHTS) {
  local GF31 lds[LDS_BYTES / sizeof(GF31)];
  GF31 u[NW];

  u32 const g = get_group_id(0);
  u32 const me = get_local_id(0);
  P(GF31) fieldOut = (P(GF31))(out + DISTGF31);
  TrigGF31 fieldTrig = (TrigGF31)(smallTrig + DISTWTRIGGF31);
  in += g * WIDTH;

  u32 exponent;
  Z31 weight = rieselFieldStartingWeight(THREAD_WEIGHTS, me, g, false,
                                          &exponent);
  for (u32 i = 0; i < NW; ++i) {
    u32 const p = G_W * i + me;
    u32 oddExponent = exponent;
    Z31 const oddWeight = advanceRieselFieldForward(
      weight, RIESEL_DELTA_ONE, (Z31)RF_FWD_ONE, &oddExponent);
    u[i] = U2(mul(make_Z31_word(in[p].x), weight),
              mul(make_Z31_word(in[p].y), oddWeight));
    weight = advanceRieselFieldForward(weight, RIESEL_DELTA_X,
                                       (Z31)RF_FWD_X, &exponent);
  }

  fft_WIDTH(lds, u, fieldTrig, 1, me);
  writeCarryFusedLine(u, fieldOut, g, me);
}


/**************************************************************************/
/*          Similar to above, but for an NTT based on GF(M61^2)           */
/**************************************************************************/

#elif FFT_TYPE == FFT61 || (FFT_TYPE == FFT31R2 && RIESEL_FOUR && NTT_GF61 && !NTT_GF31)

// fftPremul: weight words with IBDWT weights followed by FFT-width.
KERNEL(G_W) fftP(P(GF61) out, CP(Word2) in, TrigGF61 smallTrig) {
  local GF61 lds[LDS_BYTES / sizeof(GF61)];
  GF61 u[NW];

  u32 g = get_group_id(0);
  u32 me = get_local_id(0);

#if RIESEL_FOUR
  out += DISTGF61;
  smallTrig += DISTWTRIGGF61;
#endif

  in += g * WIDTH;

  u32 word_index = (me * BIG_HEIGHT + g) * 2;

  // Weight is 2^[ceil(qj / n) - qj/n] where j is the word index, q is the Mersenne exponent, and n is the number of words.
  // Weights can be applied with shifts because 2 is the 60th root GF61.
  // Let s be the shift amount for word 1.  The shift amount for word x is ceil(x * (s - 1) + num_big_words_less_than_x) % 61.
  const u32 log2_root_two = (u32) (((1ULL << 60) / NWORDS) % 61);
  const u32 bigword_weight_shift = (NWORDS - EXP % NWORDS) * log2_root_two % 61;
  const u32 bigword_weight_shift_minus1 = (bigword_weight_shift + 60) % 61;

  // Derive the big vs. little flags from the fractional number of bits in each word.
  // Create a 64-bit counter that tracks both weight shifts and frac_bits (adding 0xFFFFFFFF to effect the ceil operation required for weight shift).
  union { uint2 a; u64 b; } combo;
#define frac_bits	combo.a[0]
#define weight_shift	combo.a[1]
#define combo_counter	combo.b

  const u64 combo_step = make_u64(bigword_weight_shift_minus1, FRAC_BPW_HI);
  const u64 combo_bigstep = (comboFracBits(G_W * BIG_HEIGHT * 2 - 1) + make_u64((G_W * BIG_HEIGHT * 2 - 1) * bigword_weight_shift_minus1, 0)) % (61ULL << 32);
  combo_counter = comboFracBits(word_index) + make_u64(word_index * bigword_weight_shift_minus1, 0xFFFFFFFF);
  weight_shift = weight_shift % 61;

  for (u32 i = 0; i < NW; ++i) {
    u32 p = G_W * i + me;
    // Generate the second weight shift
    u32 weight_shift0 = weight_shift;
    combo_counter += combo_step;
    if (weight_shift > 61) weight_shift -= 61;
    u32 weight_shift1 = weight_shift;
    // Convert and weight input
    u[i] = U2(shl(make_Z61(in[p].x), weight_shift0), shl(make_Z61(in[p].y), weight_shift1));      // Form a GF61 from each pair of input words
    // Generate weight shifts and frac_bits for next pair
    combo_counter += combo_bigstep;
    if (weight_shift > 61) weight_shift -= 61;
  }

  fft_WIDTH(lds, u, smallTrig, 1, me);

  writeCarryFusedLine(u, out, g, me);
}


/**************************************************************************/
/*    Similar to above, but for a hybrid FFT based on FP64 & GF(M31^2)    */
/**************************************************************************/

#elif FFT_TYPE == FFT6431

// fftPremul: weight words with IBDWT weights followed by FFT-width.
KERNEL(G_W) fftP(P(T2) out, CP(Word2) in, Trig smallTrig, BigTab THREAD_WEIGHTS) {
  local T2 lds[LDS_BYTES / sizeof(T2)];
  local GF31 *lds31 = (local GF31 *) lds;
  T2 u[NW];
  GF31 u31[NW];

  u32 g = get_group_id(0);
  u32 me = get_local_id(0);

  P(GF31) out31 = (P(GF31)) (out + DISTGF31);
  TrigGF31 smallTrig31 = (TrigGF31) (smallTrig + DISTWTRIGGF31);

  in += g * WIDTH;

  T base = optionalHalve(fancyMul(THREAD_WEIGHTS[me].y, THREAD_WEIGHTS[G_W + g].y));

  u32 word_index = (me * BIG_HEIGHT + g) * 2;

  // Weight is 2^[ceil(qj / n) - qj/n] where j is the word index, q is the Mersenne exponent, and n is the number of words.
  // Let s be the shift amount for word 1.  The shift amount for word x is ceil(x * (s - 1) + num_big_words_less_than_x) % 31.
  const u32 log2_root_two = (u32) (((1ULL << 30) / NWORDS) % 31);
  const u32 bigword_weight_shift = (NWORDS - EXP % NWORDS) * log2_root_two % 31;
  const u32 bigword_weight_shift_minus1 = (bigword_weight_shift + 30) % 31;

  // Derive the big vs. little flags from the fractional number of bits in each word.
  // Create a 64-bit counter that tracks both weight shifts and frac_bits (adding 0xFFFFFFFF to effect the ceil operation required for weight shift).
  union { uint2 a; u64 b; } combo;
#define frac_bits           combo.a[0]
#define weight_shift        combo.a[1]
#define combo_counter       combo.b

  const u64 combo_step = make_u64(bigword_weight_shift_minus1, FRAC_BPW_HI);
  const u64 combo_bigstep = (comboFracBits(G_W * BIG_HEIGHT * 2 - 1) + make_u64((G_W * BIG_HEIGHT * 2 - 1) * bigword_weight_shift_minus1, 0)) % (31ULL << 32);
  combo_counter = comboFracBits(word_index) + make_u64(word_index * bigword_weight_shift_minus1, 0xFFFFFFFF);
  weight_shift = weight_shift % 31;

  for (u32 i = 0; i < NW; ++i) {
    u32 p = G_W * i + me;
    // Generate the FP64 weights and the second GF31 weight shift
    T w1 = i == 0 ? base : optionalHalve(fancyMul(base, fweightStep(i)));
    T w2 = optionalHalve(fancyMul(w1, WEIGHT_STEP));
    u32 weight_shift0 = weight_shift;
    combo_counter += combo_step;
    if (weight_shift > 31) weight_shift -= 31;
    u32 weight_shift1 = weight_shift;
    // Convert and weight input
    u[i] = U2(in[p].x * w1, in[p].y * w2);
    u31[i] = U2(shl(make_Z31(in[p].x), weight_shift0), shl(make_Z31(in[p].y), weight_shift1));      // Form a GF31 from each pair of input words
    // Generate weight shifts and frac_bits for next pair
    combo_counter += combo_bigstep;
    if (weight_shift > 31) weight_shift -= 31;
  }

  fft_WIDTH(lds, u, smallTrig, 1, me);
  writeCarryFusedLine(u, out, g, me);

  fft_WIDTH(lds31, u31, smallTrig31, 1, me);
  writeCarryFusedLine(u31, out31, g, me);
}


/**************************************************************************/
/*    Similar to above, but for a hybrid FFT based on FP32 & GF(M31^2)    */
/**************************************************************************/

#elif FFT_TYPE == FFT3231

// fftPremul: weight words with IBDWT weights followed by FFT-width.
KERNEL(G_W) fftP(P(T2) out, CP(Word2) in, Trig smallTrig, BigTabFP32 THREAD_WEIGHTS) {
  local F2 ldsF2[LDS_BYTES / sizeof(F2)];
  local GF31 *lds31 = (local GF31 *) ldsF2;
  F2 uF2[NW];
  GF31 u31[NW];

  u32 g = get_group_id(0);
  u32 me = get_local_id(0);

  P(F2) outF2 = (P(F2)) out;
  TrigFP32 smallTrigF2 = (TrigFP32) smallTrig;
  P(GF31) out31 = (P(GF31)) (out + DISTGF31);
  TrigGF31 smallTrig31 = (TrigGF31) (smallTrig + DISTWTRIGGF31);

  in += g * WIDTH;

  u32 word_index = (me * BIG_HEIGHT + g) * 2;

  u32 me_frac_bits = fracBits(me * BIG_HEIGHT * 2);
  u32 line_frac_bits = fracBits(g * 2);
  u32 base_frac_bits = me_frac_bits + line_frac_bits;
  F base = optionalHalve(fancyMul(THREAD_WEIGHTS[me].y, THREAD_WEIGHTS[G_W + g].y), base_frac_bits > line_frac_bits);

  // Weight is 2^[ceil(qj / n) - qj/n] where j is the word index, q is the Mersenne exponent, and n is the number of words.
  // Let s be the shift amount for word 1.  The shift amount for word x is ceil(x * (s - 1) + num_big_words_less_than_x) % 31.
  const u32 log2_root_two = (u32) (((1ULL << 30) / NWORDS) % 31);
  const u32 bigword_weight_shift = (NWORDS - EXP % NWORDS) * log2_root_two % 31;
  const u32 bigword_weight_shift_minus1 = (bigword_weight_shift + 30) % 31;

  // Derive the big vs. little flags from the fractional number of bits in each word.
  // Create a 64-bit counter that tracks both weight shifts and frac_bits (adding 0xFFFFFFFF to effect the ceil operation required for weight shift).
  union { uint2 a; u64 b; } combo;
#define frac_bits           combo.a[0]
#define weight_shift        combo.a[1]
#define combo_counter       combo.b

  const u64 combo_step = make_u64(bigword_weight_shift_minus1, FRAC_BPW_HI);
  const u64 combo_bigstep = (comboFracBits(G_W * BIG_HEIGHT * 2 - 1) + make_u64((G_W * BIG_HEIGHT * 2 - 1) * bigword_weight_shift_minus1, 0)) % (31ULL << 32);
  combo_counter = comboFracBits(word_index) + make_u64(word_index * bigword_weight_shift_minus1, 0xFFFFFFFF);
  weight_shift = weight_shift % 31;

  for (u32 i = 0; i < NW; ++i) {
    u32 p = G_W * i + me;
    // Generate the FP32 weights and the second GF31 weight shift
    F w1 = i == 0 ? base : optionalHalve(fancyMul(base, fweightStep(i)), frac_bits > base_frac_bits);
    F w2 = optionalHalve(fancyMul(w1, WEIGHT_STEP), frac_bits + FRAC_BPW_HI > FRAC_BPW_HI);
    u32 weight_shift0 = weight_shift;
    combo_counter += combo_step;
    if (weight_shift > 31) weight_shift -= 31;
    u32 weight_shift1 = weight_shift;
    // Convert and weight input
    uF2[i] = U2(in[p].x * w1, in[p].y * w2);
    u31[i] = U2(shl(make_Z31(in[p].x), weight_shift0), shl(make_Z31(in[p].y), weight_shift1));      // Form a GF31 from each pair of input words

    // Generate weight shifts and frac_bits for next pair
    combo_counter += combo_bigstep;
    if (weight_shift > 31) weight_shift -= 31;
  }

  fft_WIDTH(ldsF2, uF2, smallTrigF2, 1, me);
  writeCarryFusedLine(uF2, outF2, g, me);

  fft_WIDTH(lds31, u31, smallTrig31, 1, me);
  writeCarryFusedLine(u31, out31, g, me);
}


/**************************************************************************/
/*    Similar to above, but for a hybrid FFT based on FP32 & GF(M61^2)    */
/**************************************************************************/

#elif FFT_TYPE == FFT3261

// fftPremul: weight words with IBDWT weights followed by FFT-width.
KERNEL(G_W) fftP(P(T2) out, CP(Word2) in, Trig smallTrig, BigTabFP32 THREAD_WEIGHTS) {
  local GF61 lds61[LDS_BYTES / sizeof(GF61)];
  local F2 *ldsF2 = (local F2 *) lds61;
  F2 uF2[NW];
  GF61 u61[NW];

  u32 g = get_group_id(0);
  u32 me = get_local_id(0);

  P(F2) outF2 = (P(F2)) out;
  TrigFP32 smallTrigF2 = (TrigFP32) smallTrig;
  P(GF61) out61 = (P(GF61)) (out + DISTGF61);
  TrigGF61 smallTrig61 = (TrigGF61) (smallTrig + DISTWTRIGGF61);

  CP(Word2) naturalIn = in;
  in += g * WIDTH;

  u32 word_index = (me * BIG_HEIGHT + g) * 2;

  u32 me_frac_bits = fracBits(me * BIG_HEIGHT * 2);
  u32 line_frac_bits = fracBits(g * 2);
  u32 base_frac_bits = me_frac_bits + line_frac_bits;
  F base = optionalHalve(fancyMul(THREAD_WEIGHTS[me].y, THREAD_WEIGHTS[G_W + g].y), base_frac_bits > line_frac_bits);

  // Weight is 2^[ceil(qj / n) - qj/n] where j is the word index, q is the Mersenne exponent, and n is the number of words.
  // Let s be the shift amount for word 1.  The shift amount for word x is ceil(x * (s - 1) + num_big_words_less_than_x) % 61.
#if GOOD_THOMAS9
  // NWORDS=9*2^19 and 30*NWORDS == 1 (mod 61), so 2^30 is
  // the Crandall--Fagin N-th root of two in M61.
  const u32 log2_root_two = 30;
#else
  const u32 log2_root_two = (u32) (((1ULL << 60) / NWORDS) % 61);
#endif
  const u32 bigword_weight_shift = (NWORDS - EXP % NWORDS) * log2_root_two % 61;
  const u32 bigword_weight_shift_minus1 = (bigword_weight_shift + 60) % 61;

  // Derive the big vs. little flags from the fractional number of bits in each word.
  // Create a 64-bit counter that tracks both weight shifts and frac_bits (adding 0xFFFFFFFF to effect the ceil operation required for weight shift).
  union { uint2 a; u64 b; } combo;
#define frac_bits           combo.a[0]
#define weight_shift        combo.a[1]
#define combo_counter       combo.b

  const u64 combo_step = make_u64(bigword_weight_shift_minus1, FRAC_BPW_HI);
  const u64 combo_bigstep = (comboFracBits(G_W * BIG_HEIGHT * 2 - 1) + make_u64((G_W * BIG_HEIGHT * 2 - 1) * bigword_weight_shift_minus1, 0)) % (61ULL << 32);
  combo_counter = comboFracBits(word_index) + make_u64(word_index * bigword_weight_shift_minus1, 0xFFFFFFFF);
  weight_shift = weight_shift % 61;

  for (u32 i = 0; i < NW; ++i) {
    u32 p = G_W * i + me;
    // Generate the FP32 weights and the second GF61 weight shift
    F w1 = i == 0 ? base : optionalHalve(fancyMul(base, fweightStep(i)), frac_bits > base_frac_bits);
    F w2 = optionalHalve(fancyMul(w1, WEIGHT_STEP), frac_bits + FRAC_BPW_HI > FRAC_BPW_HI);
    u32 weight_shift0 = weight_shift;
    combo_counter += combo_step;
    if (weight_shift > 61) weight_shift -= 61;
    u32 weight_shift1 = weight_shift;
    // Convert and weight input
    uF2[i] = U2(in[p].x * w1, in[p].y * w2);
#if GOOD_THOMAS9
    // The floating transform stays in natural order.  The M61 transform uses
    // the CRT map Z/(9*2^19) -> Z/9 x Z/2^19.  One GF61 value contains two
    // scalar words from the same channel; because 512 == -1 (mod 9), the
    // natural middle coordinates are y+4*c and y+4*c+5.  57 is 9^-1 mod 512.
    const u32 gt_y = g % SMALL_HEIGHT;
    const u32 gt_channel = g / SMALL_HEIGHT;
    const u32 gt_m_even = (gt_y + 4 * gt_channel) % 9;
    const u32 gt_m_odd = (gt_m_even + 5) % 9;
    const u32 gt_x_even = ((p + WIDTH - gt_m_even) * 57u) & (WIDTH - 1);
    const u32 gt_x_odd = ((p + WIDTH - gt_m_odd) * 57u) & (WIDTH - 1);
    const u32 gt_line_even = gt_m_even * SMALL_HEIGHT + gt_y;
    const u32 gt_line_odd = gt_m_odd * SMALL_HEIGHT + gt_y;
    const Word gt_even = naturalIn[gt_line_even * WIDTH + gt_x_even].x;
    const Word gt_odd = naturalIn[gt_line_odd * WIDTH + gt_x_odd].y;
    const u32 gt_word_even = 2 * (gt_x_even * BIG_HEIGHT + gt_line_even);
    const u32 gt_word_odd = 2 * (gt_x_odd * BIG_HEIGHT + gt_line_odd) + 1;
    union { uint2 a; u64 b; } gt_even_combo, gt_odd_combo;
    gt_even_combo.b = comboFracBits(gt_word_even) +
      make_u64(gt_word_even * bigword_weight_shift_minus1, 0xFFFFFFFF);
    gt_odd_combo.b = comboFracBits(gt_word_odd) +
      make_u64(gt_word_odd * bigword_weight_shift_minus1, 0xFFFFFFFF);
    u61[i] = U2(shl(make_Z61(gt_even), gt_even_combo.a[1] % 61),
                  shl(make_Z61(gt_odd), gt_odd_combo.a[1] % 61));
#else
    u61[i] = U2(shl(make_Z61(in[p].x), weight_shift0), shl(make_Z61(in[p].y), weight_shift1));      // Form a GF61 from each pair of input words
#endif
    // Generate weight shifts and frac_bits for next pair
    combo_counter += combo_bigstep;
    if (weight_shift > 61) weight_shift -= 61;
  }

  fft_WIDTH(ldsF2, uF2, smallTrigF2, 1, me);
  writeCarryFusedLine(uF2, outF2, g, me);

  fft_WIDTH(lds61, u61, smallTrig61, 1, me);
  writeCarryFusedLine(u61, out61, g, me);
}


/**************************************************************************/
/*    Similar to above, but for an NTT based on GF(M31^2)*GF(M61^2)       */
/**************************************************************************/

#elif FFT_TYPE == FFT3161

// fftPremul: weight words with IBDWT weights followed by FFT-width.
KERNEL(G_W) fftP(P(T2) out, CP(Word2) in, Trig smallTrig
#if GOLD_PAIR
                  , BigTab THREAD_WEIGHTS
#endif
                  ) {
  local GF61 lds61[LDS_BYTES / sizeof(GF61)];
  local GF31 *lds31 = (local GF31 *) lds61;
  GF31 u31[NW];
  GF61 u61[NW];

  u32 g = get_group_id(0);
  u32 me = get_local_id(0);

  P(GF31) out31 = (P(GF31)) (out + DISTGF31);
  TrigGF31 smallTrig31 = (TrigGF31) (smallTrig + DISTWTRIGGF31);
  P(GF61) out61 = (P(GF61)) (out + DISTGF61);
  TrigGF61 smallTrig61 = (TrigGF61) (smallTrig + DISTWTRIGGF61);

  in += g * WIDTH;

  u32 word_index = (me * BIG_HEIGHT + g) * 2;

  // Weight is 2^[ceil(qj / n) - qj/n] where j is the word index, q is the Mersenne exponent, and n is the number of words.
  // Weights can be applied with shifts because 2 is the 60th root GF61.
  // Let s be the shift amount for word 1.  The shift amount for word x is ceil(x * (s - 1) + num_big_words_less_than_x) % 61.
#if GOOD_THOMAS3
  // NWORDS=3*2^20.  Use NWORDS^-1 in the exponent groups modulo
  // 31 and 61; the historical integer-division shortcut is valid only when
  // NWORDS is a power of two.
  const u32 m31_log2_root_two = 21;
#else
  const u32 m31_log2_root_two = (u32) (((1ULL << 30) / NWORDS) % 31);
#endif
  const u32 m31_bigword_weight_shift = (NWORDS - EXP % NWORDS) * m31_log2_root_two % 31;
  const u32 m31_bigword_weight_shift_minus1 = (m31_bigword_weight_shift + 30) % 31;
#if !GOLD_PAIR
#if GOOD_THOMAS3
  const u32 m61_log2_root_two = 45;
#else
  const u32 m61_log2_root_two = (u32) (((1ULL << 60) / NWORDS) % 61);
#endif
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
  const u64 m31_combo_bigstep = (comboFracBits(G_W * BIG_HEIGHT * 2 - 1) + make_u64((G_W * BIG_HEIGHT * 2 - 1) * m31_bigword_weight_shift_minus1, 0)) % (31ULL << 32);
  m31_combo_counter = comboFracBits(word_index) + make_u64(word_index * m31_bigword_weight_shift_minus1, 0xFFFFFFFF);
  m31_weight_shift = m31_weight_shift % 31;
#if GOLD_PAIR
  u32 goldExponent;
  Z61 goldWeight = goldStartingWeight(THREAD_WEIGHTS, me, g, false,
                                      &goldExponent);
#else
  const u64 m61_combo_step = make_u64(m61_bigword_weight_shift_minus1, FRAC_BPW_HI);
  const u64 m61_combo_bigstep = (comboFracBits(G_W * BIG_HEIGHT * 2 - 1) + make_u64((G_W * BIG_HEIGHT * 2 - 1) * m61_bigword_weight_shift_minus1, 0)) % (61ULL << 32);
  m61_combo_counter = comboFracBits(word_index) + make_u64(word_index * m61_bigword_weight_shift_minus1, 0xFFFFFFFF);
  m61_weight_shift = m61_weight_shift % 61;
#endif

  for (u32 i = 0; i < NW; ++i) {
    u32 p = G_W * i + me;
    // Generate the second weight shifts
    u32 m31_weight_shift0 = m31_weight_shift;
    m31_combo_counter += m31_combo_step;
    m31_weight_shift = adjust_m31_weight_shift(m31_weight_shift);
    u32 m31_weight_shift1 = m31_weight_shift;
    // Convert and weight input
    u31[i] = U2(shl(make_Z31(in[p].x), m31_weight_shift0), shl(make_Z31(in[p].y), m31_weight_shift1));      // Form a GF31 from each pair of input words
#if GOLD_PAIR
    u32 oddExponent = goldExponent;
    Z61 const oddWeight = advanceGoldForward(
      goldWeight, GOLD_DELTA_ONE, (Z61)GOLD_FWD_ONE, &oddExponent);
    u61[i] = U2(mul(make_Z61_word(in[p].x), goldWeight),
                  mul(make_Z61_word(in[p].y), oddWeight));
    goldWeight = advanceGoldForward(goldWeight, GOLD_DELTA_X,
                                    (Z61)GOLD_FWD_X, &goldExponent);
#else
    u32 m61_weight_shift0 = m61_weight_shift;
    m61_combo_counter += m61_combo_step;
    m61_weight_shift = adjust_m61_weight_shift(m61_weight_shift);
    u32 m61_weight_shift1 = m61_weight_shift;
    u61[i] = U2(shl(make_Z61(in[p].x), m61_weight_shift0), shl(make_Z61(in[p].y), m61_weight_shift1));      // Form a GF61 from each pair of input words
#endif

// Generate weight shifts and frac_bits for next pair
    m31_combo_counter += m31_combo_bigstep;
    m31_weight_shift = adjust_m31_weight_shift(m31_weight_shift);
#if !GOLD_PAIR
    m61_combo_counter += m61_combo_bigstep;
    m61_weight_shift = adjust_m61_weight_shift(m61_weight_shift);
#endif
  }

  fft_WIDTH(lds31, u31, smallTrig31, 1, me);
  writeCarryFusedLine(u31, out31, g, me);

  fft_WIDTH(lds61, u61, smallTrig61, 1, me);
  writeCarryFusedLine(u61, out61, g, me);
}


/******************************************************************************/
/*  Similar to above, but for a hybrid FFT based on FP32*GF(M31^2)*GF(M61^2)  */
/******************************************************************************/

#elif FFT_TYPE == FFT323161

// fftPremul: weight words with IBDWT weights followed by FFT-width.
KERNEL(G_W) fftP(P(T2) out, CP(Word2) in, Trig smallTrig, BigTabFP32 THREAD_WEIGHTS) {
  local GF61 lds61[LDS_BYTES / sizeof(GF61)];
  local F2 *ldsF2 = (local F2 *) lds61;
  local GF31 *lds31 = (local GF31 *) lds61;
  F2 uF2[NW];
  GF31 u31[NW];
  GF61 u61[NW];

  u32 g = get_group_id(0);
  u32 me = get_local_id(0);

  P(F2) outF2 = (P(F2)) out;
  TrigFP32 smallTrigF2 = (TrigFP32) smallTrig;
  P(GF31) out31 = (P(GF31)) (out + DISTGF31);
  TrigGF31 smallTrig31 = (TrigGF31) (smallTrig + DISTWTRIGGF31);
  P(GF61) out61 = (P(GF61)) (out + DISTGF61);
  TrigGF61 smallTrig61 = (TrigGF61) (smallTrig + DISTWTRIGGF61);

  CP(Word2) naturalIn = in;
  in += g * WIDTH;

  u32 word_index = (me * BIG_HEIGHT + g) * 2;

  u32 me_frac_bits = fracBits(me * BIG_HEIGHT * 2);
  u32 line_frac_bits = fracBits(g * 2);
  u32 base_frac_bits = me_frac_bits + line_frac_bits;
  F base = optionalHalve(fancyMul(THREAD_WEIGHTS[me].y, THREAD_WEIGHTS[G_W + g].y), base_frac_bits > line_frac_bits);

  // Weight is 2^[ceil(qj / n) - qj/n] where j is the word index, q is the Mersenne exponent, and n is the number of words.
  // Weights can be applied with shifts because 2 is the 60th root GF61.
  // Let s be the shift amount for word 1.  The shift amount for word x is ceil(x * (s - 1) + num_big_words_less_than_x) % 61.
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

  const u64 m31_combo_step = make_u64(m31_bigword_weight_shift_minus1, FRAC_BPW_HI);
  const u64 m31_combo_bigstep = (comboFracBits(G_W * BIG_HEIGHT * 2 - 1) + make_u64((G_W * BIG_HEIGHT * 2 - 1) * m31_bigword_weight_shift_minus1, 0)) % (31ULL << 32);
  m31_combo_counter = comboFracBits(word_index) + make_u64(word_index * m31_bigword_weight_shift_minus1, 0xFFFFFFFF);
  m31_weight_shift = m31_weight_shift % 31;
  const u64 m61_combo_step = make_u64(m61_bigword_weight_shift_minus1, FRAC_BPW_HI);
  const u64 m61_combo_bigstep = (comboFracBits(G_W * BIG_HEIGHT * 2 - 1) + make_u64((G_W * BIG_HEIGHT * 2 - 1) * m61_bigword_weight_shift_minus1, 0)) % (61ULL << 32);
  m61_combo_counter = comboFracBits(word_index) + make_u64(word_index * m61_bigword_weight_shift_minus1, 0xFFFFFFFF);
  m61_weight_shift = m61_weight_shift % 61;

  for (u32 i = 0; i < NW; ++i) {
    u32 p = G_W * i + me;
    // Generate the FP32 weights and the second GF31 and GF61 weight shift
    F w1 = i == 0 ? base : optionalHalve(fancyMul(base, fweightStep(i)), frac_bits > base_frac_bits);
    F w2 = optionalHalve(fancyMul(w1, WEIGHT_STEP), frac_bits + FRAC_BPW_HI > FRAC_BPW_HI);
    u32 m31_weight_shift0 = m31_weight_shift;
    m31_combo_counter += m31_combo_step;
    m31_weight_shift = adjust_m31_weight_shift(m31_weight_shift);
    u32 m31_weight_shift1 = m31_weight_shift;
    u32 m61_weight_shift0 = m61_weight_shift;
    m61_combo_counter += m61_combo_step;
    m61_weight_shift = adjust_m61_weight_shift(m61_weight_shift);
    u32 m61_weight_shift1 = m61_weight_shift;
    // Convert and weight input.  FP32 retains the ordinary natural-order
    // 3M transform.  The exact fields use the Good--Thomas channel layout:
    // channel=n mod 3 and base=n mod 2^20.  A packed channel pair therefore
    // gathers its even and odd components from two different natural pairs.
    uF2[i] = U2(in[p].x * w1, in[p].y * w2);
#if GOOD_THOMAS3
    const u32 gt_y = g % SMALL_HEIGHT;
    const u32 gt_channel = g / SMALL_HEIGHT;
    const u32 gt_m_even = (gt_channel + gt_y) % 3;
    const u32 gt_m_odd = (gt_channel + gt_y + 2) % 3;
    // 683 is 3^-1 modulo 1024.  GOOD_THOMAS3 currently admits only the
    // production-sized WIDTH=1024 middle-3 geometry.
    const u32 gt_x_even = ((p + WIDTH - gt_m_even) * 683u) & (WIDTH - 1);
    const u32 gt_x_odd = ((p + WIDTH - gt_m_odd) * 683u) & (WIDTH - 1);
    const u32 gt_line_even = gt_m_even * SMALL_HEIGHT + gt_y;
    const u32 gt_line_odd = gt_m_odd * SMALL_HEIGHT + gt_y;
    const Word gt_even = naturalIn[gt_line_even * WIDTH + gt_x_even].x;
    const Word gt_odd = naturalIn[gt_line_odd * WIDTH + gt_x_odd].y;
    const u32 gt_word_even = 2 * (gt_x_even * BIG_HEIGHT + gt_line_even);
    const u32 gt_word_odd = 2 * (gt_x_odd * BIG_HEIGHT + gt_line_odd) + 1;
    union { uint2 a; u64 b; } gt_m31_even, gt_m31_odd;
    union { uint2 a; u64 b; } gt_m61_even, gt_m61_odd;
    gt_m31_even.b = comboFracBits(gt_word_even) +
      make_u64(gt_word_even * m31_bigword_weight_shift_minus1, 0xFFFFFFFF);
    gt_m31_odd.b = comboFracBits(gt_word_odd) +
      make_u64(gt_word_odd * m31_bigword_weight_shift_minus1, 0xFFFFFFFF);
    gt_m61_even.b = comboFracBits(gt_word_even) +
      make_u64(gt_word_even * m61_bigword_weight_shift_minus1, 0xFFFFFFFF);
    gt_m61_odd.b = comboFracBits(gt_word_odd) +
      make_u64(gt_word_odd * m61_bigword_weight_shift_minus1, 0xFFFFFFFF);
    u31[i] = U2(shl(make_Z31(gt_even), gt_m31_even.a[1] % 31),
                shl(make_Z31(gt_odd), gt_m31_odd.a[1] % 31));
    u61[i] = U2(shl(make_Z61(gt_even), gt_m61_even.a[1] % 61),
                shl(make_Z61(gt_odd), gt_m61_odd.a[1] % 61));
#else
    u31[i] = U2(shl(make_Z31(in[p].x), m31_weight_shift0), shl(make_Z31(in[p].y), m31_weight_shift1));      // Form a GF31 from each pair of input words
    u61[i] = U2(shl(make_Z61(in[p].x), m61_weight_shift0), shl(make_Z61(in[p].y), m61_weight_shift1));      // Form a GF61 from each pair of input words
#endif

// Generate weight shifts and frac_bits for next pair
    m31_combo_counter += m31_combo_bigstep;
    m31_weight_shift = adjust_m31_weight_shift(m31_weight_shift);
    m61_combo_counter += m61_combo_bigstep;
    m61_weight_shift = adjust_m61_weight_shift(m61_weight_shift);
  }

  fft_WIDTH(ldsF2, uF2, smallTrigF2, 1, me);
  writeCarryFusedLine(uF2, outF2, g, me);

  fft_WIDTH(lds31, u31, smallTrig31, 1, me);
  writeCarryFusedLine(u31, out31, g, me);

  fft_WIDTH(lds61, u61, smallTrig61, 1, me);
  writeCarryFusedLine(u61, out61, g, me);
}


#else
error - missing FFTp kernel implementation
#endif
