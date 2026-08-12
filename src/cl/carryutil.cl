// Copyright (C) Mihai Preda

#if CARRY64
typedef i64 CFcarry;
#else
typedef i32 CFcarry;
#endif

// The carry for the non-fused CarryA, CarryB, CarryM kernels.
// Simply use largest possible carry always as the split kernels are slow anyway (and seldomly used normally).
#if FFT_TYPE != FFT32 && FFT_TYPE != FFT31
typedef i64 CarryABM;
#else
typedef i32 CarryABM;
#endif

/********************************/
/*       Helper routines        */
/********************************/

// Return unsigned low bits (number of bits must be between 1 and 31)
#if defined(__has_builtin) && __has_builtin(__builtin_amdgcn_ubfe)
u32 OVERLOAD ulowBits(i32 u, u32 bits) { return __builtin_amdgcn_ubfe(u, 0, bits); }
#elif HAS_PTX >= 700        // szext instruction requires sm_70 support or higher
u32 OVERLOAD ulowBits(i32 u, u32 bits) { u32 res; __asm("szext.clamp.u32 %0, %1, %2;" : "=r"(res) : "r"(u), "r"(bits)); return res; }
#else
u32 OVERLOAD ulowBits(i32 u, u32 bits) { return (((u32) u << (32 - bits)) >> (32 - bits)); }
#endif
u32 OVERLOAD ulowBits(u32 u, u32 bits) { return ulowBits((i32) u, bits); }
// Return unsigned low bits (number of bits must be between 1 and 63)
u64 OVERLOAD ulowBits(i64 u, u32 bits) { return (((u64) u << (64 - bits)) >> (64 - bits)); }
u64 OVERLOAD ulowBits(u64 u, u32 bits) { return ulowBits((i64) u, bits); }

// Return unsigned low bits where number of bits is known at compile time (number of bits can be 0 to 32)
u32 OVERLOAD ulowFixedBits(i32 u, const u32 bits) { if (bits == 32) return u; return u & ((1 << bits) - 1); }
u32 OVERLOAD ulowFixedBits(u32 u, const u32 bits) { return ulowFixedBits((i32) u, bits); }
// Return unsigned low bits where number of bits is known at compile time (number of bits can be 0 to 64)
u64 OVERLOAD ulowFixedBits(i64 u, const u32 bits) { return u & ((1LL << bits) - 1); }
u64 OVERLOAD ulowFixedBits(u64 u, const u32 bits) { return ulowFixedBits((i64) u, bits); }

// Return signed low bits (number of bits must be between 1 and 31)
#if defined(__has_builtin) && __has_builtin(__builtin_amdgcn_sbfe)
i32 OVERLOAD lowBits(i32 u, u32 bits) { return __builtin_amdgcn_sbfe(u, 0, bits); }
#elif HAS_PTX >= 700        // szext instruction requires sm_70 support or higher
i32 OVERLOAD lowBits(i32 u, u32 bits) { i32 res; __asm("szext.clamp.s32 %0, %1, %2;" : "=r"(res) : "r"(u), "r"(bits)); return res; }
#else
i32 OVERLOAD lowBits(i32 u, u32 bits) { return ((u << (32 - bits)) >> (32 - bits)); }
#endif
i32 OVERLOAD lowBits(u32 u, u32 bits) { return lowBits((i32)u, bits); }
// Return signed low bits (number of bits must be between 1 and 63)
i64 OVERLOAD lowBits(i64 u, u32 bits) { return ((u << (64 - bits)) >> (64 - bits)); }
i64 OVERLOAD lowBits(u64 u, u32 bits) { return lowBits((i64)u, bits); }

// Return signed low bits (number of bits must be between 1 and 32)
#if HAS_PTX                 // szext does not return result we are looking for if bits = 32
i32 OVERLOAD lowBitsSafe32(i32 u, u32 bits) { return lowBits(u, bits); }
#else
i32 OVERLOAD lowBitsSafe32(i32 u, u32 bits) { return lowBits((u64)u, bits); }
#endif
i32 OVERLOAD lowBitsSafe32(u32 u, u32 bits) { return lowBitsSafe32((i32)u, bits); }

// Return signed low bits where number of bits is known at compile time (number of bits can be 0 to 32)
#if defined(__has_builtin) && __has_builtin(__builtin_amdgcn_sbfe)
i32 OVERLOAD lowFixedBits(i32 u, const u32 bits) { if (bits == 32) return u; return __builtin_amdgcn_sbfe(u, 0, bits); }
#elif HAS_PTX >= 700        // szext instruction requires sm_70 support or higher
i32 OVERLOAD lowFixedBits(i32 u, const u32 bits) { if (bits == 32) return u; i32 res; __asm("szext.clamp.s32 %0, %1, %2;" : "=r"(res) : "r"(u), "r"(bits)); return res; }
#else
i32 OVERLOAD lowFixedBits(i32 u, const u32 bits) { if (bits == 32) return u; return (u << (32 - bits)) >> (32 - bits); }
#endif
i32 OVERLOAD lowFixedBits(u32 u, const u32 bits) { return lowFixedBits((i32)u, bits); }
// Return signed low bits where number of bits is known at compile time (number of bits can be 1 to 63).  The two versions are the same speed on TitanV.
i64 OVERLOAD lowFixedBits(i64 u, const u32 bits) { if (bits <= 32) return lowFixedBits((i32) u, bits); return ((u << (64 - bits)) >> (64 - bits)); }
//i64 OVERLOAD lowFixedBits(i64 u, const u32 bits) { if (bits <= 32) return lowFixedBits((i32) u, bits); return (i64) ulowFixedBits(u, bits - 1) - (u & (1LL << (bits - 1))); }
i64 OVERLOAD lowFixedBits(u64 u, const u32 bits) { return lowFixedBits((i64)u, bits); }

// Extract 32 bits from a 64-bit value (starting bit offset can be 0 to 31)
#if defined(__has_builtin) && __has_builtin(__builtin_amdgcn_alignbit)
i32 xtract32(i64 x, u32 bits) { return __builtin_amdgcn_alignbit(as_int2(x).y, as_int2(x).x, bits); }
#elif HAS_PTX >= 320        // shf instruction requires sm_32 support or higher
i32 xtract32(i64 x, u32 bits) { i32 res; __asm("shf.r.clamp.b32 %0, %1, %2, %3;" : "=r"(res) : "r"(as_uint2(x).x), "r"(as_uint2(x).y), "r"(bits)); return res; }
#else
i32 xtract32(i64 x, u32 bits) { return x >> bits; }
#endif

// Extract 32 bits from a 64-bit value (starting bit offset can be 0 to 32)
#if HAS_PTX >= 320        // shf instruction requires sm_32 support or higher
i32 xtractSafe32(i64 x, u32 bits) { i32 res; __asm("shf.r.clamp.b32 %0, %1, %2, %3;" : "=r"(res) : "r"(as_uint2(x).x), "r"(as_uint2(x).y), "r"(bits)); return res; }
#else
i32 xtractSafe32(i64 x, u32 bits) { return x >> bits; }
#endif

u32 bitlen(bool b) { return EXP / NWORDS + b; }
bool test(u32 bits, u32 pos) { return (bits >> pos) & 1; }

#if FFT_FP64
// Rounding constant: 3 * 2^51, See https://stackoverflow.com/questions/17035464
#define RNDVAL (3.0 * (1ull << 51))

// Convert a double to long efficiently.  Double must be in RNDVAL+integer format.
i64 RNDVALdoubleToLong(double d) {
  int2 words = as_int2(d);
#if EXP / NWORDS >= 19
  // We extend the range to 52 bits instead of 51 by taking the sign from the negation of bit 51
  words.y ^= 0x00080000u;
  words.y = lowBits(words.y, 20);
#else
  // Take the sign from bit 50 (i.e. use lower 51 bits).
  words.y = lowBits(words.y, 19);
#endif
  return as_long(words);
}

#elif FFT_FP32
// Rounding constant: 3 * 2^22
#define RNDVAL (3.0f * (1 << 22))

// Convert a float to int efficiently.  Float must be in RNDVAL+integer format.
i32 RNDVALfloatToInt(float d) {
  int w = as_int(d);
//#if 0
// We extend the range to 23 bits instead of 22 by taking the sign from the negation of bit 22
//  w ^= 0x00800000u;
//  w = lowBits(words.y, 23);
//#else
//  // Take the sign from bit 21 (i.e. use lower 22 bits).
  w = lowBits(w, 22);
//#endif
  return w;
}

// Round an FP32 estimate of a multiple of M61 to a signed integer.  The
// historical RNDVALfloatToInt path deliberately retains only 22 signed bits.
// For the FP32+M61 hybrid near its upper useful range the quotient itself can
// require all 24 exactly representable signed-integer bits.  Round the
// magnitude with one fused add of 2^23, which keeps the result in [2^23,2^24)
// where binary32 has unit spacing, then restore the sign.  Unlike a separate
// multiply followed by cvt.rni this preserves the exact-product rounding
// supplied by the FMA.
i32 roundM61Quotient(float work, float* signedRoundoff) {
#if FP32_WIDE_QUOTIENT
  const float invM61 = 4.3368086899420177360298112034798e-19f;
  const float roundBase = 8388608.0f;  // 2^23
  float roundedMagnitude = fma(fabs(work), invM61, roundBase);
  i32 magnitude = as_int(roundedMagnitude) - as_int(roundBase);
  i32 quotient = as_uint(work) & 0x80000000u ? -magnitude : magnitude;
  *signedRoundoff = fma(work, invM61, -(float) quotient);
  return quotient;
#else
  float rounded = fma(work, 4.3368086899420177360298112034798e-19f, RNDVAL);
  i32 quotient = RNDVALfloatToInt(rounded);
  *signedRoundoff = fma(work, 4.3368086899420177360298112034798e-19f,
                        RNDVAL - rounded);
  return quotient;
#endif
}
#endif

// map abs(carry) to floats, with 2^32 corresponding to 1.0
// So that the maximum CARRY32 abs(carry), 2^31, is mapped to 0.5 (the same as the maximum ROE)
float OVERLOAD boundCarry(i32 c) { return ldexp(fabs((float) c), -32); }
float OVERLOAD boundCarry(i64 c) { return ldexp(fabs((float) (i32) (c >> 8)), -24); }

#if STATS || ROE
void updateStats(local u32 *lds, u32 num_threads, u32 num_blocks, global uint *bufROE, u32 posROE, float roundMax) {
  assert(roundMax >= 0);
  u32 me = get_local_id(0);
  u32 u32RoundMax = as_uint(roundMax);

  // Reduce to a handful of roundMax values
  // We could use shfl_down_sync (and AMD's equivalent) instead of LDS memory once num_threads < WAVEFRONT
  // (see https://github.com/mahmoudmaftah/MaxReduction-Cuda/blob/main/code/reduction_benchmarks.cu)
  while (num_threads > 8) {
    // Write roundMax for high half of threads to local memory.  Ignore threads not participating in the reduction.
    if (num_threads > WAVEFRONT) bar();
    if (me >= num_threads / 2 && me < num_threads) lds[me - num_threads / 2] = u32RoundMax;
    if (num_threads > WAVEFRONT) {
      bar();                             // work around a weird CUDA NVCC bug where two bar() calls are required???! (Titan V, CUDA 13.0, WMUL=2)
      bar();
    }
    // Low half of threads do a max
    if (me < num_threads / 2) {
      u32 highHalfMax = lds[me];
      if (u32RoundMax < highHalfMax) u32RoundMax = highHalfMax;
    }
    // Cut num threads in half, loop
    num_threads /= 2;
  }

  // The bufROE entry to update is stored in the first bufROE entry.  This value used to be passed into carryFused as an argument.
  // CUDA graphs don't allow arguments to change.  Thus, calculating posROE and storing it in bufROE works better.
  if (me < num_threads) {
    posROE = bufROE[0];
    atomic_max(bufROE + posROE + 2, u32RoundMax);

    // The second bufRoe entry is a count of the number atomic_maxes performed.  When the last atomic_max is done, increment posROE and clear the counter.
    if (me == 0) {
      u32 old_value = atomic_add(bufROE + 1, 1);
      if (old_value == num_blocks - 1) {
        bufROE[0] = posROE + 1;
        bufROE[1] = 0;
      }
    }
  }
}

#if ROE && ROE_COUNT
// Diagnostic counterpart to updateStats.  Each workgroup contributes the
// number of coefficient pairs whose reconstruction error crossed the selected
// ROE_COUNT/1000 threshold.  The output slots deliberately reuse bufROE; the
// host converts their raw uint representation before computing statistics.
void updateCountStats(local u32 *lds, u32 num_threads, u32 num_blocks, global uint *bufROE, u32 posROE, u32 count) {
  u32 me = get_local_id(0);

  while (num_threads > 1) {
    if (num_threads > WAVEFRONT) bar();
    if (me >= num_threads / 2 && me < num_threads) lds[me - num_threads / 2] = count;
    if (num_threads > WAVEFRONT) {
      bar();
      bar();
    }
    if (me < num_threads / 2) count += lds[me];
    num_threads /= 2;
  }

  if (me == 0) {
    posROE = bufROE[0];
    atomic_add(bufROE + posROE + 2, count);
    u32 old_value = atomic_add(bufROE + 1, 1);
    if (old_value == num_blocks - 1) {
      bufROE[0] = posROE + 1;
      bufROE[1] = 0;
    }
  }
}
#endif
#endif

#if 0
// Check for round off errors above a threshold (default is 0.43)
void ROUNDOFF_CHECK(double x) {
#if DEBUG
#ifndef ROUNDOFF_LIMIT
#define ROUNDOFF_LIMIT 0.43
#endif
  float error = fabs(x - rint(x));
  if (error > ROUNDOFF_LIMIT) printf("Roundoff: %g %30.2f\n", error, x);
#endif
}
#endif


/***************************************************************************/
/*  From the FFT data, construct a value to normalize and carry propagate  */
/***************************************************************************/

#if FFT_TYPE == FFT64

// Apply inverse weight, add in optional carry, calculate roundoff error, convert to integer. Handle MUL3.
i64 weightAndCarryOne(T u, T invWeight, i64 inCarry, float* maxROE, int sloppy_result_is_acceptable) {

#if !MUL3

  // Convert carry into RNDVAL + carry.
  int2 tmp = as_int2(inCarry); tmp.y += as_int2(RNDVAL).y;
  double RNDVALCarry = as_double(tmp);

  // Apply inverse weight and RNDVAL+carry
  double d = fma(u, invWeight, RNDVALCarry);

  // Optionally calculate roundoff error
  float roundoff = fabs((float) fma(u, invWeight, RNDVALCarry - d));
  *maxROE = max(*maxROE, roundoff);

  // Convert to long (for CARRY32 case we don't need to strip off the RNDVAL bits)
  if (sloppy_result_is_acceptable) return as_long(d);
  else return RNDVALdoubleToLong(d);

#else  // We cannot add in the carry until after the mul by 3

  // Apply inverse weight and RNDVAL
  double d = fma(u, invWeight, RNDVAL);

  // Optionally calculate roundoff error
  float roundoff = fabs((float) fma(u, -invWeight, d - RNDVAL));
  *maxROE = max(*maxROE, roundoff);

  // Convert to long, mul by 3, and add carry
  return RNDVALdoubleToLong(d) * 3 + inCarry;

#endif
}

/**************************************************************************/
/*            Similar to above, but for an FFT based on FP32              */
/**************************************************************************/

#elif FFT_TYPE == FFT32

// Apply inverse weight, add in optional carry, calculate roundoff error, convert to integer.  Handle MUL3.
i32 weightAndCarryOne(F u, F invWeight, i32 inCarry, float* maxROE, int sloppy_result_is_acceptable) {

#if !MUL3

  // Convert carry into RNDVAL + carry.
  float RNDVALCarry = as_float(as_int(RNDVAL) + inCarry);                       // GWBUG - just the float arithmetic?  s.b. fast

  // Apply inverse weight and RNDVAL+carry
  float d = fma(u, invWeight, RNDVALCarry);

  // Optionally calculate roundoff error
  float roundoff = fabs(fma(u, invWeight, RNDVALCarry - d));
  *maxROE = max(*maxROE, roundoff);

  // Convert to int
  return RNDVALfloatToInt(d);

#else  // We cannot add in the carry until after the mul by 3

  // Apply inverse weight and RNDVAL
  float d = fma(u, invWeight, RNDVAL);

  // Optionally calculate roundoff error
  float roundoff = fabs(fma(u, -invWeight, d - RNDVAL));
  *maxROE = max(*maxROE, roundoff);

  // Convert to int, mul by 3, and add carry
  return RNDVALfloatToInt(d) * 3 + inCarry;

#endif
}

/**************************************************************************/
/*          Similar to above, but for an NTT based on GF(M31^2)           */
/**************************************************************************/

#elif FFT_TYPE == FFT31

// Apply inverse weight, add in optional carry, calculate roundoff error, convert to integer. Handle MUL3.
i64 weightAndCarryOne(Z61 u, u32 invWeight, i32 inCarry, u32* maxROE) {

  // Apply inverse weight
  u = shr(u, invWeight);

  // Convert input to balanced representation
  i32 value = get_balanced_Z31(u);

  // Optionally calculate roundoff error as proximity to M31/2.
  u32 roundoff = (u32) abs(value);
  *maxROE = max(*maxROE, roundoff);

  // Mul by 3 and add carry
#if MUL3
  return (i64)value * 3 + inCarry;
#endif
  return value + inCarry;
}

/**************************************************************************/
/*          Similar to above, but for an NTT based on GF(M61^2)           */
/**************************************************************************/

#elif FFT_TYPE == FFT61

// Apply inverse weight, add in optional carry, calculate roundoff error, convert to integer. Handle MUL3.
i64 weightAndCarryOne(Z61 u, u32 invWeight, i64 inCarry, u32* maxROE) {

  // Apply inverse weight
  u = shr(u, invWeight);

  // Convert input to balanced representation
  i64 value = get_balanced_Z61(u);

  // Optionally calculate roundoff error as proximity to M61/2.  28 bits of accuracy should be sufficient.
  u32 roundoff = (u32) abs((i32) hi32(value));
  *maxROE = max(*maxROE, roundoff);

  // Mul by 3 and add carry
#if MUL3
  value *= 3;
#endif
  return value + inCarry;
}

/**************************************************************************/
/*    Similar to above, but for a hybrid FFT based on FP64 & GF(M31^2)    */
/**************************************************************************/

#elif FFT_TYPE == FFT6431

// Apply inverse weight, add in optional carry, calculate roundoff error, convert to integer. Handle MUL3.
i96 weightAndCarryOne(T u, Z31 u31, T invWeight, u32 m31_invWeight, bool hasInCarry, i64 inCarry, float* maxROE) {

  // Apply inverse weight and get the Z31 data
  u31 = shr(u31, m31_invWeight);
  u32 n31 = get_Z31(u31);

  // The final result must be n31 mod M31.  Use FP64 data to calculate this value.
  u = fma(u, invWeight, - (double) n31);                               // This should be close to a multiple of M31
  double uInt = fma(u, 4.656612875245796924105750827168e-10, RNDVAL);  // Divide by M31 and round to int
  i64 n64 = RNDVALdoubleToLong(uInt);

  // Optionally calculate roundoff error
  float roundoff = (float) fabs(fma(u, 4.656612875245796924105750827168e-10, RNDVAL - uInt));
  *maxROE = max(*maxROE, roundoff);

  // Compute the value using i96 math
  i64 vhi = n64 >> 33;
  u64 vlo = ((u64)n64 << 31) | n31;
  i96 value = make_i96(vhi, vlo);                   // (n64 << 31) + n31
  value = sub(value, n64);                          // n64 * M31 + n31

  // Mul by 3 and add carry
#if MUL3
  value = add(value, add(value, value));
#endif
  if (hasInCarry) value = add(value, inCarry);
  return value;
}

/**************************************************************************/
/*    Similar to above, but for a hybrid FFT based on FP32 & GF(M31^2)    */
/**************************************************************************/

#elif FFT_TYPE == FFT3231

// Apply inverse weight, add in optional carry, calculate roundoff error, convert to integer. Handle MUL3.
i64 weightAndCarryOne(float uF2, Z31 u31, float F2_invWeight, u32 m31_invWeight, i32 inCarry, float* maxROE) {

  // Apply inverse weight and get the Z31 data
  u31 = shr(u31, m31_invWeight);
  u32 n31 = get_Z31(u31);

  // The final result must be n31 mod M31.  Use FP32 data to calculate this value.
  uF2 = fma(uF2, F2_invWeight, - (float) n31);                           // This should be close to a multiple of M31
  float uF2int = fma(uF2, 4.656612875245796924105750827168e-10f, RNDVAL);   // Divide by M31 and round to int
  i32 nF2 = RNDVALfloatToInt(uF2int);

  i64 v = (((i64) nF2 << 31) | n31) - nF2;         // nF2 * M31 + n31

  // Optionally calculate roundoff error
  float roundoff = fabs(fma(uF2, 4.656612875245796924105750827168e-10f, RNDVAL - uF2int));
  *maxROE = max(*maxROE, roundoff);

  // Mul by 3 and add carry
#if MUL3
  v = v * 3;
#endif
  return v + inCarry;
}

/**************************************************************************/
/*    Similar to above, but for a hybrid FFT based on FP32 & GF(M61^2)    */
/**************************************************************************/

#elif FFT_TYPE == FFT3261

#if PARITY_SQUARE
// For a weighted cyclic square, all off-diagonal products occur twice and
// vanish modulo two.  For output word 2*k, exactly one of the two possible
// diagonal inputs has an odd Crandall--Fagin multiplier.  Output word 2*k+1
// has no diagonal and is therefore even.  NWORDS is a power of two and EXP is
// odd, so the fractional comparison needs only the low log2(NWORDS) product
// bits rather than a division.
u32 squareCoefficientParity(CP(u32) parity, u32 outputPair) {
#if PARITY_PACKED
  return (parity[outputPair >> 5] >> (outputPair & 31u)) & 1u;
#elif PARITY_PHYSICAL
  u32 sourceWord = outputPair;
  u32 frac = ((u32) EXP * sourceWord) & (NWORDS - 1);
  if (frac != 0 && frac <= NWORDS / 2) sourceWord += NWORDS / 2;
  u32 const sourcePair = sourceWord >> 1;
  u32 const sourceX = sourcePair / BIG_HEIGHT;
  u32 const sourceLine = sourcePair - sourceX * BIG_HEIGHT;
  u32 const packed = parity[sourceLine * WIDTH + sourceX];
  return (packed >> (sourceWord & 1)) & 1u;
#elif PARITY_PREPARED
  return parity[outputPair] & 1u;
#else
  u32 sourceWord = outputPair;
  u32 frac = ((u32) EXP * sourceWord) & (NWORDS - 1);
  if (frac != 0 && frac <= NWORDS / 2) sourceWord += NWORDS / 2;
  u32 packed = parity[sourceWord >> 1];
  return (packed >> (sourceWord & 1)) & 1u;
#endif
}
#endif

// Apply inverse weight, add in optional carry, calculate roundoff error, convert to integer. Handle MUL3.
i96 weightAndCarryOne(float uF2, Z61 u61, float F2_invWeight, u32 m61_invWeight,
#if PARITY_SQUARE
#if PARITY_LAZY
                      CP(u32) parity, u32 parityIndex, bool parityCanBeOdd,
#else
                      u32 expectedParity,
#endif
#endif
                      bool hasInCarry, i64 inCarry, float* maxROE) {

  // Apply inverse weight and get the Z61 data
  u61 = shr(u61, m61_invWeight);
#if GOOD_THOMAS9
  // The inverse base-field DFT is deliberately left unnormalized in the
  // middle kernel.  Remove its factor nine once per scalar at the exact
  // reconstruction boundary.
  u61 = mul(u61, (Z61)2049638230412172401ul);
#endif
  u64 n61 = get_Z61(u61);

  // The final result mod M61 must be n61.  Use FP32 data to calculate how many multiples of M61 need to be added to n61.
  float n61f = (float)hi32(n61) * -4294967296.0f;                             // Estimate -n61 as a float.
  uF2 = fma(uF2, F2_invWeight, n61f);                                         // This should be close to an integer multiple of M61
  float signedRoundoff;
  i32 nF2 = roundM61Quotient(uF2, &signedRoundoff);

  // Optionally calculate roundoff error
  float roundoff = fabs(signedRoundoff);
  *maxROE = max(*maxROE, roundoff);

#if PARITY_SQUARE
  // Adjacent quotient candidates differ by the odd modulus M61.  Parity tells
  // whether the rounded candidate is wrong; the signed residual tells which
  // adjacent integer lies on the other side of the rounding boundary.
  if (roundoff >= 0.49f) {
#if PARITY_LAZY
    u32 const expectedParity = parityCanBeOdd ?
      squareCoefficientParity(parity, parityIndex) : 0u;
#endif
    if ((((u32) nF2 ^ (u32) n61) & 1u) != expectedParity) {
      nF2 += signedRoundoff < 0 ? -1 : 1;
#if ROE_COUNT
      *maxROE = max(*maxROE, 0.75f);
#endif
    }
  }
#endif

  // Compute the value using i96 math
  i32 vhi = nF2 >> 3;
  u64 vlo = ((u64)nF2 << 61) | n61;
  i96 value = make_i96(vhi, vlo);                // (nF2 << 61) + n61
  value = sub(value, nF2);                       // nF2 * M61 + n61

  // Mul by 3 and add carry
#if MUL3
  value = add(value, add(value, value));
#endif
  if (hasInCarry) value = add(value, inCarry);
  return value;
}

/**************************************************************************/
/*    Similar to above, but for an NTT based on GF(M31^2)*GF(M61^2)       */
/**************************************************************************/

#elif FFT_TYPE == FFT3161

#if GOLD_PAIR

i96 goldTimesM31Multiplier(Z31 multiplier) {
  // Gold*k = k*2^64 - k*2^32 + k.
  i96 value = make_i96((i32)multiplier, (u64)0);
  value = sub(value, make_i96((i32)0, (u64)multiplier << 32));
  return add(value, make_i96((i32)0, (u64)multiplier));
}

// Apply inverse weights and reconstruct the balanced M31*Gold coefficient.
i96 weightAndCarryOne(Z31 u31, Z61 uGold, u32 m31_invWeight,
                      Z61 gold_invWeight, bool hasInCarry, i64 inCarry,
                      u32* maxROE) {
  u31 = shr(u31, m31_invWeight);
  uGold = mul(uGold, gold_invWeight);

  Z31 const n31 = get_Z31(u31);
  Z61 const nGold = get_Z61(uGold);
  Z31 const delta = sub(n31, modM31(nGold));
  Z31 const multiplier = mul(delta, (Z31)0x55555555u);  // inverse of 3 mod M31

  i96 value = add(goldTimesM31Multiplier(multiplier),
                  make_i96((i32)0, (u64)nGold));
  bool const aboveHalf = i96_hi32(value) > 0x3fffffffu ||
    (i96_hi32(value) == 0x3fffffffu &&
     i96_lo64(value) > 0x40000000bfffffffull);
  if (aboveHalf) {
    value = sub(value, make_i96((i32)0x7ffffffe,
                                0x800000017fffffffull));
  }

  // Preserve the incumbent integer ROE diagnostic scale.  This is a range
  // indicator for an exact NTT, not a floating roundoff estimate.
  u32 const balancedMultiplier = multiplier > M31 / 2 ?
                                 M31 - multiplier : multiplier;
  *maxROE = max(*maxROE, balancedMultiplier);

#if MUL3
  value = add(value, add(value, value));
#endif
  if (hasInCarry) value = add(value, inCarry);
  return value;
}

#else

// Apply inverse weight, add in optional carry, calculate roundoff error, convert to integer. Handle MUL3.
i96 weightAndCarryOne(Z31 u31, Z61 u61, u32 m31_invWeight, u32 m61_invWeight, bool hasInCarry, i64 inCarry, u32* maxROE) {

  // Apply inverse weights
  u31 = shr(u31, m31_invWeight);
  u61 = shr(u61, m61_invWeight);

  // Use chinese remainder theorem to create a 92-bit result.  Loosely copied from Yves Gallot's mersenne2 program.
  u32 n31 = get_Z31(u31);
  u61 += make_u64(hi32(M61), lo32(M61) - n31);   // u61 - u31
  u61 += shl(u61, 31);                           // u61 + (u61 << 31)

  // The resulting value will be get_Z61(u61) * M31 + n31 and if larger than ~M31*M61/2 return a negative value by subtracting M31 * M61.
  // We can save a little work by determining if the result will be large using just u61 and returning (get_Z61(u61) - M61) * M31 + n31.
  // This simplifies to get_balanced_Z61(u61) * M31 + n31.
  i64 n61 = get_balanced_Z61(modM61(u61));

  // Optionally calculate roundoff error as proximity to M61/2.  28 bits of accuracy should be sufficient.
  u32 roundoff = (u32) abs((i32)hi32(n61));
  *maxROE = max(*maxROE, roundoff);

  // Compute the value using i96 math
  i64 vhi = n61 >> 1;
  u32 vlo = ((u32)n61 << 31) | n31;
  i96 value = make_i96(vhi, vlo);                // (n61 << 31) + n31
  value = sub(value, n61);                       // n61 * M31 + n31

  // Mul by 3 and add carry
#if MUL3
  value = add(value, add(value, value));
#endif
  if (hasInCarry) value = add(value, inCarry);
  return value;
}

#endif

/**************************************************************************/
/* Exact CRT for three independent 31-bit residue planes (93 bits total) */
/**************************************************************************/

#elif FFT_TYPE == FFT31R2

i96 weightAndCarryOne(Z31 u31, Z31 u0, Z31 u1, u32 m31_invWeight,
                      Z31 riesel0_invWeight, Z31 riesel1_invWeight,
                      bool hasInCarry, i64 inCarry, u32* maxROE) {
  u31 = shr(u31, m31_invWeight);

  u32 const n31 = get_Z31(u31);
  u32 const n0 = riesel0Decode(riesel0Mul(u0, riesel0_invWeight));
  u32 const n1 = riesel1Decode(riesel1Mul(u1, riesel1_invWeight));

  // Garner reconstruction.  q0 < q1 and the constant is q0^-1 (mod q1)
  // in Montgomery form, so a Montgomery multiply of two normal operands
  // deliberately returns a normal residue here.
  u32 const delta = rieselSubExplicit(n1, n0, RIESEL_Q1);
#if RIESEL_LAZY
  u32 const t1 = riesel1Mul(delta, 346029397u);
#else
  u32 const t1 = riesel1Mul(delta, 713730645u);
#endif
  u64 const x01 = (u64)n0 + (u64)RIESEL_Q0 * t1;

  // q0*q1 inverse modulo M31.
#if RIESEL_LAZY
  u32 const t2 = mul(sub(n31, modM31(x01)), 1972039331u);
#else
  u32 const t2 = mul(sub(n31, modM31(x01)), 1564229429u);
#endif

  // q0*q1 * t2 as an unsigned 96-bit integer, then add x01.
#if RIESEL_LAZY
  u64 const q01 = 1071100245244379137ULL;
#else
  u64 const q01 = 4476934267141619713ULL;
#endif
  u64 const p0 = (u64)lo32(q01) * t2;
  u64 const p1 = (u64)hi32(q01) * t2;
  u64 const productLo = p0 + (p1 << 32);
  u32 const productHi = (u32)(p1 >> 32) + (productLo < p0);
  i96 value = add(make_i96((i32)productHi, productLo), make_i96((i64)x01));

  // Select the balanced representative modulo P=M31*q0*q1.
#if RIESEL_LAZY
  u32 const halfHi = 62346239u;
  u64 const halfLo = 15688667536053764095ULL;
#else
  u32 const halfHi = 260591871u;
  u64 const halfLo = 11664144917195653119ULL;
#endif
  bool const aboveHalf = i96_hi32(value) > halfHi ||
                         (i96_hi32(value) == halfHi && i96_lo64(value) > halfLo);
#if RIESEL_LAZY
  if (aboveHalf) value = sub(value, make_i96((i32)124692479, 12930590998397976575ULL));
#else
  if (aboveHalf) value = sub(value, make_i96((i32)521183743, 4881545760681754623ULL));
#endif

  // This is exact arithmetic.  The modulus product, rather than floating
  // roundoff, is the sole coefficient-range limit.
  *maxROE = max(*maxROE, 0u);

#if MUL3
  value = add(value, add(value, value));
#endif
  if (hasInCarry) value = add(value, inCarry);
  return value;
}

/******************************************************************************/
/*  Similar to above, but for a hybrid FFT based on FP32*GF(M31^2)*GF(M61^2)  */
/******************************************************************************/

#elif FFT_TYPE == FFT323161

#if PARITY_SQUARE
u32 squareCoefficientParity(CP(u32) parity, u32 outputPair) {
#if PARITY_PACKED
  return (parity[outputPair >> 5] >> (outputPair & 31u)) & 1u;
#elif PARITY_PHYSICAL
  u32 sourceWord = outputPair;
  u32 frac = ((u32) EXP * sourceWord) & (NWORDS - 1);
  if (frac != 0 && frac <= NWORDS / 2) sourceWord += NWORDS / 2;
  u32 const sourcePair = sourceWord >> 1;
  u32 const sourceX = sourcePair / BIG_HEIGHT;
  u32 const sourceLine = sourcePair - sourceX * BIG_HEIGHT;
  u32 const packed = parity[sourceLine * WIDTH + sourceX];
  return (packed >> (sourceWord & 1)) & 1u;
#elif PARITY_PREPARED
  return parity[outputPair] & 1u;
#else
  u32 sourceWord = outputPair;
  u32 frac = ((u32) EXP * sourceWord) & (NWORDS - 1);
  if (frac != 0 && frac <= NWORDS / 2) sourceWord += NWORDS / 2;
  u32 packed = parity[sourceWord >> 1];
  return (packed >> (sourceWord & 1)) & 1u;
#endif
}
#endif

// Apply inverse weight, add in optional carry, calculate roundoff error, convert to integer. Handle MUL3.
i128 weightAndCarryOne(float uF2, Z31 u31, Z61 u61, float F2_invWeight, u32 m31_invWeight, u32 m61_invWeight,
#if PARITY_SQUARE
                       u32 expectedParity,
#endif
                       bool hasInCarry, i64 inCarry, float* maxROE) {

  // Apply inverse weights
  u31 = shr(u31, m31_invWeight);
  u61 = shr(u61, m61_invWeight);
#if QUOTIENT_STATS || (PARITY_SQUARE && ROE_COUNT)
  u64 originalN61 = get_Z61(u61);
  float m61Estimate = (float)hi32(originalN61) * -4294967296.0f;
  float m61Work = fma(uF2, F2_invWeight, m61Estimate);
  float diagSignedRoundoff;
  i32 estimatedM61Quotient = roundM61Quotient(m61Work, &diagSignedRoundoff);
#endif
  // Use chinese remainder theorem to create a 92-bit result.  Loosely copied from Yves Gallot's mersenne2 program.
  u32 n31 = get_Z31(u31);
  u61 += make_u64(hi32(M61), lo32(M61) - n31);       // u61 - u31
  u61 += shl(u61, 31);                               // u61 + (u61 << 31)
  u64 n61 = get_Z61(modM61(u61));
  // Let's call the 92-bit CRT result n3161.  At this point, n3161 = n61 * M31 + n31.

  // The final result mod M31*M61 must be n3161.  Use FP32 data to calculate how many multiples of M31*M61 need to be added to n3161.
  float n3161f = (float)hi32(n61) * -9223372036854775808.0f;                 // Estimate -n3161 as a float.  -n61 << 31 should be close enough.
  uF2 = fma(uF2, F2_invWeight, n3161f);                                      // This should be close to an integer multiple of M31*M61
  float uF2int = fma(uF2, 2.0194839183061857038255724444152e-28f, RNDVAL);   // Divide by M31*M61 and round to int
  i32 nF2 = RNDVALfloatToInt(uF2int);

  // The final result will be nF2 * M31*M61 + n3161.  Rearranging to use as few 128-bit and 64-bit ops as possible:
  // = nF2 * M61 * M31 + n61 * M31 + n31
  // = (nF2 * M61 + n61) * M31 + n31
  // = ((nF2 << 61) - nF2 + n61) * M31 + n31
  // = (((nF2 << 61) - nF2 + n61) << 31) - ((nF2 << 61) - nF2 + n61) + n31
  // = (nF2 << 92) + ((n61 - nF2) << 31) - (nF2 << 61) - (n61 - nF2) + n31
  // = (nF2 << 92) - (nF2 << 61) + ((n61 - nF2) << 31) - (n61 - nF2) + n31
  // = (((nF2 << 31) - nF2) << 61) + ((n61 - nF2) << 31) - (n61 - nF2) + n31
  // = (((nF2 << 32) - nF2*2) << 60) + ((n61 - nF2) << 31) - (n61 - nF2) + n31

  // Compute x = (n61 - nF2)
  i64 x = (i64)n61 - nF2;
  // Compute y = ((n61 - nF2) << 31) + n31
  i128 y = make_i128(x >> 33, (x << 31) | n31);
  // Compute z = ((nF2 << 32) - nF2*2) << 60
  i64 tmp = make_i64(nF2, 0) - (i64)(nF2 + nF2);
  i128 z = make_i128(tmp >> 4, tmp << 60);

  // Put the parts together
  i128 v = sub(add(z, y), x);

  // Optionally calculate roundoff error
  float roundoff = fabs(fma(uF2, 2.0194839183061857038255724444152e-28f, RNDVAL - uF2int));
#if COEFF_RANGE_STATS
  i64 coeffHi = i128_hi64(v);
  u64 coeffLo = i128_lo64(v);
#if COEFF_RANGE_STATS == 23
  // Balanced range of (2^23-1)*(2^61-1):
  // floor(P/2) = 0x7ffffefffffffffc00000.
  bool coeffTooPositive = coeffHi > 0x7ffffLL ||
                          (coeffHi == 0x7ffffLL && coeffLo > 0xefffffffffc00000ULL);
  bool coeffTooNegative = coeffHi < -0x80000LL ||
                          (coeffHi == -0x80000LL && coeffLo < 0x1000000000400000ULL);
  if (coeffTooPositive || coeffTooNegative) roundoff = max(roundoff, 4.0f);
#elif COEFF_RANGE_STATS == 24
  // Exact balanced range of qC*(2^61-1), qC=7*2^21-1=14680063.
  // floor(P/2) = 0xdffffefffffffff900000.
  bool coeffTooPositive = coeffHi > 0xdffffLL ||
                          (coeffHi == 0xdffffLL && coeffLo > 0xefffffffff900000ULL);
  bool coeffTooNegative = coeffHi < -0xe0000LL ||
                          (coeffHi == -0xe0000LL && coeffLo < 0x1000000000700000ULL);
  if (coeffTooPositive || coeffTooNegative) roundoff = max(roundoff, 4.0f);
#elif COEFF_RANGE_STATS >= 65 && COEFF_RANGE_STATS <= 126
  i64 coeffLimitHi = (i64)1 << (COEFF_RANGE_STATS - 64);
  bool coeffOutsidePower = coeffHi >= coeffLimitHi || coeffHi < -coeffLimitHi ||
                           (coeffHi == -coeffLimitHi && coeffLo == 0);
  if (coeffOutsidePower) roundoff = max(roundoff, 4.0f);
#endif
#endif
#if PARITY_SQUARE && ROE_COUNT
  u32 exactParity = ((u32) nF2 ^ (u32) n61 ^ n31) & 1u;
  u32 estimatedParity = ((u32) estimatedM61Quotient ^ (u32) originalN61) & 1u;
  // Recover the FP32+M61 quotient error with the otherwise redundant M31
  // channel.  ROE_COUNT=700 counts odd errors; 701 counts cases where the
  // FP residual would choose the wrong adjacent quotient.
  Z31 diagOriginalM61Mod31 = modM31(originalN61);
  Z31 diagEstimatedMod31 = add(modM31((i64) estimatedM61Quotient * ((1LL << 30) - 1)), diagOriginalM61Mod31);
  Z31 diagErrorMod31 = mul(sub(diagEstimatedMod31, n31), M31 - 2);
  i32 diagError = diagErrorMod31 > M31 / 2 ? (i32) (diagErrorMod31 - M31) : (i32) diagErrorMod31;
  bool oddMismatch = estimatedParity != expectedParity;
  i32 proposedDelta = diagSignedRoundoff < 0 ? -1 : 1;
  bool wrongDirection = oddMismatch && proposedDelta != -diagError;
#if ROE_COUNT == 701
  if (wrongDirection) roundoff = max(roundoff, 4.0f);
#elif ROE_COUNT == 702
  // Parity-free folded-decoder gate: raw quotient errors must be ternary.
  if (abs(diagError) > 1) roundoff = max(roundoff, 4.0f);
#elif ROE_COUNT == 703
  // Every nonzero raw quotient error must lie in the candidate window.
  if (diagError != 0 && fabs(diagSignedRoundoff) < 0.49f)
    roundoff = max(roundoff, 4.0f);
#elif ROE_COUNT == 704
  // Count the parity-free decoder's complete residual candidate population.
  if (fabs(diagSignedRoundoff) >= 0.49f) roundoff = max(roundoff, 4.0f);
#elif ROE_COUNT == 705
  if (abs(diagError) > 2) roundoff = max(roundoff, 4.0f);
#elif ROE_COUNT == 706
  if (abs(diagError) > 3) roundoff = max(roundoff, 4.0f);
#elif ROE_COUNT >= 710 && ROE_COUNT < 800
  if (wrongDirection && fabs(diagSignedRoundoff) >= (float) (ROE_COUNT - 700) * .001f) roundoff = max(roundoff, 4.0f);
#elif ROE_COUNT >= 810 && ROE_COUNT < 900
  if (oddMismatch && fabs(diagSignedRoundoff) < (float) (ROE_COUNT - 800) * .001f) roundoff = max(roundoff, 4.0f);
#else
  // expectedParity was separately validated against exactParity.  Count how
  // often the cheaper FP32+M61 reconstruction needs an odd correction.
  if (estimatedParity != expectedParity) roundoff = max(roundoff, 4.0f);
#endif
#endif
#if QUOTIENT_STATS
  // Determine q_est-q_exact without dividing a 96-bit coefficient.  Mod M31,
  // M61 is 2^30-1 and its inverse is -2.  The exact M31 residue therefore
  // recovers the small quotient error directly.
  Z31 originalM61Mod31 = modM31(originalN61);
  Z31 estimatedMod31 = add(modM31((i64) estimatedM61Quotient * ((1LL << 30) - 1)), originalM61Mod31);
  Z31 quotientErrorMod31 = mul(sub(estimatedMod31, n31), M31 - 2);
  i32 quotientError = quotientErrorMod31 > M31 / 2 ? (i32) (quotientErrorMod31 - M31) : (i32) quotientErrorMod31;
  roundoff = max(roundoff, (float) abs(quotientError));
#endif
  *maxROE = max(*maxROE, roundoff);

  // Mul by 3 and add carry
#if MUL3
  v = add(v, add(v, v));
#endif
  if (hasInCarry) v = add(v, inCarry);
  return v;
}

#else
error - missing weightAndCarryOne implementation
#endif


/************************************************************************/
/*   Split a value + carryIn into a big-or-little word and a carryOut   */
/************************************************************************/

Word OVERLOAD carryStep(i128 x, i64 *outCarry, bool isBigWord) {
  u32 nBits = bitlen(isBigWord);
  i64 w = lowBits(i128_lo64(x), nBits);
  *outCarry = i128_shrlo64(x, nBits) + (w < 0);
  return w;
}

Word OVERLOAD carryStep(i96 x, i64 *outCarry, bool isBigWord) {
  u32 nBits = bitlen(isBigWord);
  u32 nBitsLess32 = bitlen(isBigWord) - 32;

// This code can be tricky because we must not shift i32 or u32 variables by 32.
#if EXP / NWORDS >= 33
  i32 whi = lowBits(i96_mid32(x), nBitsLess32);
  *outCarry = ((i64)i96_hi64(x) - (i64)whi) >> nBitsLess32;
  return as_ulong((uint2)(i96_lo32(x), (u32)whi));
#elif EXP / NWORDS == 32
  i32 whi = xtract32(i96_lo64(x), nBitsLess32) >> 31;
  *outCarry = ((i64)i96_hi64(x) - (i64)whi) >> nBitsLess32;
  return as_ulong((uint2)(i96_lo32(x), (u32)whi));
#elif EXP / NWORDS == 31
  i32 w = lowBitsSafe32(i96_lo32(x), nBits);
  *outCarry = as_long((int2)(xtractSafe32(i96_lo64(x), nBits), xtractSafe32(i96_hi64(x), nBits))) + (w < 0);
  return w;
//  i64 w = lowBits(i96_lo64(x), nBits);
//  *outCarry = ((i96_hi64(x) << (32 - nBits)) | ((i96_lo32(x) >> 16) >> (nBits - 16))) + (w < 0);
//  return w;
#else
  i32 w = lowBits(i96_lo32(x), nBits);
  *outCarry = as_long((int2)(xtract32(i96_lo64(x), nBits), xtract32(i96_hi64(x), nBits))) + (w < 0);
  return w;
#endif
}

Word OVERLOAD carryStep(i64 x, i64 *outCarry, bool isBigWord) {
  u32 nBits = bitlen(isBigWord);
#if EXP / NWORDS >= 33
  i32 xhi = hi32(x);
  i32 whi = lowBits(xhi, nBits - 32);
  *outCarry = (xhi - whi) >> (nBits - 32);
  return (Word) as_long((int2)(lo32(x), whi));
#elif EXP / NWORDS == 32
  i32 xhi = hi32(x);
  i64 w = lowBits(x, nBits);
  xhi -= (i32)hi32(w);
  *outCarry = xhi >> (nBits - 32);
  return w;
#elif EXP / NWORDS == 31
  i32 w = lowBitsSafe32(lo32(x), nBits);
  *outCarry = (x - w) >> nBits;
  return w;
#else
  Word w = lowBits(lo32(x), nBits);
  *outCarry = (x - w) >> nBits;
  return w;
#endif
}

Word OVERLOAD carryStep(i64 x, i32 *outCarry, bool isBigWord) {
  u32 nBits = bitlen(isBigWord);
#if EXP / NWORDS >= 33
  i32 xhi = hi32(x);
  i32 w = lowBits(xhi, nBits - 32);
  *outCarry = (xhi >> (nBits - 32)) + (w < 0);
  return as_long((int2)(lo32(x), w));
#elif EXP / NWORDS == 32
  i32 xhi = hi32(x);
  i64 w = lowBits(x, nBits);
  *outCarry = (xhi >> (nBits - 32)) + (w < 0);
  return w;
#elif EXP / NWORDS == 31
  i32 w = lowBitsSafe32(lo32(x), nBits);
  *outCarry = xtractSafe32(x, nBits) + (w < 0);
  return w;
#else
  i32 w = lowBits(x, nBits);
  *outCarry = xtract32(x, nBits) + (w < 0);
  return w;
#endif
}

Word OVERLOAD carryStep(i32 x, i32 *outCarry, bool isBigWord) {
  u32 nBits = bitlen(isBigWord);
  Word w = lowBits(x, nBits);
  *outCarry = (x - w) >> nBits;
  return w;
}

/*****************************************************************/
/*  Same as CarryStep but returns a faster unsigned result.      */
/*  Used on first word of pair in carryFused.                    */
/* CarryFinal will later turn this into a balanced signed value. */
/*****************************************************************/

Word OVERLOAD carryStepUnsignedSloppy(i128 x, i64 *outCarry, bool isBigWord) {
  const u32 bigwordBits = EXP / NWORDS + 1;
  u32 nBits = bitlen(isBigWord);

// Return a Word using the big word size.  Big word size is a constant which allows for more optimization.
  u64 w = ulowFixedBits(i128_lo64(x), bigwordBits);
  x = i128_masklo64(x, ~((u64)1 << (bigwordBits - 1)));
  *outCarry = i128_shrlo64(x, nBits);
  return w;
}

Word OVERLOAD carryStepUnsignedSloppy(i96 x, i64 *outCarry, bool isBigWord) {
  const u32 bigwordBits = EXP / NWORDS + 1;
  u32 nBits = bitlen(isBigWord);

// Return a Word using the big word size.  Big word size is a constant which allows for more optimization.
#if EXP / NWORDS >= 32                                  // nBits is 32 or more
  i64 xhi = as_ulong((uint2)(i96_mid32(x) & ~((1 << (bigwordBits - 32)) - 1), i96_hi32(x)));
  *outCarry = xhi >> (nBits - 32);
  return as_ulong((uint2)(i96_lo32(x), ulowFixedBits(i96_mid32(x), bigwordBits - 32)));
#elif EXP / NWORDS == 31 || EXP / NWORDS >= 22          // nBits = 31 or 32, fastest version. Should also work on smaller nBits.
  *outCarry = i96_hi64(x) << (32 - nBits);
  return i96_lo32(x);                                   // ulowBits(x, bigwordBits = 32);
#else                                                   // nBits less than 32
  u32 w = ulowFixedBits(i96_lo32(x), bigwordBits);
  *outCarry = as_long((int2)(xtract32(as_long((int2)(i96_lo32(x) - w, i96_mid32(x))), nBits), xtract32(i96_hi64(x), nBits)));
  return w;
#endif
}

Word OVERLOAD carryStepUnsignedSloppy(i64 x, i64 *outCarry, bool isBigWord) {
  u32 nBits = bitlen(isBigWord);
  *outCarry = x >> nBits;
  return ulowBits(x, nBits);
}

Word OVERLOAD carryStepUnsignedSloppy(i64 x, i32 *outCarry, bool isBigWord) {
  u32 nBits = bitlen(isBigWord);
  *outCarry = xtract32(x, nBits);
  return ulowBits(x, nBits);
}

Word OVERLOAD carryStepUnsignedSloppy(i32 x, i32 *outCarry, bool isBigWord) {
  u32 nBits = bitlen(isBigWord);
  *outCarry = x >> nBits;
  return ulowBits(x, nBits);
}

/**********************************************************************/
/*  Same as CarryStep but may return a faster big word signed result. */
/*  Used on second word of pair in carryFused when not near max BPW.  */
/*  Also used on first word in carryFinal when not near max BPW.      */
/**********************************************************************/

// We only allow sloppy results when not near the maximum bits-per-word.  For now, this is defined as 1.1 bits below maxbpw.
// No studies have been done on reducing this 1,1 value since this is a rather minor optimization.  Since the preprocessor can't
// handle floats, the MAXBPW value passed in is 100 * maxbpw.
#define SLOPPY_MAXBPW   (MAXBPW - 110)
#define ACTUAL_BPW      (EXP / (NWORDS / 100))

Word OVERLOAD carryStepSignedSloppy(i128 x, i64 *outCarry, bool isBigWord) {
#if ACTUAL_BPW > SLOPPY_MAXBPW
  return carryStep(x, outCarry, isBigWord);
#else

//GW:  Need to compare to simple carryStep
  
// Return a Word using the big word size.  Big word size is a constant which allows for more optimization.
  const u32 bigwordBits = EXP / NWORDS + 1;
  u32 nBits = bitlen(isBigWord);
  u64 xlo = i128_lo64(x);
  u64 xlo_topbit = xlo & ((u64)1 << (bigwordBits - 1));
  i64 w = ulowFixedBits(xlo, bigwordBits - 1) - xlo_topbit;
  *outCarry = i128_shrlo64(add(x, xlo_topbit), nBits);
  return w;
#endif
}

Word OVERLOAD carryStepSignedSloppy(i96 x, i64 *outCarry, bool isBigWord) {
#if ACTUAL_BPW > SLOPPY_MAXBPW
  return carryStep(x, outCarry, isBigWord);
#else

// Return a Word using the big word size.  Big word size is a constant which allows for more optimization.
  const u32 bigwordBits = EXP / NWORDS + 1;
  u32 nBits = bitlen(isBigWord);
#if EXP / NWORDS >= 32                                  // nBits is 32 or more
  return carryStep(x, outCarry, isBigWord);             // Should be just as fast as code below
//  u32 xmid_topbit = i96_mid32(x) & (1 << (bigwordBits - 32 - 1));
//  i32 whi = ulowFixedBits(i96_mid32(x), bigwordBits - 32 - 1) - xmid_topbit;
//  i64 xhi = i96_hi64(x) + xmid_topbit;
//  *outCarry = xhi >> (nBits - 32);
//  return as_long((int2)(i96_lo32(x), whi));
#elif EXP / NWORDS == 31 || (SLOPPY_MAXBPW >= 3200 && EXP / NWORDS >= 22) // nBits = 31 or 32, bigwordBits = 32 (or allowed to create 32-bit word for better performance)
  i32 w = i96_lo32(x);                                  // lowBits(x, bigwordBits = 32);
  *outCarry = (i96_hi64(x) + (w < 0)) << (32 - nBits);
  return w;
#else                                                   // nBits less than 32
  return carryStep(x, outCarry, isBigWord);             // Should be faster than code below
//  i32 w = lowFixedBits(i96_lo32(x), bigwordBits);
//  *outCarry = (as_long((int2)(xtract32(i96_lo64(x), bigwordBits), xtract32(i96_hi64(x), bigwordBits))) + (w < 0)) << (bigwordBits - nBits);
//  return w;
#endif
#endif
}

Word OVERLOAD carryStepSignedSloppy(i64 x, i64 *outCarry, bool isBigWord) {
#if ACTUAL_BPW > SLOPPY_MAXBPW
  return carryStep(x, outCarry, isBigWord);
#else

  // We're unlikely to find code that is better than carryStep
  return carryStep(x, outCarry, isBigWord);
#endif
}

Word OVERLOAD carryStepSignedSloppy(i64 x, i32 *outCarry, bool isBigWord) {
#if ACTUAL_BPW > SLOPPY_MAXBPW
  return carryStep(x, outCarry, isBigWord);
#else

//GW: I need to look at PTX code generated by the code below vs. carryStep

// Return a Word using the big word size.  Big word size is a constant which allows for more optimization.
  const u32 bigwordBits = EXP / NWORDS + 1;
  u32 nBits = bitlen(isBigWord);
#if EXP / NWORDS >= 32                                  // nBits is 32 or more
  u64 x_topbit = x & ((u64)1 << (bigwordBits - 1));
  i64 w = ulowFixedBits(x, bigwordBits - 1) - x_topbit;
  i32 xhi = (i32)hi32(x) + (i32)hi32(x_topbit);
  *outCarry = xhi >> (nBits - 32);
  return w;
// nBits = 31 or 32, bigwordBits = 32 (or allowed to create 32-bit word for better performance).  For reasons I don't fully understand the sloppy
// case fails if BPW is too low.  Probably something to do with a small BPW with sloppy 32-bit values would require CARRY_LONG to work properly.
// Not a major concern as end users should avoid small BPW as there is probably a more efficient NTT that could be used.
#elif EXP / NWORDS == 31 || (EXP / NWORDS >= 23 && SLOPPY_MAXBPW >= 3200)        
  i32 w = x;                                            // lowBits(x, bigwordBits = 32);
  *outCarry = ((i32)hi32(x) + (w < 0)) << (32 - nBits);
  return w;
#else                                                   // nBits less than 32         //GWBUG - is there a faster version?  Is this faster than plain old carryStep? No
//  u32 x_topbit = (u32) x & (1 << (bigwordBits - 1));
//  i32 w = ulowFixedBits((u32) x, bigwordBits - 1) - x_topbit;
//  *outCarry = (i64)(x + x_topbit) >> nBits;
//  return w;
  return carryStep(x, outCarry, isBigWord);
#endif
#endif
}

Word OVERLOAD carryStepSignedSloppy(i32 x, i32 *outCarry, bool isBigWord) {
  return carryStep(x, outCarry, isBigWord);
}



// Carry propagation from word and carry.  Used by carryB.cl.
Word2 carryWord(Word2 a, CarryABM* carry, bool b1, bool b2) {
  a.x = carryStep(a.x + *carry, carry, b1);
  a.y = carryStep(a.y + *carry, carry, b2);
  return a;
}

/**************************************************************************/
/*     Do this last, it depends on weightAndCarryOne defined above        */
/**************************************************************************/

/* Support both 32-bit and 64-bit carries */

#if WordSize <= 4
#define iCARRY i32
#include "carryinc.cl"
#undef iCARRY
#endif

#if FFT_TYPE != FFT32 && FFT_TYPE != FFT31
#define iCARRY i64
#include "carryinc.cl"
#undef iCARRY
#endif
