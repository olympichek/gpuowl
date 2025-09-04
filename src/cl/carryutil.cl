// Copyright (C) Mihai Preda

#include "base.cl"
#include "math.cl"

#if CARRY64
typedef i64 CFcarry;
#else
typedef i32 CFcarry;
#endif

// The carry for the non-fused CarryA, CarryB, CarryM kernels.
// Simply use large carry always as the split kernels are slow anyway (and seldomly used normally).
typedef i64 CarryABM;

#if defined(__has_builtin) && __has_builtin(__builtin_amdgcn_sbfe)
i32 lowBits(i32 u, u32 bits) { return __builtin_amdgcn_sbfe(u, 0, bits); }
#else
i32 lowBits(i32 u, u32 bits) { return ((u << (32 - bits)) >> (32 - bits)); }
#endif

#if defined(__has_builtin) && __has_builtin(__builtin_amdgcn_ubfe)
i32 ulowBits(i32 u, u32 bits) { return __builtin_amdgcn_ubfe(u, 0, bits); }
#else
i32 ulowBits(i32 u, u32 bits) { u32 uu = (u32) u; return ((uu << (32 - bits)) >> (32 - bits)); }
#endif

#if defined(__has_builtin) && __has_builtin(__builtin_amdgcn_alignbit)
i32 xtract32(i64 x, u32 bits) { return __builtin_amdgcn_alignbit(as_int2(x).y, as_int2(x).x, bits); }
#else
i32 xtract32(i64 x, u32 bits) { return x >> bits; }
#endif

u32 bitlen(bool b) { return EXP / NWORDS + b; }
bool test(u32 bits, u32 pos) { return (bits >> pos) & 1; }

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

Word OVERLOAD carryStep(i96 x, i64 *outCarry, bool isBigWord) {
  u32 nBits = bitlen(isBigWord);
  i64 xhi = i96_hi64(x);
#if EXP / NWORDS > 32
  i64 w = lowBits(xhi, nBits - 32);
  xhi -= w;
  *outCarry = xhi >> (nBits - 32);
  return (w << 32) | i96_lo32(x);
#else
  const u32 upshift = 33 - EXP / NWORDS;                // Shift up 33 - EXP / NWORDS so that we can work on the hi64 bits
  const u32 upshift_mask = (u32) ((1ULL << upshift) - 1);
  xhi = (xhi << upshift) + (i96_lo32(x) >> (32 - upshift));
  i64 w = lowBits(xhi, nBits - upshift);
  xhi -= w;
  *outCarry = xhi >> (nBits - upshift);
  return (w << upshift) | (i96_lo32(x) & upshift_mask);
#endif
}

Word OVERLOAD carryStep(i64 x, i64 *outCarry, bool isBigWord) {
  u32 nBits = bitlen(isBigWord);
  Word w = lowBits(x, nBits);
  x -= w;
  *outCarry = x >> nBits;
  return w;
}

Word OVERLOAD carryStep(i64 x, i32 *outCarry, bool isBigWord) {
  u32 nBits = bitlen(isBigWord);
  Word w = lowBits(x, nBits);
  *outCarry = xtract32(x, nBits) + (w < 0);
  return w;
}

Word OVERLOAD carryStep(i32 x, i32 *outCarry, bool isBigWord) {
  u32 nBits = bitlen(isBigWord);
  Word w = lowBits(x, nBits);
  *outCarry = (x - w) >> nBits;
  return w;
}

Word OVERLOAD carryStepSloppy(i96 x, i64 *outCarry, bool isBigWord) {
  u32 nBits = bitlen(isBigWord);
  i64 xhi = i96_hi64(x);
#if EXP / NWORDS >= 32                                  // nBits is 32 or more
  *outCarry = xhi >> (nBits - 32);
  u64 xlo = i96_lo64(x);
  return (xlo << (64 - nBits)) >> (64 - nBits);         // ulowBits(xlo, nBits);
#elif EXP / NWORDS == 31                                // nBits = 31 or 32
  *outCarry = isBigWord ? xhi : (xhi << 1) + ((u32) i96_lo32(x) >> 31);
  u64 xlo = i96_lo64(x);
  return (xlo << (64 - nBits)) >> (64 - nBits);         // ulowBits(xlo, nBits);
#else                                                   // nBits less than 32
  *outCarry = (xhi << (32 - nBits)) + ((u32) i96_lo32(x) >> nBits);
  return ulowBits(i96_lo32(x), nBits);
#endif
}

Word OVERLOAD carryStepSloppy(i64 x, i64 *outCarry, bool isBigWord) {
  u32 nBits = bitlen(isBigWord);
  *outCarry = x >> nBits;
  return ulowBits(x, nBits);
}

Word OVERLOAD carryStepSloppy(i64 x, i32 *outCarry, bool isBigWord) {
  u32 nBits = bitlen(isBigWord);
  *outCarry = xtract32(x, nBits);
  return ulowBits(x, nBits);
}

Word OVERLOAD carryStepSloppy(i32 x, i32 *outCarry, bool isBigWord) {
  u32 nBits = bitlen(isBigWord);
  *outCarry = x >> nBits;
  return ulowBits(x, nBits);
}

// Carry propagation from word and carry.
Word2 carryWord(Word2 a, CarryABM* carry, bool b1, bool b2) {
  a.x = carryStep(a.x + *carry, carry, b1);
  a.y = carryStep(a.y + *carry, carry, b2);
  return a;
}

// map abs(carry) to floats, with 2^32 corresponding to 1.0
// So that the maximum CARRY32 abs(carry), 2^31, is mapped to 0.5 (the same as the maximum ROE)
float OVERLOAD boundCarry(i32 c) { return ldexp(fabs((float) c), -32); }
float OVERLOAD boundCarry(i64 c) { return ldexp(fabs((float) (i32) (c >> 8)), -24); }

#if STATS || ROE
void updateStats(global uint *bufROE, u32 posROE, float roundMax) {
  assert(roundMax >= 0);
  // work_group_reduce_max() allocates an additional 256Bytes LDS for a 64lane workgroup, so avoid it.
  // u32 groupRound = work_group_reduce_max(as_uint(roundMax));
  // if (get_local_id(0) == 0) { atomic_max(bufROE + posROE, groupRound); }

  // Do the reduction directly over global mem.
  atomic_max(bufROE + posROE, as_uint(roundMax));
}
#endif


#if FFT_FP64 & !COMBO_FFT

// Rounding constant: 3 * 2^51, See https://stackoverflow.com/questions/17035464
#define RNDVAL (3.0 * (1l << 51))

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

// Apply inverse weight, add in optional carry, calculate roundoff error, convert to integer. Handle MUL3.
i64 weightAndCarryOne(T u, T invWeight, i64 inCarry, float* maxROE, int sloppy_result_is_acceptable) {

#if !MUL3

  // Convert carry into RNDVAL + carry.
  int2 tmp = as_int2(inCarry); tmp.y += as_int2(RNDVAL).y;
  double RNDVALCarry = as_double(tmp);

  // Apply inverse weight and RNDVAL+carry
  double d = fma(u, invWeight, RNDVALCarry);

  // Optionally calculate roundoff error
  float roundoff = fabs((float) fma(u, -invWeight, d - RNDVALCarry));
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
/*          Similar to above, but for an NTT based on GF(M31^2)           */
/**************************************************************************/

#elif NTT_GF31 & !COMBO_FFT

// Apply inverse weight, add in optional carry, calculate roundoff error, convert to integer. Handle MUL3.
i64 weightAndCarryOne(Z61 u, u32 invWeight, i64 inCarry, u32* maxROE) {

  // Apply inverse weight
  u = shr(u, invWeight);

  // Convert input to balanced representation
  i32 value = get_balanced_Z31(u);

  // Optionally calculate roundoff error as proximity to M31/2.
  u32 roundoff = (u32) abs(value);
  *maxROE = max(*maxROE, roundoff);

  // Mul by 3 and add carry
#if MUL3
  value *= 3;
#endif
  return value + inCarry;
}


/**************************************************************************/
/*          Similar to above, but for an NTT based on GF(M61^2)           */
/**************************************************************************/

#elif NTT_GF61 & !COMBO_FFT

// Apply inverse weight, add in optional carry, calculate roundoff error, convert to integer. Handle MUL3.
i64 weightAndCarryOne(Z61 u, u32 invWeight, i64 inCarry, u32* maxROE) {

  // Apply inverse weight
  u = shr(u, invWeight);

  // Convert input to balanced representation
  i64 value = get_balanced_Z61(u);

  // Optionally calculate roundoff error as proximity to M61/2.  28 bits of accuracy should be sufficient.
  u32 roundoff = (u32) abs((i32) (value >> 32));
  *maxROE = max(*maxROE, roundoff);

  // Mul by 3 and add carry
#if MUL3
  value *= 3;
#endif
  return value + inCarry;
}


/**************************************************************************/
/*    Similar to above, but for an NTT based on GF(M31^2)*GF(M61^2)       */
/**************************************************************************/

#elif NTT_GF31 & NTT_GF61

// Apply inverse weight, add in optional carry, calculate roundoff error, convert to integer. Handle MUL3.
i96 weightAndCarryOne(Z31 u31, Z61 u61, u32 m31_invWeight, u32 m61_invWeight, i64 inCarry, u32* maxROE) {

  // Apply inverse weights
  u31 = shr(u31, m31_invWeight);
  u61 = shr(u61, m61_invWeight);

  // Use chinese remainder theorem to create a 92-bit result.  Loosely copied from Yves Gallot's mersenne2 program.
  u32 uu31 = get_Z31(u31);
  u61 = sub(u61, make_Z61(uu31));                // u61 + u31
  u61 = add(u61, shl(u61, 31));                  // u61 + (u61 << 31)
  u64 u = get_Z61(u61);
  i96 value(u >> 1, ((u32) u << 31) + uu31);     // (u61<<31) + u31
  value.sub((u64) u);

  // Convert input to balanced representation by subtracting M61*M31
  if (value.hi32() value.sub(i96(0x0FFFFFFF, 0xDFFFFFFF, 0x80000001);

  // Optionally calculate roundoff error as proximity to M61*M31/2.  27 bits of accuracy should be sufficient.
  u32 roundoff = (u32) abs((i32) value.hi32());
  *maxROE = max(*maxROE, roundoff);

  // Mul by 3 and add carry
#if MUL3
  value.mul(3);
#endif
  return value.add(i96((u32)(inCarry >> 63), (u64) incarry));
}


#else
error - missing carryUtil implementation
#endif



/**************************************************************************/
/*     Do this last, it depends on weightAndCarryOne defined above        */
/**************************************************************************/

/* Support both 32-bit and 64-bit carries */

#define iCARRY i32
#include "carryinc.cl"
#undef iCARRY

#define iCARRY i64
#include "carryinc.cl"
#undef iCARRY

