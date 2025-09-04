// Copyright (C) Mihai Preda

#pragma once

#include "base.cl"

// Access parts of a 64-bit value

u32 lo32(u64 x) { return (u32) x; }
u32 hi32(u64 x) { return (u32) (x >> 32); }

// A primitive partial implementation of an i96/u96 integer type
typedef union {
  struct { u32 lo32; u32 mid32; u32 hi32; } a;
  struct { u32 lo32; u64 hi64; } b;
  struct { u64 lo64; u32 hi32; } c;
} i96;
i96 OVERLOAD make_i96(u32 h, u32 m, u32 l) { i96 val; val.a.hi32 = h, val.a.mid32 = m, val.a.lo32 = l; return val; }
i96 OVERLOAD make_i96(u64 h, u32 l) { i96 val; val.b.hi64 = h, val.b.lo32 = l; return val; }
i96 OVERLOAD make_i96(u32 h, u64 l) { i96 val;  val.c.hi32 = h, val.c.lo64 = l; return val; }
void i96_add(i96 *val, i96 x) { u64 lo64 = val->c.lo64 + x.c.lo64; val->c.hi32 += x.c.hi32 + (lo64 < val->c.lo64); val->c.lo64 = lo64; }
void OVERLOAD i96_sub(i96 *val, i96 x) { u64 lo64 = val->c.lo64 - x.c.lo64; val->c.hi32 -= x.c.hi32 + (lo64 > val->c.lo64); val->c.lo64 = lo64; }
void OVERLOAD i96_sub(i96 *val, u64 x) { i96_sub(val, make_i96(0, x)); }
void i96_mul(i96 *val, u32 x) { u64 t = (u64)val->b.lo32 * x; val->b.lo32 = (u32)t; val->b.hi64 = val->b.hi64 * x + (t >> 32); }
u32 i96_hi32(i96 val) { return val.a.hi32; }
u64 i96_hi64(i96 val) { return val.b.hi64; }
u32 i96_lo32(i96 val) { return val.a.lo32; }

// The X2 family of macros and SWAP are #defines because OpenCL does not allow pass by reference.
// With NTT support added, we need to turn these macros into overloaded routines.
#define X2(a, b)              X2_internal(&(a), &(b))                 // a = a + b, b = a - b
#define X2conjb(a, b)         X2conjb_internal(&(a), &(b))            // X2(a, conjugate(b))
#define X2mul2b_conjb(a, b)   X2mul2_conjb_internal(&(a), &(b))       // X2(a, 2*b), b = conjugate(b)
#define X2_mul_t4(a, b)       X2_mul_t4_internal(&(a), &(b))          // X2(a, b), b = mul_t4(b)
#define X2_mul_t8(a, b)       X2_mul_t8_internal(&(a), &(b))          // X2(a, b), b = mul_t8(b)
#define X2_mul_3t8(a, b)      X2_mul_3t8_internal(&(a), &(b))         // X2(a, b), b = mul_3t8(b)
#define X2_conja(a, b)        X2_conja_internal(&(a), &(b))           // X2(a, b), a = conjugate(a)                           // NOT USED
#define X2_conjb(a, b)        X2_conjb_internal(&(a), &(b))           // X2(a, b), b = conjugate(b)
#define SWAP(a, b)            SWAP_internal(&(a), &(b))               // a = b, b = a
#define SWAP_XY(a)            U2((a).y, (a).x)                        // Swap real and imaginary components of a

#if FFT_FP64

T2 OVERLOAD conjugate(T2 a) { return U2(a.x, -a.y); }

// Multiply by 2 without using floating point instructions.  This is a little sloppy as an input of zero returns 2^-1022.
T OVERLOAD mul2(T a) { int2 tmp = as_int2(a); tmp.y += 0x00100000; /* Bump exponent by 1 */ return (as_double(tmp)); }
T2 OVERLOAD mul2(T2 a) { return U2(mul2(a.x), mul2(a.y)); }

// Multiply by -2 without using floating point instructions.  This is a little sloppy as an input of zero returns -2^-1022.
T OVERLOAD mulminus2(T a) { int2 tmp = as_int2(a); tmp.y += 0x80100000; /* Bump exponent by 1, flip sign bit */ return (as_double(tmp)); }
T2 OVERLOAD mulminus2(T2 a) { return U2(mulminus2(a.x), mulminus2(a.y)); }

// a * (b + 1) == a * b + a
T OVERLOAD fancyMul(T a, T b)   { return fma(a, b, a); }
T2 OVERLOAD fancyMul(T2 a, T2 b) { return U2(fancyMul(a.x, b.x), fancyMul(a.y, b.y)); }

// Square a complex number
T2 OVERLOAD csq(T2 a) { return U2(fma(a.x, a.x, - a.y * a.y), mul2(a.x) * a.y); }
// a^2 + c
T2 OVERLOAD csqa(T2 a, T2 c) { return U2(fma(a.x, a.x, fma(a.y, -a.y, c.x)), fma(mul2(a.x), a.y, c.y)); }
// Same as csq(a), -a
T2 OVERLOAD csq_neg(T2 a) { return U2(fma(-a.x, a.x, a.y * a.y), mulminus2(a.x) * a.y); }                               // NOT USED

// Complex multiply
T2 OVERLOAD cmul(T2 a, T2 b) { return U2(fma(a.x, b.x, -a.y * b.y), fma(a.x, b.y, a.y * b.x)); }

T2 OVERLOAD cfma(T2 a, T2 b, T2 c) { return U2(fma(a.x, b.x, fma(a.y, -b.y, c.x)), fma(a.y, b.x, fma(a.x, b.y, c.y))); }

T2 OVERLOAD cmul_by_conjugate(T2 a, T2 b) { return cmul(a, conjugate(b)); }

// Multiply a by b and conjugate(b).  This saves 2 multiplies.
void OVERLOAD cmul_a_by_b_and_conjb(T2 *res1, T2 *res2, T2 a, T2 b) {
  T axbx = a.x * b.x;
  T aybx = a.y * b.x;
  res1->x = fma(a.y, -b.y, axbx), res1->y = fma(a.x,  b.y, aybx);
  res2->x = fma(a.y,  b.y, axbx), res2->y = fma(a.x, -b.y, aybx);
}

// Square a (cos,sin) complex number.  Fancy squaring returns a fancy value.  Defancy squares a fancy number returning a non-fancy number.
T2 csqTrig(T2 a) { T two_ay = mul2(a.y); return U2(fma(-two_ay, a.y, 1), a.x * two_ay); }
T2 csqTrigFancy(T2 a) { T two_ay = mul2(a.y); return U2(-two_ay * a.y, fma(a.x, two_ay, two_ay)); }
T2 csqTrigDefancy(T2 a) { T two_ay = mul2(a.y); return U2(fma (-two_ay, a.y, 1), fma(a.x, two_ay, two_ay)); }

// Cube a complex number w (cos,sin) given w^2 and w.  The squared input can be either fancy or not fancy.
// Fancy cCube takes a fancy w argument and returns a fancy value.  Defancy takes a fancy w argument and returns a non-fancy value.
T2 ccubeTrig(T2 sq, T2 w) { T tmp = mul2(sq.y); return U2(fma(tmp, -w.y, w.x), fma(tmp, w.x, -w.y)); }
T2 ccubeTrigFancy(T2 sq, T2 w) { T tmp = mul2(sq.y); return U2(fma(tmp, -w.y, w.x), fma(tmp, w.x, tmp - w.y)); }
T2 ccubeTrigDefancy(T2 sq, T2 w) { T tmp = mul2(sq.y); T wx = w.x + 1; return U2(fma(tmp, -w.y, wx), fma(tmp, wx, -w.y)); }

// Complex a * (b + 1)
// Useful for mul with twiddles of small angles, where the real part is stored with the -1 trick for increased precision
T2 cmulFancy(T2 a, T2 b) { return cfma(a, b, a); }

// Multiply a by fancy b and conjugate(fancy b).  This saves 2 FMAs.
void cmul_a_by_fancyb_and_conjfancyb(T2 *res1, T2 *res2, T2 a, T2 b) {
  T axbx = fma(a.x, b.x, a.x);
  T aybx = fma(a.y, b.x, a.y);
  res1->x = fma(a.y, -b.y, axbx), res1->y = fma(a.x,  b.y, aybx);
  res2->x = fma(a.y,  b.y, axbx), res2->y = fma(a.x, -b.y, aybx);
}

T2 OVERLOAD mul_t4(T2 a)  { return U2(-a.y, a.x); } // i.e. a * i

T2 OVERLOAD mul_t8(T2 a)  { // mul(a, U2(1, 1)) * (T)(M_SQRT1_2); }
  // One mul, two FMAs
  T ay = a.y * M_SQRT1_2;
  return U2(fma(a.x, M_SQRT1_2, -ay), fma(a.x, M_SQRT1_2, ay));
// Two adds, two muls
//  return U2(a.x - a.y, a.x + a.y) * M_SQRT1_2;
}

T2 OVERLOAD mul_3t8(T2 a) { // mul(a, U2(-1, 1)) * (T)(M_SQRT1_2); }
  // One mul, two FMAs
  T ay = a.y * M_SQRT1_2;
  return U2(fma(-a.x, M_SQRT1_2, -ay), fma(a.x, M_SQRT1_2, -ay));
// Two adds, two muls
//  return U2(-(a.x + a.y), a.x - a.y) * M_SQRT1_2;
}

// Return a+b and a-b
void OVERLOAD X2_internal(T2 *a, T2 *b) { T2 t = *a; *a = t + *b; *b = t - *b; }

// Same as X2(a, b), b = mul_t4(b)
void OVERLOAD X2_mul_t4_internal(T2 *a, T2 *b) { T2 t = *a; *a = *a + *b; t.x = t.x - b->x; b->x = b->y - t.y; b->y = t.x; }

// Same as X2(a, conjugate(b))
void OVERLOAD X2conjb_internal(T2 *a, T2 *b) { T2 t = *a; a->x = a->x + b->x; a->y = a->y - b->y; b->x = t.x - b->x; b->y = t.y + b->y; }

// Same as X2(a, b), a = conjugate(a)
void OVERLOAD X2_conja_internal(T2 *a, T2 *b) { T2 t = *a; a->x = a->x + b->x; a->y = - (a->y + b->y); *b = t - *b; }

// Same as X2(a, b), b = conjugate(b)
void OVERLOAD X2_conjb_internal(T2 *a, T2 *b) { T2 t = *a; *a = t + *b; b->x = t.x - b->x; b->y = b->y - t.y; }

void OVERLOAD SWAP_internal(T2 *a, T2 *b) { T2 t = *a; *a = *b; *b = t; }

T2 fmaT2(T a, T2 b, T2 c) { return fma(U2(a, a), b, c); }

// a = c + sin * d; b = c - sin * d;
#define fma_addsub(a, b, sin, c, d) { T2 t = c + sin * d; b = c - sin * d; a = t; }

T2 OVERLOAD addsub(T2 a) { return U2(a.x + a.y, a.x - a.y); }

// computes 2*(a.x*b.x+a.y*b.y) + i*2*(a.x*b.y+a.y*b.x)
// which happens to be the cyclical convolution (a.x, a.y)x(b.x, b.y) * 2
T2 foo2(T2 a, T2 b) { a = addsub(a); b = addsub(b); return addsub(U2(RE(a) * RE(b), IM(a) * IM(b))); }

// computes 2*[x^2+y^2 + i*(2*x*y)]. i.e. 2 * cyclical autoconvolution of (x, y)
T2 foo(T2 a) { return foo2(a, a); }

#endif



/**************************************************************************/
/*          Similar to above, but for an NTT based on GF(M31^2)           */
/**************************************************************************/

#if NTT_GF31

// bits in reduced mod M.
#define MBITS 31
#define M31 ((((Z31) 1) << MBITS) - 1)

// GWBUG:  Add asserts to make sure Z31 inputs are within 0 .. M31-1 range

Z31 OVERLOAD add(Z31 a, Z31 b) { Z31 t = a + b; return t - (t >= M31 ? M31 : 0); }     //GWBUG - an if stmt may be faster
GF31 OVERLOAD add(GF31 a, GF31 b) { return U2(add(a.x, b.x), add(a.y, b.y)); }

Z31 OVERLOAD sub(Z31 a, Z31 b) { Z31 t = a - b; return t + (t >= M31 ? M31 : 0); }     //GWBUG - an if stmt may be faster.  So might "(i64) t < 0".
GF31 OVERLOAD sub(GF31 a, GF31 b) { return U2(sub(a.x, b.x), sub(a.y, b.y)); }

Z31 OVERLOAD neg(Z31 a) { return a == 0 ? 0 : M31 - a; }                // GWBUG: Examine all callers to see if neg call can be avoided
GF31 OVERLOAD neg(GF31 a) { return U2(neg(a.x), neg(a.y)); }

Z31 OVERLOAD make_Z31(i32 a) { return (Z31) (a < 0 ? a + M31 : a); }              // Handles signed values of a
Z31 OVERLOAD make_Z31(u32 a) { return (Z31) (a); }                                // a must be in range of 0 .. M31-1
Z31 OVERLOAD make_Z31(i64 a) { if (a < 0) a = ((i64) M31 << 31) - a; return add((Z31) a, (Z31) (a >> 31)); } // Handles 62-bit a values

u32 get_Z31(Z31 a) { return a; }  // Get balanced value in range 0 to M31-1
i32 get_balanced_Z31(Z31 a) { return (a & 0xC0000000) ? (i32) a - M31 : (i32) a; }  // Get balanced value in range -M31/2 to M31/2

Z31 OVERLOAD mod(Z31 a) { return (a & M31) + (a >> MBITS); }            // GWBUG: This could be larger than M31 (unless a is result of an add), need a wesk and strong mod

// Assumes k reduced mod MBITS.
Z31 OVERLOAD shl(Z31 a, u32 k) { return ((a << k) + (a >> (31 - k))) & M31; }
GF31 OVERLOAD shl(GF31 a, u32 k) { return U2(shl(a.x, k), shl(a.y, k)); }
Z31 OVERLOAD shr(Z31 a, u32 k) { return ((a >> k) + (a << (31 - k))) & M31; }
GF31 OVERLOAD shr(GF31 a, u32 k) { return U2(shr(a.x, k), shr(a.y, k)); }

Z31 OVERLOAD mul(Z31 a, Z31 b) { u64 t = a * (u64) b; return add((Z31) (t & M31), (Z31) (t >> 31)); }

Z31 OVERLOAD fma(Z31 a, Z31 b, Z31 c) { return add(mul(a, b), c); }             // GWBUG:  Can we do better?

// Multiply by 2
Z31 OVERLOAD mul2(Z31 a) { return ((a + a) + (a >> 30)) & M31; }        // GWBUG: Can we do better?
GF31 OVERLOAD mul2(GF31 a) { return U2(mul2(a.x), mul2(a.y)); }

// Return conjugate of a
GF31 OVERLOAD conjugate(GF31 a) { return U2(a.x, neg(a.y)); }

// Complex square.  input, output 31 bits. Uses (a + i*b)^2 == ((a+b)*(a-b) + i*2*a*b).
GF31 OVERLOAD csq(GF31 a) { return U2(mul(add(a.x, a.y), sub(a.x, a.y)), mul2(mul(a.x, a.y))); }        //GWBUG: Probably faster to double a.y and have a mul that takes non-normalized inputs

// a^2 + c
GF31 OVERLOAD csqa(GF31 a, GF31 c) { return add(csq(a), c); }                                           // GWBUG: inline csq so we only "mod" after adding c??  Find a way to use fma instructions

// Complex mul
//GF31 OVERLOAD cmul(GF31 a, GF31 b) { return U2(sub(mul(a.x, b.x), mul(a.y, b.y)), add(mul(a.x, b.y), mul(a.y, b.x)));}   // GWBUG:  Is a 3 multiply complex mul faster?  See above
GF31 OVERLOAD cmul(GF31 a, GF31 b) {
  Z31 k1 = mul(b.x, add(a.x, a.y));
  Z31 k2 = mul(a.x, sub(b.y, b.x));
  Z31 k3 = mul(a.y, add(b.y, b.x));
  return U2(sub(k1, k3), add(k1, k2));
}

GF31 OVERLOAD cfma(GF31 a, GF31 b, GF31 c) { return add(cmul(a, b), c); }                               //GWBUG:  Can we do better?

GF31 OVERLOAD cmul_by_conjugate(GF31 a, GF31 b) { return cmul(a, conjugate(b)); }                       //GWBUG: We can likely eliminate a negate

// Multiply a by b and conjugate(b).  This saves 2 multiplies.
void OVERLOAD cmul_a_by_b_and_conjb(GF31 *res1, GF31 *res2, GF31 a, GF31 b) {
  Z31 axbx = mul(a.x, b.x);
  Z31 aybx = mul(a.y, b.x);
  res1->x = fma(a.y, neg(b.y), axbx), res1->y = fma(a.x,  b.y, aybx);                           //GWBUG: Can we eliminate neg?
  res2->x = fma(a.y,  b.y, axbx), res2->y = fma(a.x, neg(b.y), aybx);                           //GWBUG: Can we eliminate neg?  At least make it a tmp.
}

// mul with (0, 1). (twiddle of tau/4, sqrt(-1) aka "i").
GF31 OVERLOAD mul_t4(GF31 a) { return U2(neg(a.y), a.x); }                                              // GWBUG:  Can caller use a version that does not negate real?

// mul with (2^15, 2^15). (twiddle of tau/8 aka sqrt(i)). Note: 2 * (+/-2^15)^2 == 1 (mod M31).
GF31 OVERLOAD mul_t8(GF31 a) { return U2(shl(sub(a.x, a.y), 15), shl(add(a.x, a.y), 15)); }       // GWBUG:  Can caller use a version that does not negate real?  is shl(neg) same as shr???

// mul with (-2^15, 2^15). (twiddle of 3*tau/8).
GF31 OVERLOAD mul_3t8(GF31 a) { return U2(shl(neg(add(a.x, a.y)), 15), shl(sub(a.x, a.y), 15)); }

// Return a+b and a-b
void OVERLOAD X2_internal(GF31 *a, GF31 *b) { GF31 t = *a; *a = add(t, *b); *b = sub(t, *b); }

// Same as X2(a, conjugate(b))
void OVERLOAD X2conjb_internal(GF31 *a, GF31 *b) { GF31 t = *a; a->x = add(a->x, b->x); a->y = sub(a->y, b->y); b->x = sub(t.x, b->x); b->y = add(t.y, b->y); }

// Same as X2(a, b), b = mul_t4(b)
void OVERLOAD X2_mul_t4_internal(GF31 *a, GF31 *b) { GF31 t = *a; *a = add(*a, *b); t.x = sub(t.x, b->x); b->x = sub(b->y, t.y); b->y = t.x; }

// Same as X2(a, b), a = conjugate(a)
void OVERLOAD X2_conja_internal(GF31 *a, GF31 *b) { GF31 t = *a; a->x = add(a->x, b->x); a->y = neg(add(a->y, b->y)); *b = sub(t, *b); }

// Same as X2(a, b), b = conjugate(b)
void OVERLOAD X2_conjb_internal(GF31 *a, GF31 *b) { GF31 t = *a; *a = add(t, *b); b->x = sub(t.x, b->x); b->y = sub(b->y, t.y); }

void OVERLOAD SWAP_internal(GF31 *a, GF31 *b) { GF31 t = *a; *a = *b; *b = t; }

GF31 OVERLOAD addsub(GF31 a) { return U2(add(a.x, a.y), sub(a.x, a.y)); }
GF31 OVERLOAD foo2(GF31 a, GF31 b) { a = addsub(a); b = addsub(b); return addsub(U2(mul(RE(a), RE(b)), mul(IM(a), IM(b)))); }
GF31 OVERLOAD foo(GF31 a) { return foo2(a, a); }

#endif



/**************************************************************************/
/*          Similar to above, but for an NTT based on GF(M61^2)           */
/**************************************************************************/

#if NTT_GF61

// bits in reduced mod M.
#define MBITS 61
#define M61 ((((Z61) 1) << MBITS) - 1)

Z61 OVERLOAD make_Z61(i32 a) { return (Z61) (a < 0 ? (i64) a + M61 : (i64) a); }  // Handles all values of a
Z61 OVERLOAD make_Z61(i64 a) { return (Z61) (a < 0 ? a + M61 : a); }              // a must be in range of -M61 .. M61-1
Z61 OVERLOAD make_Z61(u64 a) { return (Z61) (a); }                                // a must be in range of 0 .. M61-1

#if 0   // Slower version that keeps results strictly in the range 0 .. M61-1

u64 OVERLOAD get_Z61(Z61 a) { return a; }  // Get value in range 0 to M61-1
i64 OVERLOAD get_balanced_Z61(Z61 a) { return (hi32(a) & 0xF0000000) ? (i64) a - (i64) M61 : (i64) a; }  // Get balanced value in range -M61/2 to M61/2

Z61 OVERLOAD neg(Z61 a) { return a == 0 ? 0 : M61 - a; }                // GWBUG: Examine all callers to see if neg call can be avoided
GF61 OVERLOAD neg(GF61 a) { return U2(neg(a.x), neg(a.y)); }

Z61 OVERLOAD add(Z61 a, Z61 b) { Z61 t = a + b; Z61 m = t - M61; return (m & 0x8000000000000000ULL) ? t : m; }
//Z61 OVERLOAD add(Z61 a, Z61 b) { Z61 t = a + b; Z61 m = t - M61; return t < m ? t : m; }      // Slower on TitanV
//Z61 OVERLOAD add(Z61 a, Z61 b) { Z61 t = a + b; return t - (t >= M61 ? M61 : 0); }            // Slower on TitanV
GF61 OVERLOAD add(GF61 a, GF61 b) { return U2(add(a.x, b.x), add(a.y, b.y)); }

Z61 OVERLOAD sub(Z61 a, Z61 b) { Z61 t = a - b; return t + (((i64) t >> 63) & 0x1FFFFFFFFFFFFFFFULL); }
//Z61 OVERLOAD sub(Z61 a, Z61 b) { Z61 t = a - b; Z61 p = t + M61; return (t & 0x8000000000000000ULL) ? p : t; }  // Better???
//Z61 OVERLOAD sub(Z61 a, Z61 b) { Z61 t = a - b; return t + (t >= M61 ? M61 : 0); }            // Slower on TitanV
//  BETTER???:   t = a - b;  carry_mask = sbb x, x; (generates 32 bits of 0 or 1; return t + make_carry_mask_64_bits
GF61 OVERLOAD sub(GF61 a, GF61 b) { return U2(sub(a.x, b.x), sub(a.y, b.y)); }

// Assumes k reduced mod MBITS.
Z61 OVERLOAD shl(Z61 a, u32 k) { return ((a << k) + (a >> (61 - k))) & M61; }                   //GWBUG: Make sure & M61 operates on just one u32
GF61 OVERLOAD shl(GF61 a, u32 k) { return U2(shl(a.x, k), shl(a.y, k)); }
Z61 OVERLOAD shr(Z61 a, u32 k) { return ((a >> k) + (a << (61 - k))) & M61; }                   //GWBUG: Make sure & M61 operates on just one u32.  & M61 not needed?
GF61 OVERLOAD shr(GF61 a, u32 k) { return U2(shr(a.x, k), shr(a.y, k)); }

ulong2 wideMul(u64 ab, u64 cd) {
  u128 r = (u128) ab * (u128) cd;
  return U2((u64) r, (u64) (r >> 64));
}

Z61 OVERLOAD mul(Z61 a, Z61 b) {
  ulong2 ab = wideMul(a, b);
  u64 lo = ab.x, hi = ab.y;
  u64 lo61 = lo & M61, hi61 = (hi << 3) + (lo >> 61);
  return add(lo61, hi61);
}

Z61 OVERLOAD fma(Z61 a, Z61 b, Z61 c) { return add(mul(a, b), c); }             // GWBUG:  Can we do better?

// Multiply by 2
Z61 OVERLOAD mul2(Z61 a) { return ((a + a) + (a >> 60)) & M61; }        // GWBUG: Make sure "+ a>>60" does an add to lower u32 without a followup adc.
GF61 OVERLOAD mul2(GF61 a) { return U2(mul2(a.x), mul2(a.y)); }

// Return conjugate of a
GF61 OVERLOAD conjugate(GF61 a) { return U2(a.x, neg(a.y)); }

// Complex square.  input, output 61 bits. Uses (a + i*b)^2 == ((a+b)*(a-b) + i*2*a*b).
GF61 OVERLOAD csq(GF61 a) { return U2(mul(add(a.x, a.y), sub(a.x, a.y)), mul2(mul(a.x, a.y))); }        //GWBUG: Probably faster to double a.y and have a mul that takes non-normalized inputs

// a^2 + c
GF61 OVERLOAD csqa(GF61 a, GF61 c) { return add(csq(a), c); }                                           // GWBUG: inline csq so we only "mod" after adding c??  Find a way to use fma instructions
// a^2 - c
GF61 OVERLOAD csqs(GF61 a, GF61 c) { return sub(csq(a), c); }                                           // GWBUG: inline csq so we only "mod" after adding c??  Find a way to use fma instructions

// Complex mul
//GF61 OVERLOAD cmul(GF61 a, GF61 b) { return U2(sub(mul(a.x, b.x), mul(a.y, b.y)), add(mul(a.x, b.y), mul(a.y, b.x)));}   // GWBUG:  Is a 3 multiply complex mul faster?  See above
GF61 OVERLOAD cmul(GF61 a, GF61 b) {
  Z61 k1 = mul(b.x, add(a.x, a.y));
  Z61 k2 = mul(a.x, sub(b.y, b.x));
  Z61 k3 = mul(a.y, add(b.y, b.x));
  return U2(sub(k1, k3), add(k1, k2));
}

GF61 OVERLOAD cfma(GF61 a, GF61 b, GF61 c) { return add(cmul(a, b), c); }                               //GWBUG:  Can we do better?

GF61 OVERLOAD cmul_by_conjugate(GF61 a, GF61 b) { return cmul(a, conjugate(b)); }                       //GWBUG: We can likely eliminate a negate

// Multiply a by b and conjugate(b).  This saves 2 multiplies.
void OVERLOAD cmul_a_by_b_and_conjb(GF61 *res1, GF61 *res2, GF61 a, GF61 b) {
  Z61 axbx = mul(a.x, b.x);
  Z61 aybx = mul(a.y, b.x);
  res1->x = fma(a.y, neg(b.y), axbx), res1->y = fma(a.x,  b.y, aybx);                           //GWBUG: Can we eliminate neg?
  res2->x = fma(a.y,  b.y, axbx), res2->y = fma(a.x, neg(b.y), aybx);                           //GWBUG: Can we eliminate neg?  At least make it a tmp.
}

// mul with (0, 1). (twiddle of tau/4, sqrt(-1) aka "i").
GF61 OVERLOAD mul_t4(GF61 a) { return U2(neg(a.y), a.x); }                                              // GWBUG:  Can caller use a version that does not negate real?

// mul with (-2^30, -2^30). (twiddle of tau/8 aka sqrt(i)). Note: 2 * (+/-2^30)^2 == 1 (mod M61).
GF61 OVERLOAD mul_t8(GF61 a) { return shl(U2(sub(a.y, a.x), neg(add(a.x, a.y))), 30); }       // GWBUG:  Can caller use a version that does not negate real?

// mul with (2^30, -2^30). (twiddle of 3*tau/8).
GF61 OVERLOAD mul_3t8(GF61 a) { return shl(U2(add(a.x, a.y), sub(a.y, a.x)), 30); }

// Return a+b and a-b
void OVERLOAD X2_internal(GF61 *a, GF61 *b) { GF61 t = *a; *a = add(t, *b); *b = sub(t, *b); }

// Same as X2(a, conjugate(b))
void OVERLOAD X2conjb_internal(GF61 *a, GF61 *b) { GF61 t = *a; a->x = add(a->x, b->x); a->y = sub(a->y, b->y); b->x = sub(t.x, b->x); b->y = add(t.y, b->y); }

// Same as X2(a, b), b = mul_t4(b)
void OVERLOAD X2_mul_t4_internal(GF61 *a, GF61 *b) { GF61 t = *a; *a = add(*a, *b); t.x = sub(t.x, b->x); b->x = sub(b->y, t.y); b->y = t.x; }

// Same as X2(a, b), b = mul_t8(b)
void OVERLOAD X2_mul_t8_internal(GF61 *a, GF61 *b) { X2(*a, *b); *b = mul_t8(*b); }

// Same as X2(a, b), b = mul_3t8(b)
void OVERLOAD X2_mul_3t8_internal(GF61 *a, GF61 *b) { X2(*a, *b); *b = mul_3t8(*b); }

// Same as X2(a, b), a = conjugate(a)
void OVERLOAD X2_conja_internal(GF61 *a, GF61 *b) { GF61 t = *a; a->x = add(a->x, b->x); a->y = neg(add(a->y, b->y)); *b = sub(t, *b); }

// Same as X2(a, b), b = conjugate(b)
void OVERLOAD X2_conjb_internal(GF61 *a, GF61 *b) { GF61 t = *a; *a = add(t, *b); b->x = sub(t.x, b->x); b->y = sub(b->y, t.y); }

void OVERLOAD SWAP_internal(GF61 *a, GF61 *b) { GF61 t = *a; *a = *b; *b = t; }

GF61 OVERLOAD addsub(GF61 a) { return U2(add(a.x, a.y), sub(a.x, a.y)); }
GF61 OVERLOAD foo2(GF61 a, GF61 b) { a = addsub(a); b = addsub(b); return addsub(U2(mul(RE(a), RE(b)), mul(IM(a), IM(b)))); }
GF61 OVERLOAD foo(GF61 a) { return foo2(a, a); }





// Philosophy: This Z61/GF61 implementation uses faster, sloppier mod M61 reduction where the end result is in the range 0..M61+epsilon.
// This implementation also handles subtractions by adding enough M61s to make a value positive.  This allows us to always deal with positive
// intermediate results.  The downside is that a caller using the sloppy/quick routines must keep track of how large unreduced values can get.
// An alternative implementation is to have Z61 be an i64 (costs us a precious bit of precision) but is surprisingly slower (at least on TitanV) because
//      mod(a - b),             where the mod routinue uses a signed right shift is slower than
//      mod(a + (M61*2 - b))    where the mod routine uses an unsigned shift right.
// However, a long string of subtracts (example, fft8 does 3 subtracts before mod M61 might be better off using negative intermediate results.
// The mul routine (and obviously csq and cmul) must use only positive values as __int128 multiply is very slow.

#elif 1   // Faster version that keeps results in the range 0 .. M61+epsilon

u64 OVERLOAD get_Z61(Z61 a) { Z61 m = a - M61; return (m & 0x8000000000000000ULL) ? a : m; }  // Get value in range 0 to M61-1
i64 OVERLOAD get_balanced_Z61(Z61 a) { return (hi32(a) & 0xF0000000) ? (i64) a - (i64) M61 : (i64) a; }  // Get balanced value in range -M61/2 to M61/2

// Internal routine to bring Z61 value into the range 0..M61+epsilon
Z61 OVERLOAD mod(Z61 a) { return (a & M61) + (a >> MBITS); }
GF61 OVERLOAD mod(GF61 a) { return U2(mod(a.x), mod(a.y)); }
// Internal routine to negate a value by adding the specified number of M61s
Z61 OVERLOAD neg(Z61 a, u32 m61_count) { return m61_count * M61 - a; }
GF61 OVERLOAD neg(GF61 a, u32 m61_count) { return U2(neg(a.x, m61_count), neg(a.y, m61_count)); }

Z61 OVERLOAD add(Z61 a, Z61 b) { return mod(a + b); }
GF61 OVERLOAD add(GF61 a, GF61 b) { return U2(add(a.x, b.x), add(a.y, b.y)); }

Z61 OVERLOAD sub(Z61 a, Z61 b) { return mod(a + neg(b, 2)); }
GF61 OVERLOAD sub(GF61 a, GF61 b) { return U2(sub(a.x, b.x), sub(a.y, b.y)); }

Z61 OVERLOAD neg(Z61 a) { return mod (neg(a, 2)); }                // GWBUG: Examine all callers to see if neg call can be avoided
GF61 OVERLOAD neg(GF61 a) { return U2(neg(a.x), neg(a.y)); }

// Assumes k reduced mod MBITS.
Z61 OVERLOAD shr(Z61 a, u32 k) { return (a >> k) + ((a << (61 - k)) & M61); }         // Return range 0..M61+2^(61-k), can handle 64-bit inputs but small k is big epsilon
GF61 OVERLOAD shr(GF61 a, u32 k) { return U2(shr(a.x, k), shr(a.y, k)); }
Z61 OVERLOAD shl(Z61 a, u32 k) { return shr(a, 61 - k); }         // Return range 0..M61+2^k, can handle 64-bit inputs but large k yields big epsilon
//Z61 OVERLOAD shl(Z61 a, u32 k) { return mod(a << k) + ((a >> (64 - k)) << 3); }         // Return range 0..M61+2^k, can handle 64-bit inputs but large k is big epsilon
//Z61 OVERLOAD shl(Z61 a, u32 k) { return mod((a << k) + ((a >> (64 - k)) << 3)); }       // Return range 0..M61+epsilon, input must be M61+epsilon a full 62-bit value can overflow
GF61 OVERLOAD shl(GF61 a, u32 k) { return U2(shl(a.x, k), shl(a.y, k)); }

ulong2 wideMul(u64 ab, u64 cd) {
  u128 r = (u128) ab * (u128) cd;
  return U2((u64) r, (u64) (r >> 64));
}

Z61 OVERLOAD weakMul(Z61 a, Z61 b) {                    // a*b must fit in 125 bits, result will as large as a*b >> 61
  ulong2 ab = wideMul(a, b);
  u64 lo = ab.x, hi = ab.y;
  u64 lo61 = lo & M61, hi61 = (hi << 3) + (lo >> 61);
  return lo61 + hi61;
}

Z61 OVERLOAD mul(Z61 a, Z61 b) { return mod(weakMul(a, b)); }

Z61 OVERLOAD fma(Z61 a, Z61 b, Z61 c) { return mod(weakMul(a, b) + c); }             // GWBUG:  Can we do better?

// Multiply by 2
Z61 OVERLOAD mul2(Z61 a) { return add(a, a); }
GF61 OVERLOAD mul2(GF61 a) { return U2(mul2(a.x), mul2(a.y)); }

// Return conjugate of a
GF61 OVERLOAD conjugate(GF61 a) { return U2(a.x, neg(a.y)); }

// Complex square. Uses (a + i*b)^2 == ((a+b)*(a-b) + i*2*a*b).
GF61 OVERLOAD csq(GF61 a) { return U2(mul(a.x + a.y, mod(a.x + neg(a.y, 2))), mul2(weakMul(a.x, a.y))); }

// a^2 + c
GF61 OVERLOAD csqa(GF61 a, GF61 c) { return U2(mod(weakMul(a.x + a.y, mod(a.x + neg(a.y, 2))) + c.x), mod(weakMul(a.x + a.x, a.y) + c.y)); }

// a^2 - c
GF61 OVERLOAD csqs(GF61 a, GF61 c) { return U2(mod(weakMul(a.x + a.y, mod(a.x + neg(a.y, 2))) + neg(c.x, 2)), mod(weakMul(a.x + a.x, a.y) + neg(c.y, 2))); }

// Complex mul
//GF61 OVERLOAD cmul(GF61 a, GF61 b) { return U2(sub(mul(a.x, b.x), mul(a.y, b.y)), add(mul(a.x, b.y), mul(a.y, b.x)));}
GF61 OVERLOAD cmul(GF61 a, GF61 b) {            // Use 2 extra bits in u64
  Z61 k1 = weakMul(b.x, a.x + a.y);             // 61+e * 62+e bits = 123+e mult = 62+e bit result
  Z61 k2 = weakMul(a.x, b.y + neg(b.x, 2));     // 61+e * 63+e bits = 63+e bit result
  Z61 k3 = weakMul(neg(a.y, 2), b.y + b.x);     // 62 * 62+e bits = 63+e bit result
  return U2(mod(k1 + k3), mod(k1 + k2));        // k1+k3 and k1+k2 are full 64-bit values
}

GF61 OVERLOAD cfma(GF61 a, GF61 b, GF61 c) { return add(cmul(a, b), c); }                               //GWBUG:  Can we do better?

GF61 OVERLOAD cmul_by_conjugate(GF61 a, GF61 b) { return cmul(a, conjugate(b)); }                       //GWBUG: We can likely eliminate a negate

// Multiply a by b and conjugate(b).  This saves 2 multiplies.
void OVERLOAD cmul_a_by_b_and_conjb(GF61 *res1, GF61 *res2, GF61 a, GF61 b) {
  Z61 axbx = mul(a.x, b.x);
  Z61 aybx = mul(a.y, b.x);
  res1->x = fma(a.y, neg(b.y), axbx), res1->y = fma(a.x,  b.y, aybx);                           //GWBUG: Can we eliminate neg?
  res2->x = fma(a.y,  b.y, axbx), res2->y = fma(a.x, neg(b.y), aybx);                           //GWBUG: Can we eliminate neg?  At least make it a tmp.
}

// mul with (0, 1). (twiddle of tau/4, sqrt(-1) aka "i").
GF61 OVERLOAD mul_t4(GF61 a) { return U2(neg(a.y), a.x); }                                      // GWBUG:  Can caller use a version that does not negate real?

// mul with (-2^30, -2^30). (twiddle of tau/8 aka sqrt(i)). Note: 2 * (+/-2^30)^2 == 1 (mod M61).
GF61 OVERLOAD mul_t8(GF61 a, u32 m61_count) { return shl(U2(a.y + neg(a.x, m61_count), neg(a.x + a.y, 2 * m61_count - 1)), 30); }
GF61 OVERLOAD mul_t8(GF61 a) { return mul_t8(a, 2); }

// mul with (2^30, -2^30). (twiddle of 3*tau/8).
GF61 OVERLOAD mul_3t8(GF61 a, u32 m61_count) { return shl(U2(a.x + a.y, a.y + neg(a.x, m61_count)), 30); }
GF61 OVERLOAD mul_3t8(GF61 a) { return mul_3t8(a, 2); }

// Return a+b and a-b
void OVERLOAD X2_internal(GF61 *a, GF61 *b) { GF61 t = *a; *a = add(t, *b); *b = sub(t, *b); }

// Same as X2(a, conjugate(b))
void OVERLOAD X2conjb_internal(GF61 *a, GF61 *b) { GF61 t = *a; a->x = add(a->x, b->x); a->y = sub(a->y, b->y); b->x = sub(t.x, b->x); b->y = add(t.y, b->y); }

// Same as X2(a, 2*b), b = conjugate(b)
void OVERLOAD X2mul2_conjb_internal(GF61 *a, GF61 *b) { GF61 t = *b + *b; b->x = sub(a->x, t.x); b->y = sub(t.y, a->y); *a = add(*a, t); }

// Same as X2(a, b), b = mul_t4(b)
void OVERLOAD X2_mul_t4_internal(GF61 *a, GF61 *b) { GF61 t = *a; *a = add(*a, *b); t.x = sub(t.x, b->x); b->x = sub(b->y, t.y); b->y = t.x; }

// Same as X2(a, b), b = mul_t8(b)
void OVERLOAD X2_mul_t8_internal(GF61 *a, GF61 *b) { GF61 t = *a; *a = add(t, *b); *b = t + neg(*b, 2); *b = mul_t8(*b, 4); }

// Same as X2(a, b), b = mul_3t8(b)
void OVERLOAD X2_mul_3t8_internal(GF61 *a, GF61 *b) { GF61 t = *a; *a = add(t, *b); *b = t + neg(*b, 2); *b = mul_3t8(*b, 4); }

// Same as X2(a, b), a = conjugate(a)
void OVERLOAD X2_conja_internal(GF61 *a, GF61 *b) { GF61 t = *a; a->x = add(a->x, b->x); a->y = neg(add(a->y, b->y)); *b = sub(t, *b); }

// Same as X2(a, b), b = conjugate(b)
void OVERLOAD X2_conjb_internal(GF61 *a, GF61 *b) { GF61 t = *a; *a = add(t, *b); b->x = sub(t.x, b->x); b->y = sub(b->y, t.y); }

void OVERLOAD SWAP_internal(GF61 *a, GF61 *b) { GF61 t = *a; *a = *b; *b = t; }

GF61 OVERLOAD addsub(GF61 a) { return U2(add(a.x, a.y), sub(a.x, a.y)); }
GF61 OVERLOAD foo2(GF61 a, GF61 b) { a = addsub(a); b = addsub(b); return addsub(U2(mul(RE(a), RE(b)), mul(IM(a), IM(b)))); }
GF61 OVERLOAD foo(GF61 a) { return foo2(a, a); }

// The following routines can be used to reduce mod M61 operations.  Caller must track how many M61s need to be added to make positive
// values for subtractions.  In X2 functions, "q" stands for quick, "s" stands for slow.

Z61 OVERLOAD sub(Z61 a, Z61 b, u32 m61_count) { return mod(a + neg(b, m61_count)); }
GF61 OVERLOAD sub(GF61 a, GF61 b, u32 m61_count) { return U2(sub(a.x, b.x, m61_count), sub(a.y, b.y, m61_count)); }

void OVERLOAD X2s(GF61 *a, GF61 *b, u32 m61_count) { GF61 t = *a; *a = add(t, *b); *b = sub(t, *b, m61_count); }
void OVERLOAD X2q(GF61 *a, GF61 *b, u32 m61_count) { GF61 t = *a; *a = t + *b; *b = t + neg(*b, m61_count); }
void OVERLOAD X2q_mul_t4(GF61 *a, GF61 *b, u32 m61_count) { GF61 t = *a; *a = t + *b; t.x = t.x + neg(b->x, m61_count); b->x = b->y + neg(t.y, m61_count); b->y = t.x; }

















#elif 0   // Abandoned version that returned results in the range 0 .. M61+epsilon,  intermediate results are signed (tyecast to i64).  Z61 data type is u64.

u64 OVERLOAD get_Z61(Z61 a) { Z61 m = a - M61; return (m & 0x8000000000000000ULL) ? a : m; }  // Get value in range 0 to M61-1
i64 OVERLOAD get_balanced_Z61(Z61 a) { return (hi32(a) & 0xF0000000) ? (i64) a - (i64) M61 : (i64) a; }  // Get balanced value in range -M61/2 to M61/2

Z61 OVERLOAD mod(Z61 a) { return (a & M61) + (a >> MBITS); }            // return a in range 0..M61+epsilon
Z61 OVERLOAD mod(i64 a) { return (a & M61) + (a >> MBITS); }            // return a in range 0..M61+epsilon

Z61 OVERLOAD add(Z61 a, Z61 b) { return mod(a + b); }
GF61 OVERLOAD add(GF61 a, GF61 b) { return U2(add(a.x, b.x), add(a.y, b.y)); }

Z61 OVERLOAD sub(Z61 a, Z61 b) { return mod((i64) (a - b)); }
GF61 OVERLOAD sub(GF61 a, GF61 b) { return U2(sub(a.x, b.x), sub(a.y, b.y)); }

Z61 OVERLOAD neg(Z61 a) { return mod (- (i64) a); }                // GWBUG: Examine all callers to see if neg call can be avoided
GF61 OVERLOAD neg(GF61 a) { return U2(neg(a.x), neg(a.y)); }

// Assumes k reduced mod MBITS.
Z61 OVERLOAD shr(Z61 a, u32 k) { return (a >> k) + ((a << (61 - k)) & M61); }
GF61 OVERLOAD shr(GF61 a, u32 k) { return U2(shr(a.x, k), shr(a.y, k)); }
Z61 OVERLOAD shl(Z61 a, u32 k) { return shr(a, 61 - k); }         // Return range 0..M61+2^k, can handle 64-bit inputs but large k yields big epsilon
//Z61 OVERLOAD shl(Z61 a, u32 k) { return mod(a << k) + ((a >> (64 - k)) << 3); }         // Return range 0..M61+2^k, can handle 64-bit inputs but large k is big epsilon
//Z61 OVERLOAD shl(Z61 a, u32 k) { return mod((a << k) + ((a >> (64 - k)) << 3)); }       // Return range 0..M61+epsilon, input must be M61+epsilon a full 62-bit value can overflow
GF61 OVERLOAD shl(GF61 a, u32 k) { return U2(shl(a.x, k), shl(a.y, k)); }

ulong2 wideMul(u64 ab, u64 cd) {
  u128 r = (u128) ab * (u128) cd;
  return U2((u64) r, (u64) (r >> 64));
}

Z61 OVERLOAD weakMul(Z61 a, Z61 b) {                    // a*b must fit in 125 bits, result will as large as a*b >> 61
  ulong2 ab = wideMul(a, b);
  u64 lo = ab.x, hi = ab.y;
  u64 lo61 = lo & M61, hi61 = (hi << 3) + (lo >> 61);
  return lo61 + hi61;
}

Z61 OVERLOAD mul(Z61 a, Z61 b) { return mod(weakMul(a, b)); }

Z61 OVERLOAD fma(Z61 a, Z61 b, Z61 c) { return mod(weakMul(a, b) + c); }             // GWBUG:  Can we do better?

// Multiply by 2
Z61 OVERLOAD mul2(Z61 a) { return add(a, a); }
GF61 OVERLOAD mul2(GF61 a) { return U2(mul2(a.x), mul2(a.y)); }

// Return conjugate of a
GF61 OVERLOAD conjugate(GF61 a) { return U2(a.x, neg(a.y)); }

// Complex square. Uses (a + i*b)^2 == ((a+b)*(a-b) + i*2*a*b).
GF61 OVERLOAD csq(GF61 a) { return U2(mul(mod(a.x + a.y), mod(a.x - a.y + M61 * 2)), mul2(weakMul(a.x, a.y))); }

// a^2 + c
GF61 OVERLOAD csqa(GF61 a, GF61 c) { return add(csq(a), c); }                                           // GWBUG: inline csq so we only "mod" after adding c??  Find a way to use fma instructions

// Complex mul
//GF61 OVERLOAD cmul(GF61 a, GF61 b) { return U2(sub(mul(a.x, b.x), mul(a.y, b.y)), add(mul(a.x, b.y), mul(a.y, b.x)));}   // GWBUG:  Is a 3 multiply complex mul faster?  See above
//GF61 OVERLOAD cmul(GF61 a, GF61 b) {
//  Z61 k1 = mul(b.x, add(a.x, a.y));
//  Z61 k2 = mul(a.x, sub(b.y, b.x));
//  Z61 k3 = mul(a.y, add(b.y, b.x));
//  return U2(sub(k1, k3), add(k1, k2));
//}
GF61 OVERLOAD cmul(GF61 a, GF61 b) {            // Use 2 extra bits in u64
  Z61 k1 = weakMul(b.x, a.x + a.y);             // 61 * 62 bits = 62+e
  Z61 k2 = weakMul(a.x, mod(b.y - b.x + M61 * 2));   // 61 * 62 bits = 62+e
  Z61 k3 = weakMul(mod(M61 * 2 - a.y), b.y + b.x);   // 62 * 62 bits = 63+e
  return U2(mod(k1 + k3), mod(k1 + k2));
}

GF61 OVERLOAD cfma(GF61 a, GF61 b, GF61 c) { return add(cmul(a, b), c); }                               //GWBUG:  Can we do better?

GF61 OVERLOAD cmul_by_conjugate(GF61 a, GF61 b) { return cmul(a, conjugate(b)); }                       //GWBUG: We can likely eliminate a negate

// Multiply a by b and conjugate(b).  This saves 2 multiplies.
void OVERLOAD cmul_a_by_b_and_conjb(GF61 *res1, GF61 *res2, GF61 a, GF61 b) {
  Z61 axbx = mul(a.x, b.x);
  Z61 aybx = mul(a.y, b.x);
  res1->x = fma(a.y, neg(b.y), axbx), res1->y = fma(a.x,  b.y, aybx);                           //GWBUG: Can we eliminate neg?
  res2->x = fma(a.y,  b.y, axbx), res2->y = fma(a.x, neg(b.y), aybx);                           //GWBUG: Can we eliminate neg?  At least make it a tmp.
}

// mul with (0, 1). (twiddle of tau/4, sqrt(-1) aka "i").
GF61 OVERLOAD mul_t4(GF61 a) { return U2(neg(a.y), a.x); }                                              // GWBUG:  Can caller use a version that does not negate real?

// mul with (-2^30, -2^30). (twiddle of tau/8 aka sqrt(i)). Note: 2 * (+/-2^30)^2 == 1 (mod M61).
GF61 OVERLOAD mul_t8(GF61 a) { return U2(shl(a.y - a.x + M61 * 2, 30), shl(M61 * 3 - a.x - a.y, 30)); }       // GWBUG:  Can caller use a version that does not negate real?  

// mul with (2^30, -2^30). (twiddle of 3*tau/8).
GF61 OVERLOAD mul_3t8(GF61 a) { return U2(shl(a.x + a.y, 30), shl(a.y - a.x + M61 * 2, 30)); }

// Return a+b and a-b
void OVERLOAD X2_internal(GF61 *a, GF61 *b) { GF61 t = *a; *a = add(t, *b); *b = sub(t, *b); }

// Same as X2(a, conjugate(b))
void OVERLOAD X2conjb_internal(GF61 *a, GF61 *b) { GF61 t = *a; a->x = add(a->x, b->x); a->y = sub(a->y, b->y); b->x = sub(t.x, b->x); b->y = add(t.y, b->y); }

// Same as X2(a, b), b = mul_t4(b)
void OVERLOAD X2_mul_t4_internal(GF61 *a, GF61 *b) { GF61 t = *a; *a = add(*a, *b); t.x = sub(t.x, b->x); b->x = sub(b->y, t.y); b->y = t.x; }

// Same as X2(a, b), a = conjugate(a)
void OVERLOAD X2_conja_internal(GF61 *a, GF61 *b) { GF61 t = *a; a->x = add(a->x, b->x); a->y = neg(add(a->y, b->y)); *b = sub(t, *b); }

// Same as X2(a, b), b = conjugate(b)
void OVERLOAD X2_conjb_internal(GF61 *a, GF61 *b) { GF61 t = *a; *a = add(t, *b); b->x = sub(t.x, b->x); b->y = sub(b->y, t.y); }

void OVERLOAD SWAP_internal(GF61 *a, GF61 *b) { GF61 t = *a; *a = *b; *b = t; }

GF61 OVERLOAD addsub(GF61 a) { return U2(add(a.x, a.y), sub(a.x, a.y)); }
GF61 OVERLOAD foo2(GF61 a, GF61 b) { a = addsub(a); b = addsub(b); return addsub(U2(mul(RE(a), RE(b)), mul(IM(a), IM(b)))); }
GF61 OVERLOAD foo(GF61 a) { return foo2(a, a); }






#else   // Abandoned start of a version that returns results in the range 0..M61+epsilon but has signed interim results.  Make Z61 data type i64.

Z61 OVERLOAD mod(Z61 a) { return (a & M61) + (a >> MBITS); }            // return a in range 0..M61+epsilon

u64 OVERLOAD get_Z61(Z61 a) { Z61 m = a - M61; return (m & 0x8000000000000000ULL) ? a : m; }  // Get value in range 0 to M61-1
i64 OVERLOAD get_balanced_Z61(Z61 a) { return (hi32(a) & 0xF0000000) ? (i64) a - (i64) M61 : (i64) a; }  // Get balanced value in range -M61/2 to M61/2

Z61 OVERLOAD add(Z61 a, Z61 b) { return mod(a + b); }
GF61 OVERLOAD add(GF61 a, GF61 b) { return U2(add(a.x, b.x), add(a.y, b.y)); }

Z61 OVERLOAD sub(Z61 a, Z61 b) { return mod(a - b); }
GF61 OVERLOAD sub(GF61 a, GF61 b) { return U2(sub(a.x, b.x), sub(a.y, b.y)); }

Z61 OVERLOAD neg(Z61 a) { return mod(-a); }                // GWBUG: Examine all callers to see if neg call can be avoided
GF61 OVERLOAD neg(GF61 a) { return U2(neg(a.x), neg(a.y)); }

// Assumes k reduced mod MBITS.
Z61 OVERLOAD shr(Z61 a, u32 k) { return (a >> k) + ((a << (61 - k)) & M61); }  // Return range 0..M61+2^(61-k), can handle 64-bit inputs but small k is bad for epsilon
GF61 OVERLOAD shr(GF61 a, u32 k) { return U2(shr(a.x, k), shr(a.y, k)); }
Z61 OVERLOAD shl(Z61 a, u32 k) { return shr(a, 61 - k); }         // Return range 0..M61+2^k, can handle 64-bit inputs but large k is big epsilon
//Z61 OVERLOAD shl(Z61 a, u32 k) { return mod(a << k) + ((a >> (64 - k)) << 3); }         // Return range 0..M61+2^k, can handle 64-bit inputs but large k is big epsilon
//Z61 OVERLOAD shl(Z61 a, u32 k) { return mod((a << k) + ((a >> (64 - k)) << 3)); }       // Return range 0..M61+epsilon, input must be M61+epsilon a full 62-bit value can overflow
GF61 OVERLOAD shl(GF61 a, u32 k) { return U2(shl(a.x, k), shl(a.y, k)); }

//GWBUG - need version that handles shifted 3 bits of high?
//GWBUG - need version mods a,b?  CALLER MUST mod possible negative numbers 
Z61 OVERLOAD weakMul(Z61 a, Z61 b) {                    // a*b must fit in 125 bits, result will as large as a*b >> 61
  u128 r = (u128) a * (u128) b;
  u64 lo = (u64) r;
  u64 hi = (u64) (r >> 64);
  return (lo & M61) + ((hi << 3) | (lo >> 61));
}

Z61 OVERLOAD mul(Z61 a, Z61 b) { return mod(weakMul(a, b)); }

Z61 OVERLOAD fma(Z61 a, Z61 b, Z61 c) { return mod(weakMul(a, b) + c); }             // GWBUG:  Can we do better?

// Multiply by 2
Z61 OVERLOAD mul2(Z61 a) { return add(a, a); }
GF61 OVERLOAD mul2(GF61 a) { return U2(mul2(a.x), mul2(a.y)); }

// Return conjugate of a
GF61 OVERLOAD conjugate(GF61 a) { return U2(a.x, neg(a.y)); }

// Complex square. Uses (a + i*b)^2 == ((a+b)*(a-b) + i*2*a*b).
//GWBUG??: GF61 OVERLOAD csq(GF61 a) { return U2(mul(a.x + a.y, a.x - a.y), mul(a.x + a.x, a.y)); }
GF61 OVERLOAD csq(GF61 a) { return U2(mul(mod(a.x + a.y), mod(a.x - a.y)), mul2(weakMul(a.x, a.y))); }          //GWBUG: study mods

// a^2 + c
GF61 OVERLOAD csqa(GF61 a, GF61 c) { return add(csq(a), c); }                                           // GWBUG: inline csq so we only "mod" after adding c??  Find a way to use fma instructions

// Complex mul
//GF61 OVERLOAD cmul(GF61 a, GF61 b) { return U2(sub(mul(a.x, b.x), mul(a.y, b.y)), add(mul(a.x, b.y), mul(a.y, b.x)));}   // GWBUG:  Is a 3 multiply complex mul faster?  See above
GF61 OVERLOAD cmul(GF61 a, GF61 b) {            // Use 2 extra bits in i64 
  Z61 k1 = weakMul(b.x, mod(a.x + a.y));             // 61 * 62 bits => 62+e                    MODS need study!!
  Z61 k2 = weakMul(a.x, mod(b.y - b.x));        // 61 * 61 bits => 62+e
  Z61 k3 = weakMul(a.y, mod(b.y + b.x));             // 61 * 62 bits => 63+e
  return U2(mod(k1 - k3), mod(k1 + k2));
}

GF61 OVERLOAD cfma(GF61 a, GF61 b, GF61 c) { return add(cmul(a, b), c); }                               //GWBUG:  Can we do better?

GF61 OVERLOAD cmul_by_conjugate(GF61 a, GF61 b) { return cmul(a, conjugate(b)); }                       //GWBUG: We can likely eliminate a negate

// Multiply a by b and conjugate(b).  This saves 2 multiplies.
void OVERLOAD cmul_a_by_b_and_conjb(GF61 *res1, GF61 *res2, GF61 a, GF61 b) {
  Z61 axbx = mul(a.x, b.x);
  Z61 aybx = mul(a.y, b.x);
  res1->x = fma(a.y, neg(b.y), axbx), res1->y = fma(a.x,  b.y, aybx);                           //GWBUG: Can we eliminate neg?
  res2->x = fma(a.y,  b.y, axbx), res2->y = fma(a.x, neg(b.y), aybx);                           //GWBUG: Can we eliminate neg?  At least make it a tmp.
}

// mul with (0, 1). (twiddle of tau/4, sqrt(-1) aka "i").
GF61 OVERLOAD mul_t4(GF61 a) { return U2(neg(a.y), a.x); }                                              // GWBUG:  Can caller use a version that does not negate real?

// mul with (-2^30, -2^30). (twiddle of tau/8 aka sqrt(i)). Note: 2 * (+/-2^30)^2 == 1 (mod M61).
GF61 OVERLOAD mul_t8(GF61 a) { return U2(shl(a.y - a.x, 30), shl(-(a.x + a.y), 30)); }       // GWBUG:  Can caller use a version that does not negate real?  

// mul with (2^30, -2^30). (twiddle of 3*tau/8).
GF61 OVERLOAD mul_3t8(GF61 a) { return U2(shl(a.x + a.y, 30), shl(a.y - a.x, 30)); }

// Return a+b and a-b
void OVERLOAD X2_internal(GF61 *a, GF61 *b) { GF61 t = *a; *a = add(t, *b); *b = sub(t, *b); }

// Same as X2(a, conjugate(b))
void OVERLOAD X2conjb_internal(GF61 *a, GF61 *b) { GF61 t = *a; a->x = add(a->x, b->x); a->y = sub(a->y, b->y); b->x = sub(t.x, b->x); b->y = add(t.y, b->y); }

// Same as X2(a, b), b = mul_t4(b)
void OVERLOAD X2_mul_t4_internal(GF61 *a, GF61 *b) { GF61 t = *a; *a = add(*a, *b); t.x = sub(t.x, b->x); b->x = sub(b->y, t.y); b->y = t.x; }

// Same as X2(a, b), a = conjugate(a)
void OVERLOAD X2_conja_internal(GF61 *a, GF61 *b) { GF61 t = *a; a->x = add(a->x, b->x); a->y = neg(add(a->y, b->y)); *b = sub(t, *b); }

// Same as X2(a, b), b = conjugate(b)
void OVERLOAD X2_conjb_internal(GF61 *a, GF61 *b) { GF61 t = *a; *a = add(t, *b); b->x = sub(t.x, b->x); b->y = sub(b->y, t.y); }

void OVERLOAD SWAP_internal(GF61 *a, GF61 *b) { GF61 t = *a; *a = *b; *b = t; }

GF61 OVERLOAD addsub(GF61 a) { return U2(add(a.x, a.y), sub(a.x, a.y)); }
GF61 OVERLOAD foo2(GF61 a, GF61 b) { a = addsub(a); b = addsub(b); return addsub(U2(mul(RE(a), RE(b)), mul(IM(a), IM(b)))); }
GF61 OVERLOAD foo(GF61 a) { return foo2(a, a); }

#endif

#endif
