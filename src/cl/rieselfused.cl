// Helpers for an FFT31R2 fused carry kernel.  Unlike RIESEL_FIELD builds,
// both Riesel fields are live in this translation unit, so every operation
// takes its modulus explicitly and residue planes remain independent uint2s.

#pragma once

#if FFT_TYPE == FFT31R2

#define RIESEL0_R2 1130207172u
#define RIESEL1_R2 1409360092u
#define RIESEL0_T8 1657223104u
#define RIESEL1_T8 240631641u

u32 OVERLOAD rfAdd(u32 a, u32 b, u32 q) {
  u32 const value = a + b;
  return value >= q ? value - q : value;
}

u32 OVERLOAD rfSub(u32 a, u32 b, u32 q) {
  return a >= b ? a - b : q - (b - a);
}

u32 rfNeg(u32 a, u32 q) { return a == 0 ? 0 : q - a; }

GF31 OVERLOAD rfAdd(GF31 a, GF31 b, u32 q) {
  return U2(rfAdd(a.x, b.x, q), rfAdd(a.y, b.y, q));
}

GF31 OVERLOAD rfSub(GF31 a, GF31 b, u32 q) {
  return U2(rfSub(a.x, b.x, q), rfSub(a.y, b.y, q));
}

GF31 rfCmul(GF31 a, GF31 b, u32 q, u32 negInv) {
  u32 const k1 = rieselMontMulExplicit(b.x, rfAdd(a.x, a.y, q), q, negInv);
  u32 const k2 = rieselMontMulExplicit(a.x, rfSub(b.y, b.x, q), q, negInv);
  u32 const k3 = rieselMontMulExplicit(a.y, rfAdd(b.y, b.x, q), q, negInv);
  return U2(rfSub(k1, k3, q), rfAdd(k1, k2, q));
}

GF31 rfMulT4(GF31 a, u32 q) { return U2(rfNeg(a.y, q), a.x); }

GF31 rfMulT8(GF31 a, u32 t8, u32 q, u32 negInv) {
  return rfCmul(a, U2(t8, t8), q, negInv);
}

GF31 rfMul3T8(GF31 a, u32 t8, u32 q, u32 negInv) {
  return rfCmul(a, U2(rfNeg(t8, q), t8), q, negInv);
}

void rfX2(GF31 *a, GF31 *b, u32 q) {
  GF31 const old = *a;
  *a = rfAdd(old, *b, q);
  *b = rfSub(old, *b, q);
}

void rfX2MulT4(GF31 *a, GF31 *b, u32 q) {
  rfX2(a, b, q);
  *b = rfMulT4(*b, q);
}

void rfX2MulT8(GF31 *a, GF31 *b, u32 t8, u32 q, u32 negInv) {
  rfX2(a, b, q);
  *b = rfMulT8(*b, t8, q, negInv);
}

void rfX2Mul3T8(GF31 *a, GF31 *b, u32 t8, u32 q, u32 negInv) {
  rfX2(a, b, q);
  *b = rfMul3T8(*b, t8, q, negInv);
}

void rfFft4Core(GF31 *u, u32 q) {
  rfX2(&u[0], &u[2], q);
  rfX2MulT4(&u[1], &u[3], q);
  rfX2(&u[0], &u[1], q);
  rfX2(&u[2], &u[3], q);
}

void rfFft4(GF31 *u, u32 q) {
  rfFft4Core(u, q);
  GF31 const tmp = u[1]; u[1] = u[2]; u[2] = tmp;
}

void rfFft8(GF31 *u, u32 t8, u32 q, u32 negInv) {
  rfX2(&u[0], &u[4], q);
  rfX2MulT8(&u[1], &u[5], t8, q, negInv);
  rfX2MulT4(&u[2], &u[6], q);
  rfX2Mul3T8(&u[3], &u[7], t8, q, negInv);
  rfFft4Core(u, q);
  rfFft4Core(u + 4, q);
  GF31 tmp = u[1]; u[1] = u[4]; u[4] = tmp;
  tmp = u[3]; u[3] = u[6]; u[6] = tmp;
}

void rfFftRadix(GF31 *u, u32 t8, u32 q, u32 negInv) {
#if NW == 4
  rfFft4(u, q);
#elif NW == 8
  rfFft8(u, t8, q, negInv);
#else
#error FFT31R2 fused carry supports NW=4 or NW=8
#endif
}

void rfTabMul(TrigGF31 trig, GF31 *u, u32 f, u32 me,
              u32 q, u32 negInv) {
  u32 const p = me & ~(f - 1);
#if TABMUL_CHAIN31
  GF31 const w = TFLOAD(&trig[p]);
  GF31 power = w;
  for (u32 i = 1; i < NW; ++i) {
    u[i] = rfCmul(u[i], power, q, negInv);
    power = rfCmul(power, w, q, negInv);
  }
#else
  for (u32 i = 1; i < NW; ++i) {
    u[i] = rfCmul(u[i], TFLOAD(&trig[(i - 1) * G_W + p]), q, negInv);
  }
#endif
}

void rfFftWidth(local GF31 *lds, GF31 *u, TrigGF31 trig,
                u32 numWG, u32 lowMe, u32 t8, u32 q, u32 negInv) {
#if !UNROLL_W
  __attribute__((opencl_unroll_hint(1)))
#endif
  for (u32 s = 1; s < G_W; s *= NW) {
    rfFftRadix(u, t8, q, negInv);
    rfTabMul(trig, u, s, lowMe, q, negInv);
    shufl32((local F2 *)lds, (F2 *)u, s, numWG, lowMe);
  }
  rfFftRadix(u, t8, q, negInv);
}

u32 rfMakeWord(i64 a, u32 q, u32 negInv, u32 r2) {
  bool const negative = a < 0;
  u64 value = negative ? (u64)(-(a + 1)) + 1 : (u64)a;
#if EXP / NWORDS <= 32
  if (value >= q) value -= q;
  if (value >= q) value -= q;
  if (value >= q) value -= q;
  if (value >= q) value -= q;
#else
  value %= q;
#endif
  u32 normal = (u32)value;
  if (negative && normal != 0) normal = q - normal;
  return rieselMontMulExplicit(normal, r2, q, negInv);
}

GF31 rfReadCarryValue(CP(GF31) in, u32 line, u32 me, u32 i) {
#if INPLACE
  u32 const middle = line / SMALL_HEIGHT;
  line %= SMALL_HEIGHT;
  in += (me / 16 * SIZEW32) + (middle * SIZEM32) +
        (line % 16 * SIZEBLK32) + SWIZ32(line % 16, line / 16) * 16 +
        (me % 16);
  return FFTLOAD(&in[i * G_W / 16 * SIZEW32]);
#else
  u32 const sizeY = OUT_WG / OUT_SIZEX;
  u32 const fftMiddleOutX = line % SMALL_HEIGHT;
  u32 const chunkX = fftMiddleOutX / OUT_SIZEX;
#if PAD_SIZE > 0
  u32 const bigPadSize = (PAD_SIZE / 2 + 1) * PAD_SIZE;
  in += chunkX * (WIDTH * MIDDLE * OUT_SIZEX +
                  WIDTH / sizeY * PAD_SIZE + bigPadSize);
#else
  in += chunkX * MIDDLE * WIDTH * OUT_SIZEX;
#endif
  in += (fftMiddleOutX % OUT_SIZEX) * sizeY;
  in += (line / SMALL_HEIGHT) * OUT_WG;
  in += me % sizeY;
  u32 const chunkY = me / sizeY + i * (G_W / sizeY);
#if PAD_SIZE > 0
  return FFTLOAD(&in[chunkY * (MIDDLE * OUT_WG + PAD_SIZE)]);
#else
  return FFTLOAD(&in[chunkY * MIDDLE * OUT_WG]);
#endif
#endif
}

#endif
