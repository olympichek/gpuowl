// Reconstruction of the lost q7/M61 radix-7 gates in ONE process with
// matched alternation (2026-08-12; audit shortlist item 6).  Sol's rejection
// subtracted numbers from two different benchmarks: binary-core saving
// 18.67-20.30 us (m31_m61_q15_overlap_bench) vs radix-7 edge increment
// 20.128-20.288 us (m61_q7_radix7_bench) -- a wash decided across
// instruments, with the candidate's 12.5% smaller state priced only via the
// core arm and the edge's co-run hideability never measured.
//
// This bench measures, in one process, q65-protocol alternated:
//   S      = co-run core saving: (M31ch7|M61ch7 @2^21) - (q7ch6|M61ch6 @7*2^18)
//   E_seq  = edge increment, Sol's mode: q7edge+M61edge run back-to-back on
//            one stream, minus the same-mode square-only controls
//   E_conc = edge increment with the two field edge kernels CO-RUN on two
//            streams (as a production schedule would overlap them)
// Verdict: E_conc vs S.  The direct DFT-7 here overcharges the edge vs
// Sol's Rader/Good-Thomas graph; calibrate with E_seq against Sol's
// clock-scaled 20.2 us before reading E_conc.
//
//   q7 = 4,151,312,383 = 15836*2^18 - 1, R=2^32 Montgomery:
//   -q^-1 = 0xf7700001, R mod q = 143654913, R^-1 = 4012462336.
// Exactness: full-array host oracles for every q7 kernel and for the edge
// kernels of both fields (cyclic-7 self-convolution); 1-element oracles for
// the M31/M61 chains.

#include <cuda_runtime.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <random>
#include <string>
#include <vector>

namespace {

using u32 = std::uint32_t;
using u64 = std::uint64_t;
using u128 = unsigned __int128;

constexpr u32 M31 = 0x7fffffffu;
constexpr u64 M61 = (u64{1} << 61) - 1;
constexpr u32 Q7 = 4151312383u;               // 15836*2^18 - 1
constexpr u32 Q7_NEG_INV = 0xf7700001u;       // -q^-1 mod 2^32
constexpr u32 Q7_R = 143654913u;              // 2^32 mod q
constexpr u32 Q7_R_INV = 4012462336u;         // (2^32)^-1 mod q
constexpr u32 CONTROL_COUNT = 1u << 21;       // 4M architecture population
constexpr u32 CAND_COUNT = 7u << 18;          // 3.5M architecture population
constexpr u32 EDGE_GROUPS = CAND_COUNT / 7;   // 2^18 groups of 7
constexpr int SAMPLES = 21, REPEATS = 20;

static_assert(u32(u64(Q7) * Q7_NEG_INV) == 0xffffffffu);
static_assert(u32(u128(Q7_R) * Q7_R_INV % Q7) == 1);

struct GF31 { u32 x, y; };
struct GF61 { u64 x, y; };
struct GQ7 { u32 x, y; };

void check(cudaError_t result, char const* expression, int line) {
  if (result == cudaSuccess) return;
  std::fprintf(stderr, "%s:%d: %s\n", expression, line, cudaGetErrorString(result));
  std::exit(1);
}
#define CUDA_CHECK(expression) check((expression), #expression, __LINE__)

// ---------------- M31 / M61 (as in the surviving template) ----------------
__device__ __forceinline__ u32 add31(u32 a, u32 b) {
  u32 r = a + b; r = (r & M31) + (r >> 31); return r == M31 ? 0 : r;
}
__device__ __forceinline__ u32 sub31(u32 a, u32 b) { return a >= b ? a - b : a + M31 - b; }
__device__ __forceinline__ u32 mul31(u32 a, u32 b) {
  u64 const p = u64(a) * b;
  u64 r = (p & M31) + (p >> 31);
  r = (r & M31) + (r >> 31);
  return r == M31 ? 0 : u32(r);
}
__device__ __forceinline__ GF31 cmul31(GF31 a, GF31 b) {
  u32 const k1 = mul31(b.x, add31(a.x, a.y));
  u32 const k2 = mul31(a.x, sub31(b.y, b.x));
  u32 const k3 = mul31(a.y, add31(b.y, b.x));
  return {sub31(k1, k3), add31(k1, k2)};
}
__device__ __forceinline__ u64 add61(u64 a, u64 b) {
  u64 r = a + b; return r >= M61 ? r - M61 : r;
}
__device__ __forceinline__ u64 sub61(u64 a, u64 b) { return a >= b ? a - b : M61 - (b - a); }
__device__ __forceinline__ u64 mul61(u64 a, u64 b) {
  u64 const lo = a * b;
  u64 const hi = __umul64hi(a, b);
  u64 r = (lo & M61) + (lo >> 61) + (hi << 3);
  if (r >= M61) r -= M61;
  if (r >= M61) r -= M61;
  return r;
}
__device__ __forceinline__ GF61 cmul61(GF61 a, GF61 b) {
  u64 const k1 = mul61(b.x, add61(a.x, a.y));
  u64 const k2 = mul61(a.x, sub61(b.y, b.x));
  u64 const k3 = mul61(a.y, add61(b.y, b.x));
  return {sub61(k1, k3), add61(k1, k2)};
}
__device__ __forceinline__ GF61 cscale61(GF61 a, u64 s) { return {mul61(a.x, s), mul61(a.y, s)}; }

// ---------------- q7: 32-bit canonical Montgomery, R=2^32 ----------------
// REDC keeps the carry out of the cancelled 65-bit low sum (lo != 0).
__device__ __forceinline__ u32 addQ7(u32 a, u32 b) {
  u32 const r = a + b;
  return (r >= Q7 || r < a) ? r - Q7 : r;
}
__device__ __forceinline__ u32 subQ7(u32 a, u32 b) { return a >= b ? a - b : a + Q7 - b; }
__device__ __forceinline__ u32 mulQ7(u32 a, u32 b) {
  u64 const t = u64(a) * b;
  u32 const m = u32(t) * Q7_NEG_INV;
  u64 const mq = u64(m) * Q7;
  // The cancelled low sum is 65 bits; fold in u64 (result < 2*Q7 > 2^32).
  u64 r = (t >> 32) + (mq >> 32) + (u32(t) != 0);
  if (r >= Q7) r -= Q7;
  return u32(r);
}
__device__ __forceinline__ GQ7 cmulQ7(GQ7 a, GQ7 b) {
  u32 const k1 = mulQ7(b.x, addQ7(a.x, a.y));
  u32 const k2 = mulQ7(a.x, subQ7(b.y, b.x));
  u32 const k3 = mulQ7(a.y, addQ7(b.y, b.x));
  return {subQ7(k1, k3), addQ7(k1, k2)};
}
__device__ __forceinline__ GQ7 cscaleQ7(GQ7 a, u32 s) { return {mulQ7(a.x, s), mulQ7(a.y, s)}; }

// ---------------- chain kernels (binary-core arms) ----------------
template<int Rounds>
__global__ void m31Chain(GF31 const* input, GF31 const* roots, GF31* output, u32 count) {
  u32 const i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= count) return;
  GF31 v = input[i], r = roots[i];
#pragma unroll
  for (int round = 0; round != Rounds; ++round) v = cmul31(v, r);
  output[i] = v;
}
template<int Rounds>
__global__ void m61Chain(GF61 const* input, GF61 const* roots, GF61* output, u32 count) {
  u32 const i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= count) return;
  GF61 v = input[i], r = roots[i];
#pragma unroll
  for (int round = 0; round != Rounds; ++round) v = cmul61(v, r);
  output[i] = v;
}
template<int Rounds>
__global__ void q7Chain(GQ7 const* input, GQ7 const* roots, GQ7* output, u32 count) {
  u32 const i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= count) return;
  GQ7 v = input[i], r = roots[i];
#pragma unroll
  for (int round = 0; round != Rounds; ++round) v = cmulQ7(v, r);
  output[i] = v;
}

// ---------------- radix-7 edge kernels ----------------
// One thread owns a 7-value group (stride EDGE_GROUPS): direct DFT-7 with
// w7 powers, pointwise complex square, inverse DFT-7, scale 1/7.
__global__ void q7Edge(GQ7 const* input, GQ7* output, GQ7 const* w7, GQ7 const* w7i, u32 inv7) {
  u32 const g = blockIdx.x * blockDim.x + threadIdx.x;
  if (g >= EDGE_GROUPS) return;
  GQ7 v[7], t[7];
#pragma unroll
  for (int j = 0; j != 7; ++j) v[j] = input[g + j * EDGE_GROUPS];
#pragma unroll
  for (int k = 0; k != 7; ++k) {
    GQ7 acc = v[0];
#pragma unroll
    for (int n = 1; n != 7; ++n) {
      GQ7 const term = cmulQ7(v[n], w7[(n * k) % 7]);
      acc = {addQ7(acc.x, term.x), addQ7(acc.y, term.y)};
    }
    t[k] = cmulQ7(acc, acc);
  }
#pragma unroll
  for (int n = 0; n != 7; ++n) {
    GQ7 acc = t[0];
#pragma unroll
    for (int k = 1; k != 7; ++k) {
      GQ7 const term = cmulQ7(t[k], w7i[(n * k) % 7]);
      acc = {addQ7(acc.x, term.x), addQ7(acc.y, term.y)};
    }
    output[g + n * EDGE_GROUPS] = cscaleQ7(acc, inv7);
  }
}
__global__ void q7SquareOnly(GQ7 const* input, GQ7* output) {
  u32 const g = blockIdx.x * blockDim.x + threadIdx.x;
  if (g >= EDGE_GROUPS) return;
#pragma unroll
  for (int j = 0; j != 7; ++j) {
    GQ7 const v = input[g + j * EDGE_GROUPS];
    output[g + j * EDGE_GROUPS] = cmulQ7(v, v);
  }
}
__global__ void m61Edge(GF61 const* input, GF61* output, GF61 const* w7, GF61 const* w7i, u64 inv7) {
  u32 const g = blockIdx.x * blockDim.x + threadIdx.x;
  if (g >= EDGE_GROUPS) return;
  GF61 v[7], t[7];
#pragma unroll
  for (int j = 0; j != 7; ++j) v[j] = input[g + j * EDGE_GROUPS];
#pragma unroll
  for (int k = 0; k != 7; ++k) {
    GF61 acc = v[0];
#pragma unroll
    for (int n = 1; n != 7; ++n) {
      GF61 const term = cmul61(v[n], w7[(n * k) % 7]);
      acc = {add61(acc.x, term.x), add61(acc.y, term.y)};
    }
    t[k] = cmul61(acc, acc);
  }
#pragma unroll
  for (int n = 0; n != 7; ++n) {
    GF61 acc = t[0];
#pragma unroll
    for (int k = 1; k != 7; ++k) {
      GF61 const term = cmul61(t[k], w7i[(n * k) % 7]);
      acc = {add61(acc.x, term.x), add61(acc.y, term.y)};
    }
    output[g + n * EDGE_GROUPS] = cscale61(acc, inv7);
  }
}
__global__ void m61SquareOnly(GF61 const* input, GF61* output) {
  u32 const g = blockIdx.x * blockDim.x + threadIdx.x;
  if (g >= EDGE_GROUPS) return;
#pragma unroll
  for (int j = 0; j != 7; ++j) {
    GF61 const v = input[g + j * EDGE_GROUPS];
    output[g + j * EDGE_GROUPS] = cmul61(v, v);
  }
}

// ---------------- host reference ----------------
u64 hMulMod(u64 a, u64 b, u64 modulus) { return u64(u128(a) * b % modulus); }
struct HC { u64 x, y; };
HC hCmul(HC a, HC b, u64 modulus) {
  u64 const ac = hMulMod(a.x, b.x, modulus), bd = hMulMod(a.y, b.y, modulus);
  u64 const real = ac >= bd ? ac - bd : modulus - (bd - ac);
  u64 const imag = (hMulMod(a.x, b.y, modulus) + hMulMod(a.y, b.x, modulus)) % modulus;
  return {real, imag};
}
u64 hPow(u64 base, u64 e, u64 modulus) {
  u64 r = 1;
  while (e) {
    if (e & 1) r = hMulMod(r, base, modulus);
    base = hMulMod(base, base, modulus);
    e >>= 1;
  }
  return r;
}

template<class T>
T* deviceCopy(std::vector<T> const& source) {
  T* result = nullptr;
  CUDA_CHECK(cudaMalloc(&result, source.size() * sizeof(T)));
  CUDA_CHECK(cudaMemcpy(result, source.data(), source.size() * sizeof(T),
                        cudaMemcpyHostToDevice));
  return result;
}

struct Arm {
  std::string name;
  std::function<void(cudaStream_t)> launchA, launchB;
  std::array<float, SAMPLES> samples{};
};

}  // namespace

int main() {
  std::mt19937_64 random(0x7'2026'08'12ull);

  // 7th roots: base-field scalars of order 7 in each field.
  u64 const w61 = [&] {
    for (;;) {
      u64 const h = 2 + random() % (M61 - 3);
      u64 const w = hPow(h, (M61 - 1) / 7, M61);
      if (w != 1) return w;
    }
  }();
  u32 const w7q = [&] {
    for (;;) {
      u64 const h = 2 + random() % (Q7 - 3);
      u64 const w = hPow(h, (Q7 - 1) / 7, Q7);
      if (w != 1) return u32(w);
    }
  }();
  u64 const inv7_61 = hPow(7, M61 - 2, M61);
  u32 const inv7_q = u32(hPow(7, Q7 - 2, Q7));
  if (hPow(w61, 7, M61) != 1 || hPow(w7q, 7, Q7) != 1) { std::fprintf(stderr, "root FAILED\n"); return 3; }

  std::vector<GF61> w61Tab(7), w61iTab(7);
  std::vector<GQ7> wqTab(7), wqiTab(7);
  for (u32 j = 0; j != 7; ++j) {
    u64 const p61 = hPow(w61, j, M61);
    u64 const pq = hPow(w7q, j, Q7);
    w61Tab[j] = {p61, 0};
    w61iTab[j] = {hPow(p61, M61 - 2, M61), 0};
    // Montgomery form for q7 (value * R)
    wqTab[j] = {u32(hMulMod(pq, Q7_R, Q7)), 0};
    wqiTab[j] = {u32(hMulMod(hPow(pq, Q7 - 2, Q7), Q7_R, Q7)), 0};
  }
  u32 const inv7_qMont = u32(hMulMod(inv7_q, Q7_R, Q7));

  // ---- data ----
  std::vector<GF31> c31(CONTROL_COUNT), c31R(CONTROL_COUNT);
  std::vector<GF61> c61(CONTROL_COUNT), c61R(CONTROL_COUNT);
  std::vector<GF61> k61(CAND_COUNT), k61R(CAND_COUNT);
  std::vector<GQ7> kq(CAND_COUNT), kqR(CAND_COUNT);
  for (u32 i = 0; i != CONTROL_COUNT; ++i) {
    c31[i] = {u32(random() % M31), u32(random() % M31)};
    c31R[i] = {u32(random() % M31), u32(random() % M31)};
    c61[i] = {random() % M61, random() % M61};
    c61R[i] = {random() % M61, random() % M61};
  }
  for (u32 i = 0; i != CAND_COUNT; ++i) {
    k61[i] = {random() % M61, random() % M61};
    k61R[i] = {random() % M61, random() % M61};
    kq[i] = {u32(random() % Q7), u32(random() % Q7)};    // Montgomery-form residues
    kqR[i] = {u32(random() % Q7), u32(random() % Q7)};
  }

  GF31 *c31D = deviceCopy(c31), *c31RD = deviceCopy(c31R), *c31O{};
  GF61 *c61D = deviceCopy(c61), *c61RD = deviceCopy(c61R), *c61O{};
  GF61 *k61D = deviceCopy(k61), *k61RD = deviceCopy(k61R), *k61O{};
  GQ7 *kqD = deviceCopy(kq), *kqRD = deviceCopy(kqR), *kqO{};
  GF61 *w61D = deviceCopy(w61Tab), *w61iD = deviceCopy(w61iTab);
  GQ7 *wqD = deviceCopy(wqTab), *wqiD = deviceCopy(wqiTab);
  CUDA_CHECK(cudaMalloc(&c31O, CONTROL_COUNT * sizeof(GF31)));
  CUDA_CHECK(cudaMalloc(&c61O, CONTROL_COUNT * sizeof(GF61)));
  CUDA_CHECK(cudaMalloc(&k61O, CAND_COUNT * sizeof(GF61)));
  CUDA_CHECK(cudaMalloc(&kqO, CAND_COUNT * sizeof(GQ7)));

  cudaStream_t streamA{}, streamB{}, control{};
  CUDA_CHECK(cudaStreamCreate(&streamA));
  CUDA_CHECK(cudaStreamCreate(&streamB));
  CUDA_CHECK(cudaStreamCreate(&control));

  dim3 const block(256);
  dim3 const gridCtl((CONTROL_COUNT + 255) / 256), gridCand((CAND_COUNT + 255) / 256),
       gridEdge((EDGE_GROUPS + 255) / 256);

  auto l31c7 = [&](cudaStream_t s) { m31Chain<7><<<gridCtl, block, 0, s>>>(c31D, c31RD, c31O, CONTROL_COUNT); };
  auto l61c7 = [&](cudaStream_t s) { m61Chain<7><<<gridCtl, block, 0, s>>>(c61D, c61RD, c61O, CONTROL_COUNT); };
  auto lq6 = [&](cudaStream_t s) { q7Chain<6><<<gridCand, block, 0, s>>>(kqD, kqRD, kqO, CAND_COUNT); };
  auto l61c6 = [&](cudaStream_t s) { m61Chain<6><<<gridCand, block, 0, s>>>(k61D, k61RD, k61O, CAND_COUNT); };
  auto lqe = [&](cudaStream_t s) { q7Edge<<<gridEdge, block, 0, s>>>(kqD, kqO, wqD, wqiD, inv7_qMont); };
  auto lqs = [&](cudaStream_t s) { q7SquareOnly<<<gridEdge, block, 0, s>>>(kqD, kqO); };
  auto l61e = [&](cudaStream_t s) { m61Edge<<<gridEdge, block, 0, s>>>(k61D, k61O, w61D, w61iD, inv7_61); };
  auto l61s = [&](cudaStream_t s) { m61SquareOnly<<<gridEdge, block, 0, s>>>(k61D, k61O); };

  // ---- exactness ----
  {
    // q7 chain: full-array host oracle (Montgomery: each product removes one R)
    lq6(streamA);
    CUDA_CHECK(cudaStreamSynchronize(streamA));
    std::vector<GQ7> got(CAND_COUNT);
    CUDA_CHECK(cudaMemcpy(got.data(), kqO, CAND_COUNT * sizeof(GQ7), cudaMemcpyDeviceToHost));
    for (u32 i = 0; i != CAND_COUNT; ++i) {
      HC v{kq[i].x, kq[i].y}, r{kqR[i].x, kqR[i].y};
      for (int round = 0; round != 6; ++round) {
        v = hCmul(v, r, Q7);
        v.x = hMulMod(v.x, Q7_R_INV, Q7);
        v.y = hMulMod(v.y, Q7_R_INV, Q7);
      }
      if (got[i].x != v.x || got[i].y != v.y) { std::fprintf(stderr, "q7 chain oracle FAILED at %u\n", i); return 2; }
    }
    // q7 edge: cyclic-7 self-convolution oracle, full array.  Device values
    // are Montgomery vR; DFT+square+inverse gives conv results scaled by R^-2
    // per product pair... verify directly against a plain-domain convolution
    // of the plain values, tracking the R factors explicitly.
    lqe(streamA);
    CUDA_CHECK(cudaStreamSynchronize(streamA));
    CUDA_CHECK(cudaMemcpy(got.data(), kqO, CAND_COUNT * sizeof(GQ7), cudaMemcpyDeviceToHost));
    for (u32 g = 0; g != EDGE_GROUPS; ++g) {
      HC plain[7], conv[7];
      for (int j = 0; j != 7; ++j) {   // strip one R: plain = value * R^-1... keep in R-domain and track
        plain[j] = {kq[g + j * EDGE_GROUPS].x, kq[g + j * EDGE_GROUPS].y};
      }
      for (int j = 0; j != 7; ++j) conv[j] = {0, 0};
      for (int n = 0; n != 7; ++n)
        for (int m = 0; m != 7; ++m) {
          int const k = (n + m) % 7;
          HC const term = hCmul(plain[n], plain[m], Q7);   // = pn*pm*R^-1 in Mont terms
          conv[k].x = (conv[k].x + term.x) % Q7;
          conv[k].y = (conv[k].y + term.y) % Q7;
        }
      // Degree bookkeeping: stored values are xR.  Device: DFT keeps R-form
      // (X R); square gives X^2 R; inverse keeps R; scale by inv7*R keeps R
      // -> output = conv_x * R.  Host conv[] above = Sum (x_n R)(x_m R) mod q
      // = conv_x * R^2.  So expected output = conv[] * R^-1.
      for (int j = 0; j != 7; ++j) {
        u64 const ex = hMulMod(conv[j].x, Q7_R_INV, Q7), ey = hMulMod(conv[j].y, Q7_R_INV, Q7);
        if (got[g + j * EDGE_GROUPS].x != ex || got[g + j * EDGE_GROUPS].y != ey) {
          std::fprintf(stderr, "q7 edge oracle FAILED group %u slot %d\n", g, j);
          return 2;
        }
      }
      if (g == 4096) break;   // 4096 groups is ample and keeps host time low
    }
    // m61 edge oracle
    l61e(streamA);
    CUDA_CHECK(cudaStreamSynchronize(streamA));
    std::vector<GF61> got61(CAND_COUNT);
    CUDA_CHECK(cudaMemcpy(got61.data(), k61O, CAND_COUNT * sizeof(GF61), cudaMemcpyDeviceToHost));
    for (u32 g = 0; g != 4096; ++g) {
      HC plain[7], conv[7];
      for (int j = 0; j != 7; ++j) plain[j] = {k61[g + j * EDGE_GROUPS].x, k61[g + j * EDGE_GROUPS].y};
      for (int j = 0; j != 7; ++j) conv[j] = {0, 0};
      for (int n = 0; n != 7; ++n)
        for (int m = 0; m != 7; ++m) {
          int const k = (n + m) % 7;
          HC const term = hCmul(plain[n], plain[m], M61);
          conv[k].x = (conv[k].x + term.x) % M61;
          conv[k].y = (conv[k].y + term.y) % M61;
        }
      for (int j = 0; j != 7; ++j) {
        if (got61[g + j * EDGE_GROUPS].x != conv[j].x || got61[g + j * EDGE_GROUPS].y != conv[j].y) {
          std::fprintf(stderr, "m61 edge oracle FAILED group %u slot %d\n", g, j);
          return 2;
        }
      }
    }
    std::printf("q7 chain (full) + q7 edge (4096 groups) + m61 edge (4096 groups) oracles: PASS\n");
  }

  // ---- measurement ----
  std::vector<Arm> arms;
  arms.push_back({"core ctl (M31c7|M61c7)", l31c7, l61c7});
  arms.push_back({"core cand (q7c6|M61c6)", lq6, l61c6});
  arms.push_back({"edge seq (q7e;M61e)", [&](cudaStream_t s) { lqe(s); l61e(s); }, nullptr});
  arms.push_back({"sq seq (q7s;M61s)", [&](cudaStream_t s) { lqs(s); l61s(s); }, nullptr});
  arms.push_back({"edge conc (q7e|M61e)", lqe, l61e});
  arms.push_back({"sq conc (q7s|M61s)", lqs, l61s});

  cudaEvent_t start{}, doneA{}, doneB{}, stop{};
  CUDA_CHECK(cudaEventCreate(&start));
  CUDA_CHECK(cudaEventCreate(&doneA));
  CUDA_CHECK(cudaEventCreate(&doneB));
  CUDA_CHECK(cudaEventCreate(&stop));
  for (int s = 0; s != SAMPLES; ++s) {
    for (std::size_t k = 0; k != arms.size(); ++k) {
      Arm& arm = arms[s % 2 ? arms.size() - 1 - k : k];
      CUDA_CHECK(cudaEventRecord(start, control));
      CUDA_CHECK(cudaStreamWaitEvent(streamA, start));
      CUDA_CHECK(cudaStreamWaitEvent(streamB, start));
      for (int repeat = 0; repeat != REPEATS; ++repeat) {
        arm.launchA(streamA);
        if (arm.launchB) arm.launchB(streamB);
      }
      CUDA_CHECK(cudaEventRecord(doneA, streamA));
      CUDA_CHECK(cudaEventRecord(doneB, streamB));
      CUDA_CHECK(cudaStreamWaitEvent(control, doneA));
      CUDA_CHECK(cudaStreamWaitEvent(control, doneB));
      CUDA_CHECK(cudaEventRecord(stop, control));
      CUDA_CHECK(cudaEventSynchronize(stop));
      float sample = 0.0f;
      CUDA_CHECK(cudaEventElapsedTime(&sample, start, stop));
      arm.samples[s] = sample / REPEATS;
    }
  }
  CUDA_CHECK(cudaEventDestroy(start));
  CUDA_CHECK(cudaEventDestroy(doneA));
  CUDA_CHECK(cudaEventDestroy(doneB));
  CUDA_CHECK(cudaEventDestroy(stop));

  float med[6] = {};
  std::printf("medians (%d samples x %d repeats):\n", SAMPLES, REPEATS);
  for (std::size_t k = 0; k != arms.size(); ++k) {
    auto sorted = arms[k].samples;
    std::sort(sorted.begin(), sorted.end());
    med[k] = sorted[SAMPLES / 2];
    std::printf("  %-24s %.3f us\n", arms[k].name.c_str(), 1000.0f * med[k]);
  }
  float const S = med[0] - med[1];
  float const Eseq = med[2] - med[3];
  float const Econc = med[4] - med[5];
  std::printf("core saving S = %+.3f us (Sol 300W: 18.67-20.30)\n"
              "edge increment sequential E_seq = %+.3f us (Sol Rader, 300W: 20.13-20.29)\n"
              "edge increment co-run    E_conc = %+.3f us\n"
              "verdict margins: S - E_seq = %+.3f us | S - E_conc = %+.3f us\n",
              1000.0f * S, 1000.0f * Eseq, 1000.0f * Econc,
              1000.0f * (S - Eseq), 1000.0f * (S - Econc));
  return 0;
}
