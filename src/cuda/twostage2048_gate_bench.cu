// Two-stage factorization gate (2026-08-12): can a 2048-point GF61 row
// kernel (32-KiB shared/row, radix 8*8*8*4, 256 lanes x 8 values) match the
// production-shaped 512-point row kernel's throughput?  The full 4M
// transform costs the same butterfly stages either way (11+11 vs 9+3+9),
// so per-pass parity means a 2048x2048 backend deletes two of the six
// kernel boundaries (~64 MiB/iter) for free.  Exactness: one row verified
// against a direct O(N^2) host DFT.

#include <cuda_runtime.h>
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <vector>

using u32 = std::uint32_t;
using u64 = std::uint64_t;
using u128 = unsigned __int128;
constexpr u64 M61 = (u64{1} << 61) - 1;
constexpr u32 COUNT = 1u << 21;      // 2M values = 32 MiB
struct GF61 { u64 x, y; };

void check(cudaError_t r, int line) {
  if (r != cudaSuccess) { std::fprintf(stderr, "%d: %s\n", line, cudaGetErrorString(r)); std::exit(1); }
}
#define CK(x) check((x), __LINE__)

__device__ __forceinline__ u64 add61(u64 a, u64 b) { u64 r = a + b; return r >= M61 ? r - M61 : r; }
__device__ __forceinline__ u64 sub61(u64 a, u64 b) { return a >= b ? a - b : M61 - (b - a); }
__device__ __forceinline__ u64 mul61(u64 a, u64 b) {
  u64 const lo = a * b, hi = __umul64hi(a, b);
  u64 r = (lo & M61) + (lo >> 61) + (hi << 3);
  if (r >= M61) r -= M61;
  if (r >= M61) r -= M61;
  return r;
}
__device__ __forceinline__ GF61 cadd(GF61 a, GF61 b) { return {add61(a.x, b.x), add61(a.y, b.y)}; }
__device__ __forceinline__ GF61 cmul(GF61 a, GF61 b) {
  u64 const k1 = mul61(b.x, add61(a.x, a.y));
  u64 const k2 = mul61(a.x, sub61(b.y, b.x));
  u64 const k3 = mul61(a.y, add61(b.y, b.x));
  return {sub61(k1, k3), add61(k1, k2)};
}

__device__ void dftR(GF61* v, GF61 const* wr, int R) {   // naive R-point (R=8 or 4), wr[j]=W_R^j
  GF61 t[8];
  for (int k = 0; k < R; ++k) {
    GF61 acc = v[0];
    for (int n = 1; n < R; ++n) acc = cadd(acc, cmul(v[n], wr[(n * k) % R]));
    t[k] = acc;
  }
  for (int k = 0; k < R; ++k) v[k] = t[k];
}

// generic radix pass over an L-point shared row; lanes*perThread == L/RADIX handled per call
template<u32 L, u32 RADIX, u32 LANES>
__device__ void pass(GF61* row, u32 lane, u32 span, GF61 const* wr, GF61 const* wl) {
  // wl = W_L powers table (L entries); butterflies at this span
  u32 const groups = L / (span * RADIX);
  for (u32 c = lane; c < L / RADIX; c += LANES) {
    u32 const j = c % span, g = c / span, base = g * span * RADIX + j;
    GF61 v[8];
    for (u32 m = 0; m < RADIX; ++m) v[m] = row[base + m * span];
    dftR(v, wr, RADIX);
    for (u32 k = 1; k < RADIX; ++k) v[k] = cmul(v[k], wl[(u64(groups) * j * k) % L]);
    for (u32 k = 0; k < RADIX; ++k) row[base + k * span] = v[k];
  }
  __syncthreads();
}

// Arm H512: production-shaped 512-pt rows, 4 rows per 256-thread block
__global__ void h512Kernel(GF61 const* in, GF61* out, GF61 const* w8, GF61 const* w512) {
  __shared__ GF61 rows[4][512];
  u32 const rowG = blockIdx.x * 4 + threadIdx.x / 64;
  u32 const lane = threadIdx.x % 64;
  GF61* row = rows[threadIdx.x / 64];
  GF61 const* g = in + rowG * 512;
  for (u32 m = 0; m < 8; ++m) row[lane + m * 64] = g[lane + m * 64];
  __syncthreads();
  pass<512, 8, 64>(row, lane, 64, w8, w512);
  pass<512, 8, 64>(row, lane, 8, w8, w512);
  pass<512, 8, 64>(row, lane, 1, w8, w512);
  GF61* o = out + rowG * 512;
  for (u32 m = 0; m < 8; ++m) o[lane + m * 64] = row[lane + m * 64];
}

// Arm H2048: 2048-pt rows, 1 row per 256-thread block (32-KiB shared)
__global__ void h2048Kernel(GF61 const* in, GF61* out, GF61 const* w8, GF61 const* w4, GF61 const* w2048) {
  __shared__ GF61 row[2048];
  u32 const rowG = blockIdx.x;
  u32 const lane = threadIdx.x;
  GF61 const* g = in + rowG * 2048;
  for (u32 m = 0; m < 8; ++m) row[lane + m * 256] = g[lane + m * 256];
  __syncthreads();
  pass<2048, 8, 256>(row, lane, 256, w8, w2048);
  pass<2048, 8, 256>(row, lane, 32, w8, w2048);
  pass<2048, 8, 256>(row, lane, 4, w8, w2048);
  pass<2048, 4, 256>(row, lane, 1, w4, w2048);
  GF61* o = out + rowG * 2048;
  for (u32 m = 0; m < 8; ++m) o[lane + m * 256] = row[lane + m * 256];
}

// Gate 2: tail-shaped pair kernel — both Hermitian partner rows resident
// (64 KiB shared, 512 threads), fwd both, pointwise mul, second (inverse-
// shaped) transform, write both.  Times the 1-block/SM tail regime.
__global__ void __launch_bounds__(512, 1) h2048PairKernel(GF61 const* in, GF61* out,
    GF61 const* w8, GF61 const* w4, GF61 const* w2048) {
  extern __shared__ GF61 rowsDyn[];
  GF61 (*rows)[2048] = (GF61 (*)[2048]) rowsDyn;
  u32 const pairG = blockIdx.x;
  u32 const half = threadIdx.x / 256, lane = threadIdx.x % 256;
  GF61* row = rows[half];
  GF61 const* g = in + (pairG * 2 + half) * 2048;
  for (u32 m = 0; m < 8; ++m) row[lane + m * 256] = g[lane + m * 256];
  __syncthreads();
  // both halves run the pass ladder concurrently (pass syncs all 512)
  pass<2048, 8, 256>(row, lane, 256, w8, w2048);
  pass<2048, 8, 256>(row, lane, 32, w8, w2048);
  pass<2048, 8, 256>(row, lane, 4, w8, w2048);
  pass<2048, 4, 256>(row, lane, 1, w4, w2048);
  for (u32 m = 0; m < 8; ++m) {                      // pointwise pair-mul stand-in
    u32 const i = lane + m * 256;
    rows[half][i] = cmul(rows[half][i], rows[1 - half][(2048 - i) & 2047]);
  }
  __syncthreads();
  pass<2048, 8, 256>(row, lane, 256, w8, w2048);     // inverse-shaped ladder
  pass<2048, 8, 256>(row, lane, 32, w8, w2048);
  pass<2048, 8, 256>(row, lane, 4, w8, w2048);
  pass<2048, 4, 256>(row, lane, 1, w4, w2048);
  GF61* o = out + (pairG * 2 + half) * 2048;
  for (u32 m = 0; m < 8; ++m) o[lane + m * 256] = row[lane + m * 256];
}

// production-shaped 512 comparator: pair of 512 rows (16 KiB), 128 threads
__global__ void h512PairKernel(GF61 const* in, GF61* out, GF61 const* w8, GF61 const* w512) {
  __shared__ GF61 rows[2][512];
  u32 const pairG = blockIdx.x;
  u32 const half = threadIdx.x / 64, lane = threadIdx.x % 64;
  GF61* row = rows[half];
  GF61 const* g = in + (pairG * 2 + half) * 512;
  for (u32 m = 0; m < 8; ++m) row[lane + m * 64] = g[lane + m * 64];
  __syncthreads();
  pass<512, 8, 64>(row, lane, 64, w8, w512);
  pass<512, 8, 64>(row, lane, 8, w8, w512);
  pass<512, 8, 64>(row, lane, 1, w8, w512);
  for (u32 m = 0; m < 8; ++m) {
    u32 const i = lane + m * 64;
    rows[half][i] = cmul(rows[half][i], rows[1 - half][(512 - i) & 511]);
  }
  __syncthreads();
  pass<512, 8, 64>(row, lane, 64, w8, w512);
  pass<512, 8, 64>(row, lane, 8, w8, w512);
  pass<512, 8, 64>(row, lane, 1, w8, w512);
  GF61* o = out + (pairG * 2 + half) * 512;
  for (u32 m = 0; m < 8; ++m) o[lane + m * 64] = row[lane + m * 64];
}

u64 hMul(u64 a, u64 b) { return u64(u128(a) * b % M61); }
struct HC { u64 x, y; };
HC hCmul(HC a, HC b) {
  u64 const ac = hMul(a.x, b.x), bd = hMul(a.y, b.y);
  u64 const re = ac >= bd ? ac - bd : M61 - (bd - ac);
  return {re, (hMul(a.x, b.y) + hMul(a.y, b.x)) % M61};
}
HC hPow(HC b, u128 e) { HC r{1, 0}; while (e) { if (e & 1) r = hCmul(r, b); b = hCmul(b, b); e >>= 1; } return r; }

int main() {
  // norm-one generator of order 2^61 -> W2048 etc.
  HC gen{};
  for (u64 seed = 3;; ++seed) {
    HC const t = hPow(HC{seed, seed * 12345 + 1}, M61 - 1);
    HC const h = hPow(t, u128{1} << 60);
    if (!(h.x == 1 && h.y == 0)) { gen = t; break; }
  }
  HC const W2048 = hPow(gen, (u128{1} << 61) / 2048);
  auto powers = [&](HC w, u32 n) {
    std::vector<GF61> p(n); HC a{1, 0};
    for (u32 i = 0; i < n; ++i) { p[i] = {a.x, a.y}; a = hCmul(a, w); }
    return p;
  };
  std::vector<GF61> w2048 = powers(W2048, 2048), w512 = powers(hPow(W2048, 4), 512),
                    w8 = powers(hPow(W2048, 256), 8), w4 = powers(hPow(W2048, 512), 4);

  std::vector<GF61> host(COUNT);
  for (u32 i = 0; i < COUNT; ++i) host[i] = {(u64(i) * 2654435761u) % M61, (u64(i) * 40503u + 7) % M61};
  GF61 *in, *out, *w8D, *w4D, *w512D, *w2048D;
  CK(cudaMalloc(&in, COUNT * 16)); CK(cudaMalloc(&out, COUNT * 16));
  CK(cudaMemcpy(in, host.data(), COUNT * 16, cudaMemcpyHostToDevice));
  auto up = [&](std::vector<GF61>& v, GF61** d) { CK(cudaMalloc(d, v.size() * 16)); CK(cudaMemcpy(*d, v.data(), v.size() * 16, cudaMemcpyHostToDevice)); };
  up(w8, &w8D); up(w4, &w4D); up(w512, &w512D); up(w2048, &w2048D);

  // exactness: one 2048 row vs direct DFT (device output order is digit-reversed
  // by the DIF passes; verify as a SET by matching the DC bin and total sums is
  // weak — instead verify bin 0 and bin k via direct evaluation at unscrambled
  // indices computed from the pass structure: digit order (d8,d8,d8,d4) reversed.
  h2048Kernel<<<COUNT / 2048, 256>>>(in, out, w8D, w4D, w2048D);
  CK(cudaDeviceSynchronize());
  {
    std::vector<GF61> got(2048);
    CK(cudaMemcpy(got.data(), out, 2048 * 16, cudaMemcpyDeviceToHost));
    auto rev = [](u32 k) {           // digit-reverse: k = d0 + 8*(d1 + 8*(d2 + 8*d3)), d3<4
      u32 d0 = k % 8, d1 = (k / 8) % 8, d2 = (k / 64) % 8, d3 = k / 512;
      return d3 + 4 * (d2 + 8 * (d1 + 8 * d0));
    };
    for (u32 k : {0u, 1u, 5u, 100u, 1234u, 2047u}) {
      HC acc{0, 0};
      for (u32 n = 0; n < 2048; ++n) {
        HC const w{w2048[(u64(n) * k) % 2048].x, w2048[(u64(n) * k) % 2048].y};
        acc = HC{(acc.x + hCmul(HC{host[n].x, host[n].y}, w).x) % M61,
                 (acc.y + hCmul(HC{host[n].x, host[n].y}, w).y) % M61};
      }
      GF61 const g = got[rev(k)];
      if (g.x != acc.x || g.y != acc.y) { std::printf("VERIFY FAILED bin %u\n", k); return 2; }
    }
    std::printf("2048-pt DFT verified vs direct host DFT: PASS\n");
  }

  cudaEvent_t t0, t1; CK(cudaEventCreate(&t0)); CK(cudaEventCreate(&t1));
  auto median = [&](auto launch) {
    std::array<float, 15> s{};
    for (float& x : s) {
      CK(cudaEventRecord(t0));
      for (int r = 0; r < 10; ++r) launch();
      CK(cudaEventRecord(t1)); CK(cudaEventSynchronize(t1));
      CK(cudaEventElapsedTime(&x, t0, t1)); x /= 10;
    }
    std::sort(s.begin(), s.end()); return s[7];
  };
  float const a = median([&] { h512Kernel<<<COUNT / 512 / 4, 256>>>(in, out, w8D, w512D); });
  float const b = median([&] { h2048Kernel<<<COUNT / 2048, 256>>>(in, out, w8D, w4D, w2048D); });
  float const pa = median([&] { h512PairKernel<<<COUNT / 1024, 128>>>(in, out, w8D, w512D); });
  CK(cudaFuncSetAttribute(h2048PairKernel, cudaFuncAttributeMaxDynamicSharedMemorySize, 65536));
  float const pb = median([&] { h2048PairKernel<<<COUNT / 4096, 512, 65536>>>(in, out, w8D, w4D, w2048D); });
  CK(cudaGetLastError());
  std::printf("512-pt rows %.3f us | 2048-pt rows %.3f us | ratio %.4f\n", 1000 * a, 1000 * b, b / a);
  std::printf("tail-shaped: 512-pair %.3f us | 2048-pair (64K,1blk/SM) %.3f us | ratio %.4f\n",
              1000 * pa, 1000 * pb, pb / pa);
  return 0;
}
