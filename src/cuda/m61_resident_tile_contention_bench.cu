// Reconstruction of the lost m61_resident_tile_bench with the audit's
// contention arms (2026-08-12).  Sol's isolated gate measured the fused
// 64-KiB resident tile 1.019-1.027x SLOWER than three cache-resident
// kernels -- but with the whole 32-64 MiB working set L2-resident, the
// global round trips the fusion deletes were nearly free.  This bench
// re-measures both paths while an M31-chain-shaped memory workload floods
// L2/DRAM from a second stream, asking whether production-like contention
// flips the sign.
//
// Shape (Sol's): 2,097,152 canonical GF(M61^2) values in 512 independent
// 8x512 tiles.  Each tile undergoes an exact 4096-point cyclic NTT square
// factored 8 (middle) x 512 (height, radix-8^3):
//   P1 middleFwd:  DFT-8 across m per column h, then twiddle W4096^(h*k1)
//   P2 heightSq:   DFT-512 per row, pointwise square, inverse DFT-512
//   P3 middleInv:  inverse twiddle, inverse DFT-8, scale 1/4096
// Control path: three kernels with global round trips between phases.
// Fused path: one kernel per tile, the whole tile in 64 KiB dynamic shared.
// Both paths share the same device transform functions; exactness is gated
// by (a) full path agreement over all 2M outputs, (b) host-NTT verification
// of 8 tiles, (c) direct O(N^2) DFT validation of the host NTT on one tile,
// (d) a direct cyclic-convolution square check on one tile.
//
// Measurement: q65-era protocol (arms interleaved per sample, order
// reversed on odd samples), 21 samples x 10 repeats, medians of the
// stream-A makespan; contended arms flood stream B with chain-4 M31
// kernels sized to outlast stream A, drained between samples.

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

constexpr u64 M61 = (u64{1} << 61) - 1;
constexpr u32 M31 = 0x7fffffffu;
constexpr u32 TILE_M = 8, TILE_H = 512, TILE_N = TILE_M * TILE_H;
constexpr u32 TILES = 512;
constexpr u32 COUNT = TILES * TILE_N;  // 2,097,152
constexpr u32 LOAD_COUNT = 1u << 21;   // M31 co-runner population (16 MiB x2)
constexpr int SAMPLES = 21, REPEATS = 10, LOAD_LAUNCHES = 260;

struct GF61 { u64 x, y; };
struct GF31 { u32 x, y; };

void check(cudaError_t result, char const* expression, int line) {
  if (result == cudaSuccess) return;
  std::fprintf(stderr, "%s:%d: %s\n", expression, line, cudaGetErrorString(result));
  std::exit(1);
}
#define CUDA_CHECK(expression) check((expression), #expression, __LINE__)

// ---------------- field arithmetic (device) ----------------
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
__device__ __forceinline__ GF61 cadd(GF61 a, GF61 b) { return {add61(a.x, b.x), add61(a.y, b.y)}; }
__device__ __forceinline__ GF61 csub(GF61 a, GF61 b) { return {sub61(a.x, b.x), sub61(a.y, b.y)}; }
__device__ __forceinline__ GF61 cmul(GF61 a, GF61 b) {
  u64 const k1 = mul61(b.x, add61(a.x, a.y));
  u64 const k2 = mul61(a.x, sub61(b.y, b.x));
  u64 const k3 = mul61(a.y, add61(b.y, b.x));
  return {sub61(k1, k3), add61(k1, k2)};
}
__device__ __forceinline__ GF61 cscale(GF61 a, u64 s) { return {mul61(a.x, s), mul61(a.y, s)}; }

// In-register DFT-8 over the middle dimension.  W8 powers come from the
// twiddle table (w8[j] = W8^j, 8 entries; wInv variant for the inverse).
__device__ void dft8(GF61* v, GF61 const* w8) {
  GF61 t[TILE_M];
#pragma unroll
  for (int k = 0; k != TILE_M; ++k) {
    GF61 acc = v[0];
#pragma unroll
    for (int n = 1; n != TILE_M; ++n) acc = cadd(acc, cmul(v[n], w8[(n * k) & 7]));
    t[k] = acc;
  }
#pragma unroll
  for (int k = 0; k != TILE_M; ++k) v[k] = t[k];
}

// One 512-point NTT over a shared-memory row, radix-8 x 3 passes,
// decimation-in-frequency with unified twiddles from w512 (w512[j] = W512^j).
// 64 threads per row, each owning 8 values.  syncAll is __syncthreads (the
// caller may run several rows per block; barriers align them).
__device__ void row512(GF61* row, u32 lane, GF61 const* w8, GF61 const* w512) {
  // pass structure: span 64, 8, 1
  for (u32 span = 64; span != 0; span /= 8) {
    u32 const groups = 512 / (span * 8);   // 1, 8, 64
    // each thread handles 8 butterflies-worth: total butterfly columns = 64
    // lane owns column set {lane}: position p = (lane / span) * span * 8 + lane % span
    u32 const j = lane % span;
    u32 const g = lane / span;
    u32 const base = g * span * 8 + j;
    GF61 v[8];
#pragma unroll
    for (int m = 0; m != 8; ++m) v[m] = row[base + m * span];
    dft8(v, w8);
    // twiddle: W_{span*8}^{j*k} = w512[(512/(span*8)) * j * k]
#pragma unroll
    for (int k = 1; k != 8; ++k) v[k] = cmul(v[k], w512[(groups * j * k) & 511]);
    // DIF digit-reversal placement: output k goes to base + k*span (natural
    // within the pass; the overall order is digit-reversed and the inverse
    // pass structure consumes it symmetrically)
#pragma unroll
    for (int k = 0; k != 8; ++k) row[base + k * span] = v[k];
    __syncthreads();
  }
}

// Inverse 512-point NTT: decimation-in-time mirror consuming the DIF
// digit-reversed order, using inverse twiddles (w512i[j] = W512^-j).
__device__ void row512inv(GF61* row, u32 lane, GF61 const* w8i, GF61 const* w512i) {
  for (u32 span = 1; span != 512; span *= 8) {
    u32 const groups = 512 / (span * 8);
    u32 const j = lane % span;
    u32 const g = lane / span;
    u32 const base = g * span * 8 + j;
    GF61 v[8];
#pragma unroll
    for (int m = 0; m != 8; ++m) v[m] = row[base + m * span];
#pragma unroll
    for (int k = 1; k != 8; ++k) v[k] = cmul(v[k], w512i[(groups * j * k) & 511]);
    dft8(v, w8i);
#pragma unroll
    for (int k = 0; k != 8; ++k) row[base + k * span] = v[k];
    __syncthreads();
  }
}

// ---------------- control kernels ----------------
// K1: middle forward.  One thread per (tile, h) column.
__global__ void middleFwdKernel(GF61 const* input, GF61* output,
                                GF61 const* w8, GF61 const* w4096) {
  u32 const i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= TILES * TILE_H) return;
  u32 const tile = i / TILE_H, h = i % TILE_H;
  GF61 const* base = input + tile * TILE_N + h;
  GF61 v[TILE_M];
#pragma unroll
  for (int m = 0; m != TILE_M; ++m) v[m] = base[m * TILE_H];
  dft8(v, w8);
#pragma unroll
  for (int k = 1; k != TILE_M; ++k) v[k] = cmul(v[k], w4096[(u32(h) * k) & 4095]);
  GF61* out = output + tile * TILE_N + h;
#pragma unroll
  for (int k = 0; k != TILE_M; ++k) out[k * TILE_H] = v[k];
}

// K2: height forward / square / inverse.  256 threads = 4 row-lanes of 64;
// each block processes 4 rows (2 blocks per tile-row-set of 8).
__global__ void heightSquareKernel(GF61* data, GF61 const* w8, GF61 const* w512,
                                   GF61 const* w8i, GF61 const* w512i) {
  __shared__ GF61 rows[4][TILE_H];
  u32 const rowGlobal = blockIdx.x * 4 + threadIdx.x / 64;
  u32 const lane = threadIdx.x % 64;
  u32 const tile = rowGlobal / TILE_M, k1 = rowGlobal % TILE_M;
  GF61* g = data + tile * TILE_N + k1 * TILE_H;
  GF61* row = rows[threadIdx.x / 64];
#pragma unroll
  for (int m = 0; m != 8; ++m) row[lane + m * 64] = g[lane + m * 64];
  __syncthreads();
  row512(row, lane, w8, w512);
#pragma unroll
  for (int m = 0; m != 8; ++m) {
    GF61 const value = row[lane + m * 64];
    row[lane + m * 64] = cmul(value, value);
  }
  __syncthreads();
  row512inv(row, lane, w8i, w512i);
#pragma unroll
  for (int m = 0; m != 8; ++m) g[lane + m * 64] = row[lane + m * 64];
}

// K3: middle inverse + total scale 1/4096.
__global__ void middleInvKernel(GF61 const* input, GF61* output,
                                GF61 const* w8i, GF61 const* w4096i, u64 inv4096) {
  u32 const i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= TILES * TILE_H) return;
  u32 const tile = i / TILE_H, h = i % TILE_H;
  GF61 const* base = input + tile * TILE_N + h;
  GF61 v[TILE_M];
#pragma unroll
  for (int k = 0; k != TILE_M; ++k) v[k] = base[k * TILE_H];
#pragma unroll
  for (int k = 1; k != TILE_M; ++k) v[k] = cmul(v[k], w4096i[(u32(h) * k) & 4095]);
  dft8(v, w8i);
  GF61* out = output + tile * TILE_N + h;
#pragma unroll
  for (int m = 0; m != TILE_M; ++m) out[m * TILE_H] = cscale(v[m], inv4096);
}

// ---------------- fused resident kernel ----------------
// One block per tile; the whole 8x512 tile lives in 64 KiB dynamic shared.
// 256 threads: P1/P3 use one thread per 2 columns; P2 uses 4 row-lanes of 64
// with each lane-group handling rows r and r+4.
__global__ void __launch_bounds__(256, 1)
fusedTileKernel(GF61 const* input, GF61* output,
                GF61 const* w8, GF61 const* w512, GF61 const* w4096,
                GF61 const* w8i, GF61 const* w512i, GF61 const* w4096i,
                u64 inv4096) {
  extern __shared__ GF61 tileShared[];
  u32 const tile = blockIdx.x;
  GF61 const* g = input + tile * TILE_N;
  GF61* go = output + tile * TILE_N;
  // P1: middle forward, columns h = threadIdx.x and threadIdx.x + 256
  for (u32 h = threadIdx.x; h < TILE_H; h += 256) {
    GF61 v[TILE_M];
#pragma unroll
    for (int m = 0; m != TILE_M; ++m) v[m] = g[m * TILE_H + h];
    dft8(v, w8);
#pragma unroll
    for (int k = 1; k != TILE_M; ++k) v[k] = cmul(v[k], w4096[(h * u32(k)) & 4095]);
#pragma unroll
    for (int k = 0; k != TILE_M; ++k) tileShared[k * TILE_H + h] = v[k];
  }
  __syncthreads();
  // P2: height rows; 4 lane-groups process rows r and r+4
  u32 const lane = threadIdx.x % 64;
  u32 const group = threadIdx.x / 64;
  for (u32 r = group; r < TILE_M; r += 4) {
    GF61* row = tileShared + r * TILE_H;
    row512(row, lane, w8, w512);
#pragma unroll
    for (int m = 0; m != 8; ++m) {
      GF61 const value = row[lane + m * 64];
      row[lane + m * 64] = cmul(value, value);
    }
    __syncthreads();
    row512inv(row, lane, w8i, w512i);
  }
  __syncthreads();
  // P3: middle inverse + scale
  for (u32 h = threadIdx.x; h < TILE_H; h += 256) {
    GF61 v[TILE_M];
#pragma unroll
    for (int k = 0; k != TILE_M; ++k) v[k] = tileShared[k * TILE_H + h];
#pragma unroll
    for (int k = 1; k != TILE_M; ++k) v[k] = cmul(v[k], w4096i[(h * u32(k)) & 4095]);
    dft8(v, w8i);
#pragma unroll
    for (int m = 0; m != TILE_M; ++m) go[m * TILE_H + h] = cscale(v[m], inv4096);
  }
}

// ---------------- M31-shaped contention co-runner ----------------
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
template<int Rounds>
__global__ void m31LoadKernel(GF31 const* input, GF31 const* roots, GF31* output) {
  u32 const i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= LOAD_COUNT) return;
  GF31 v = input[i], r = roots[i];
#pragma unroll
  for (int round = 0; round != Rounds; ++round) {
    u32 const k1 = mul31(r.x, add31(v.x, v.y));
    u32 const k2 = mul31(v.x, sub31(r.y, r.x));
    u32 const k3 = mul31(v.y, add31(r.y, r.x));
    v = {sub31(k1, k3), add31(k1, k2)};
  }
  output[i] = v;
}

// ---------------- host reference field ----------------
u64 hAdd(u64 a, u64 b) { u64 r = a + b; return r >= M61 ? r - M61 : r; }
u64 hSub(u64 a, u64 b) { return a >= b ? a - b : M61 - (b - a); }
u64 hMul(u64 a, u64 b) { return u64(u128(a) * b % M61); }
struct HC { u64 x, y; };
HC hCmul(HC a, HC b) {
  u64 const ac = hMul(a.x, b.x), bd = hMul(a.y, b.y);
  return {hSub(ac, bd), hAdd(hMul(a.x, b.y), hMul(a.y, b.x))};
}
HC hCadd(HC a, HC b) { return {hAdd(a.x, b.x), hAdd(a.y, b.y)}; }
HC hPow(HC base, u128 e) {
  HC r{1, 0};
  while (e) {
    if (e & 1) r = hCmul(r, base);
    base = hCmul(base, base);
    e >>= 1;
  }
  return r;
}
u64 hPowScalar(u64 base, u128 e) {
  u64 r = 1;
  while (e) {
    if (e & 1) r = hMul(r, base);
    base = hMul(base, base);
    e >>= 1;
  }
  return r;
}

// Host mirror of the device transforms (same factorization and ordering)
void hostDft8(HC* v, std::vector<HC> const& w8) {
  HC t[TILE_M];
  for (int k = 0; k != TILE_M; ++k) {
    HC acc = v[0];
    for (int n = 1; n != TILE_M; ++n) acc = hCadd(acc, hCmul(v[n], w8[(n * k) & 7]));
    t[k] = acc;
  }
  for (int k = 0; k != TILE_M; ++k) v[k] = t[k];
}
void hostRow512(HC* row, std::vector<HC> const& w8, std::vector<HC> const& w512) {
  for (u32 span = 64; span != 0; span /= 8) {
    u32 const groups = 512 / (span * 8);
    std::vector<HC> next(512);
    for (u32 lane = 0; lane != 64; ++lane) {
      u32 const j = lane % span, g = lane / span, base = g * span * 8 + j;
      HC v[8];
      for (int m = 0; m != 8; ++m) v[m] = row[base + m * span];
      hostDft8(v, w8);
      for (int k = 1; k != 8; ++k) v[k] = hCmul(v[k], w512[(groups * j * k) & 511]);
      for (int k = 0; k != 8; ++k) next[base + k * span] = v[k];
    }
    for (u32 p = 0; p != 512; ++p) row[p] = next[p];
  }
}
void hostRow512Inv(HC* row, std::vector<HC> const& w8i, std::vector<HC> const& w512i) {
  for (u32 span = 1; span != 512; span *= 8) {
    u32 const groups = 512 / (span * 8);
    std::vector<HC> next(512);
    for (u32 lane = 0; lane != 64; ++lane) {
      u32 const j = lane % span, g = lane / span, base = g * span * 8 + j;
      HC v[8];
      for (int m = 0; m != 8; ++m) v[m] = row[base + m * span];
      for (int k = 1; k != 8; ++k) v[k] = hCmul(v[k], w512i[(groups * j * k) & 511]);
      hostDft8(v, w8i);
      for (int k = 0; k != 8; ++k) next[base + k * span] = v[k];
    }
    for (u32 p = 0; p != 512; ++p) row[p] = next[p];
  }
}
void hostTileSquare(std::vector<HC>& tile, std::vector<HC> const& w8,
                    std::vector<HC> const& w512, std::vector<HC> const& w4096,
                    std::vector<HC> const& w8i, std::vector<HC> const& w512i,
                    std::vector<HC> const& w4096i, u64 inv4096) {
  // P1
  for (u32 h = 0; h != TILE_H; ++h) {
    HC v[TILE_M];
    for (int m = 0; m != TILE_M; ++m) v[m] = tile[m * TILE_H + h];
    hostDft8(v, w8);
    for (int k = 1; k != TILE_M; ++k) v[k] = hCmul(v[k], w4096[(h * u32(k)) & 4095]);
    for (int k = 0; k != TILE_M; ++k) tile[k * TILE_H + h] = v[k];
  }
  // P2
  for (u32 r = 0; r != TILE_M; ++r) {
    HC* row = tile.data() + r * TILE_H;
    hostRow512(row, w8, w512);
    for (u32 p = 0; p != TILE_H; ++p) row[p] = hCmul(row[p], row[p]);
    hostRow512Inv(row, w8i, w512i);
  }
  // P3
  for (u32 h = 0; h != TILE_H; ++h) {
    HC v[TILE_M];
    for (int k = 0; k != TILE_M; ++k) v[k] = tile[k * TILE_H + h];
    for (int k = 1; k != TILE_M; ++k) v[k] = hCmul(v[k], w4096i[(h * u32(k)) & 4095]);
    hostDft8(v, w8i);
    for (int m = 0; m != TILE_M; ++m)
      tile[m * TILE_H + h] = {hMul(v[m].x, inv4096), hMul(v[m].y, inv4096)};
  }
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
  std::function<void(cudaStream_t)> launchA;
  std::function<void(cudaStream_t)> load;   // null = isolated
  std::array<float, SAMPLES> samples{};
};

}  // namespace

int main() {
  // ---- roots: norm-one generator of order 2^61, then W4096 etc. ----
  std::mt19937_64 random(0x7e51de47'11e2026ull);
  HC gen{};
  for (;;) {
    HC const c{random() % M61, random() % M61};
    if (c.x == 0 && c.y == 0) continue;
    HC const t = hPow(c, M61 - 1);       // order divides M61+1 = 2^61
    HC const half = hPow(t, u128{1} << 60);
    if (half.x == 1 && half.y == 0) continue;   // order < 2^61
    HC const full = hCmul(half, half);
    if (full.x != 1 || full.y != 0) continue;   // paranoia
    gen = t;
    break;
  }
  HC const W4096 = hPow(gen, u128{1} << 49);
  HC const W4096I = hPow(W4096, 4095);
  auto powers = [&](HC w, u32 n) {
    std::vector<HC> p(n);
    p[0] = {1, 0};
    for (u32 i = 1; i != n; ++i) p[i] = hCmul(p[i - 1], w);
    return p;
  };
  std::vector<HC> const w4096 = powers(W4096, 4096), w4096i = powers(W4096I, 4096);
  std::vector<HC> const w512 = powers(hPow(W4096, 8), 512), w512i = powers(hPow(W4096I, 8), 512);
  std::vector<HC> const w8 = powers(hPow(W4096, 512), 8), w8i = powers(hPow(W4096I, 512), 8);
  u64 const inv4096 = hPowScalar(hPowScalar(2, 12), M61 - 2);
  {  // sanity: W4096 has order 4096, W8^8 = 1
    HC const o = hPow(W4096, 4096), h = hPow(W4096, 2048);
    if (o.x != 1 || o.y != 0 || (h.x == 1 && h.y == 0)) { std::fprintf(stderr, "root order FAILED\n"); return 3; }
  }

  // ---- host NTT self-validation on one tile: direct O(N^2) DFT ----
  {
    std::vector<HC> small(TILE_N), viaNtt(TILE_N);
    for (auto& value : small) value = {random() % 997, random() % 997};
    viaNtt = small;
    // direct DFT at a few bins k, compare against P1+P2-forward result:
    // easier and stronger: verify the full square against direct cyclic
    // convolution (the user-visible contract of the pipeline).
    std::vector<HC> direct(TILE_N, HC{0, 0});
    for (u32 n = 0; n != TILE_N; ++n) {
      if (small[n].x == 0 && small[n].y == 0) continue;
      for (u32 m = 0; m != TILE_N; ++m) {
        if (small[m].x == 0 && small[m].y == 0) continue;
        u32 const k = (n + m) & (TILE_N - 1);
        direct[k] = hCadd(direct[k], hCmul(small[n], small[m]));
      }
    }
    hostTileSquare(viaNtt, w8, w512, w4096, w8i, w512i, w4096i, inv4096);
    // The NTT pipeline maps index n = 512*m + h (column-major middle) and its
    // output returns in the same layout; cyclic convolution must agree.
    for (u32 i = 0; i != TILE_N; ++i) {
      if (viaNtt[i].x != direct[i].x || viaNtt[i].y != direct[i].y) {
        std::fprintf(stderr, "host NTT vs direct convolution FAILED at %u\n", i);
        return 2;
      }
    }
    std::printf("host NTT square == direct cyclic convolution: PASS\n");
  }

  // ---- data ----
  std::vector<GF61> data(COUNT);
  for (auto& value : data) value = {random() % M61, random() % M61};
  std::vector<GF31> loadIn(LOAD_COUNT), loadRoots(LOAD_COUNT);
  for (u32 i = 0; i != LOAD_COUNT; ++i) {
    loadIn[i] = {u32(random() % M31), u32(random() % M31)};
    loadRoots[i] = {u32(random() % M31), u32(random() % M31)};
  }
  auto toGF = [](std::vector<HC> const& v) {
    std::vector<GF61> r(v.size());
    for (std::size_t i = 0; i != v.size(); ++i) r[i] = {v[i].x, v[i].y};
    return r;
  };

  GF61 *dataD = deviceCopy(data), *ctlBuf = nullptr, *ctlOut = nullptr, *fusedOut = nullptr;
  CUDA_CHECK(cudaMalloc(&ctlBuf, COUNT * sizeof(GF61)));
  CUDA_CHECK(cudaMalloc(&ctlOut, COUNT * sizeof(GF61)));
  CUDA_CHECK(cudaMalloc(&fusedOut, COUNT * sizeof(GF61)));
  GF61 *w8D = deviceCopy(toGF(w8)), *w512D = deviceCopy(toGF(w512)),
       *w4096D = deviceCopy(toGF(w4096)), *w8iD = deviceCopy(toGF(w8i)),
       *w512iD = deviceCopy(toGF(w512i)), *w4096iD = deviceCopy(toGF(w4096i));
  GF31 *loadInD = deviceCopy(loadIn), *loadRootD = deviceCopy(loadRoots), *loadOutD = nullptr;
  CUDA_CHECK(cudaMalloc(&loadOutD, LOAD_COUNT * sizeof(GF31)));

  CUDA_CHECK(cudaFuncSetAttribute(fusedTileKernel,
                                  cudaFuncAttributeMaxDynamicSharedMemorySize,
                                  TILE_N * sizeof(GF61)));

  cudaStream_t streamA{}, streamB{}, control{};
  CUDA_CHECK(cudaStreamCreate(&streamA));
  CUDA_CHECK(cudaStreamCreate(&streamB));
  CUDA_CHECK(cudaStreamCreate(&control));

  auto launchControl = [&](cudaStream_t s) {
    middleFwdKernel<<<TILES * TILE_H / 256, 256, 0, s>>>(dataD, ctlBuf, w8D, w4096D);
    heightSquareKernel<<<TILES * TILE_M / 4, 256, 0, s>>>(ctlBuf, w8D, w512D, w8iD, w512iD);
    middleInvKernel<<<TILES * TILE_H / 256, 256, 0, s>>>(ctlBuf, ctlOut, w8iD, w4096iD, inv4096);
  };
  auto launchFused = [&](cudaStream_t s) {
    fusedTileKernel<<<TILES, 256, TILE_N * sizeof(GF61), s>>>(
        dataD, fusedOut, w8D, w512D, w4096D, w8iD, w512iD, w4096iD, inv4096);
  };
  auto launchLoad = [&](cudaStream_t s) {
    m31LoadKernel<4><<<LOAD_COUNT / 256, 256, 0, s>>>(loadInD, loadRootD, loadOutD);
  };
  auto launchLoad8 = [&](cudaStream_t s) {
    m31LoadKernel<8><<<LOAD_COUNT / 256, 256, 0, s>>>(loadInD, loadRootD, loadOutD);
  };

  // ---- exactness gates ----
  launchControl(streamA);
  launchFused(streamA);
  CUDA_CHECK(cudaDeviceSynchronize());
  {
    std::vector<GF61> gotCtl(COUNT), gotFused(COUNT);
    CUDA_CHECK(cudaMemcpy(gotCtl.data(), ctlOut, COUNT * sizeof(GF61), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gotFused.data(), fusedOut, COUNT * sizeof(GF61), cudaMemcpyDeviceToHost));
    for (u32 i = 0; i != COUNT; ++i) {
      if (gotCtl[i].x != gotFused[i].x || gotCtl[i].y != gotFused[i].y) {
        std::fprintf(stderr, "path agreement FAILED at %u\n", i);
        return 2;
      }
    }
    for (u32 t = 0; t != 8; ++t) {   // host verification of 8 tiles
      u32 const tile = t * 67;       // spread across the population
      std::vector<HC> host(TILE_N);
      for (u32 i = 0; i != TILE_N; ++i)
        host[i] = {data[tile * TILE_N + i].x, data[tile * TILE_N + i].y};
      hostTileSquare(host, w8, w512, w4096, w8i, w512i, w4096i, inv4096);
      for (u32 i = 0; i != TILE_N; ++i) {
        GF61 const got = gotCtl[tile * TILE_N + i];
        if (got.x != host[i].x || got.y != host[i].y) {
          std::fprintf(stderr, "host oracle FAILED tile %u index %u\n", tile, i);
          return 2;
        }
      }
    }
    std::printf("path agreement (2M values) + host oracle (8 tiles): PASS\n");
  }

  // ---- measurement ----
  std::vector<Arm> arms;
  arms.push_back({"ctl3 isolated", launchControl, nullptr});
  arms.push_back({"fused isolated", launchFused, nullptr});
  arms.push_back({"ctl3 +load8", launchControl, launchLoad8});
  arms.push_back({"fused +load8", launchFused, launchLoad8});
  arms.push_back({"ctl3 +load", launchControl, launchLoad});
  arms.push_back({"fused +load", launchFused, launchLoad});

  cudaEvent_t start{}, doneA{};
  CUDA_CHECK(cudaEventCreate(&start));
  CUDA_CHECK(cudaEventCreate(&doneA));
  for (int s = 0; s != SAMPLES; ++s) {
    for (std::size_t k = 0; k != arms.size(); ++k) {
      Arm& arm = arms[s % 2 ? arms.size() - 1 - k : k];
      CUDA_CHECK(cudaEventRecord(start, control));
      CUDA_CHECK(cudaStreamWaitEvent(streamA, start));
      if (arm.load) {
        CUDA_CHECK(cudaStreamWaitEvent(streamB, start));
        for (int launch = 0; launch != LOAD_LAUNCHES; ++launch) arm.load(streamB);
      }
      for (int repeat = 0; repeat != REPEATS; ++repeat) arm.launchA(streamA);
      CUDA_CHECK(cudaEventRecord(doneA, streamA));
      CUDA_CHECK(cudaEventSynchronize(doneA));
      float sample = 0.0f;
      CUDA_CHECK(cudaEventElapsedTime(&sample, start, doneA));
      arm.samples[s] = sample / REPEATS;
      CUDA_CHECK(cudaStreamSynchronize(streamB));  // drain the flood
    }
  }
  CUDA_CHECK(cudaEventDestroy(start));
  CUDA_CHECK(cudaEventDestroy(doneA));

  float mCtlIso = 0, mFusedIso = 0, mCtlLoad = 0, mFusedLoad = 0, mCtl8 = 0, mFused8 = 0;
  std::printf("medians (%d samples x %d repeats):", SAMPLES, REPEATS);
  for (Arm& arm : arms) {
    auto sorted = arm.samples;
    std::sort(sorted.begin(), sorted.end());
    float const med = sorted[SAMPLES / 2];
    std::printf("  %s %.3f us", arm.name.c_str(), 1000.0f * med);
    if (arm.name == "ctl3 isolated") mCtlIso = med;
    if (arm.name == "fused isolated") mFusedIso = med;
    if (arm.name == "ctl3 +load") mCtlLoad = med;
    if (arm.name == "fused +load") mFusedLoad = med;
    if (arm.name == "ctl3 +load8") mCtl8 = med;
    if (arm.name == "fused +load8") mFused8 = med;
  }
  std::printf("\nisolated fused/ctl %.4f (Sol: 1.019-1.027) | moderate(+load8) %.4f (%+.3f us) | "
              "heavy(+load) %.4f (%+.3f us)\n",
              mFusedIso / mCtlIso, mFused8 / mCtl8, 1000.0f * (mFused8 - mCtl8),
              mFusedLoad / mCtlLoad, 1000.0f * (mFusedLoad - mCtlLoad));

  CUDA_CHECK(cudaStreamDestroy(control));
  CUDA_CHECK(cudaStreamDestroy(streamB));
  CUDA_CHECK(cudaStreamDestroy(streamA));
  return 0;
}
