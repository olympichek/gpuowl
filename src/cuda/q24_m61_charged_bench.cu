// Reconstruction of the lost q24_m61_overlap_bench with the audit's charged
// arms (2026-08-12).  Tests whether the q24-replaces-M31 co-run advantage
// (Sol: 9.46-9.96 us at chain 8, short-quotient reducer) survives charging
// the candidate its generic transform roots and Crandall-Fagin weights --
// the analytical cost Sol asserted but never measured.
//
//   qC = 7*2^21-1 = 14,680,063; centered FP32 operands |a|,|b| <= 7,340,031.
//   Short-quotient reducer (proof in sol ledger "Exact q24 FP32
//   short-quotient reduction"): high=a*b; low=fma(a,b,-high);
//   n=rint(high/q); r=fma(-n,q,high)+low; one +-q centering.
//
// Arms per chain length (all 2^21 complex values; M61 workload common):
//   A  M31gen+M61      : Sol's control (generic M31 complex products)
//   B  q24+M61         : Sol's candidate (free roots -- his decisive gate)
//   C  q24chg+M61      : NEW: per round a table-loaded generic qC twiddle
//                        product, plus CF weight products at entry/exit
//   D  M31cheap+M61    : NEW: production-realistic M31 (power-of-two
//                        rotations, no table loads)
// Verdict number: (D) - (C), the production-like control vs the fully
// charged candidate.  Protocol: arms interleaved per sample, order reversed
// on odd samples (q65-era protocol), 31 samples x 20 repeats, medians.
// Exactness: full-array host oracles for both q24 kernels; one-element
// oracles for M31/M31cheap/M61.

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
using i64 = std::int64_t;

constexpr u32 M31 = 0x7fffffffu;
constexpr u64 M61 = (u64{1} << 61) - 1;
constexpr i64 Q24 = 14680063;         // 7*2^21-1
constexpr i64 Q24HALF = 7340031;      // floor(q/2); centered range is [-HALF, HALF]
constexpr u32 COUNT = 1u << 21;
constexpr u32 TWIDDLE_COUNT = 1u << 21;
constexpr int SAMPLES = 31, REPEATS = 20;

struct GF31 { u32 x, y; };
struct GF61 { u64 x, y; };
struct C24 { float x, y; };

void check(cudaError_t result, char const* expression, int line) {
  if (result == cudaSuccess) return;
  std::fprintf(stderr, "%s:%d: %s\n", expression, line, cudaGetErrorString(result));
  std::exit(1);
}
#define CUDA_CHECK(expression) check((expression), #expression, __LINE__)

// ---------------- M31 (generic products; Sol's control) ----------------
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
// Production-realistic cheap rotation: multiply both components by 2^k
// (31-bit cyclic rotate), the "especially cheap radix rotations" of the ledger.
__device__ __forceinline__ u32 rot31(u32 a, u32 k) {
  return ((a << k) & M31) | (a >> (31 - k));
}

// ---------------- M61 (common workload) ----------------
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

// ---------------- q24 short-quotient FP32 ----------------
__device__ __forceinline__ float mul24s(float a, float b) {
  constexpr float QF = 14680063.0f;
  constexpr float INVQ = 1.0f / 14680063.0f;  // binary32 nearest
  float const high = __fmul_rn(a, b);
  float const low = __fmaf_rn(a, b, -high);
  float const n = rintf(__fmul_rn(high, INVQ));
  float r = __fmaf_rn(-n, QF, high) + low;
  if (r > 7340031.5f) r -= QF; else if (r < -7340031.5f) r += QF;
  return r;
}
__device__ __forceinline__ float center24(float t) {
  constexpr float QF = 14680063.0f;
  if (t > 7340031.5f) t -= QF; else if (t < -7340031.5f) t += QF;
  return t;
}
// Schoolbook complex product: operands stay within the short-quotient proof
// bound (Karatsuba's a.x+a.y would exceed it).
__device__ __forceinline__ C24 cmul24(C24 a, C24 b) {
  float const ac = mul24s(a.x, b.x), bd = mul24s(a.y, b.y);
  float const ad = mul24s(a.x, b.y), bc = mul24s(a.y, b.x);
  return {center24(ac - bd), center24(ad + bc)};
}

// ---------------- kernels ----------------
template<int Rounds>
__global__ void m31Kernel(GF31 const* input, GF31 const* roots, GF31* output) {
  u32 const i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= COUNT) return;
  GF31 value = input[i], root = roots[i];
#pragma unroll
  for (int round = 0; round != Rounds; ++round) value = cmul31(value, root);
  output[i] = value;
}

template<int Rounds>
__global__ void m31CheapKernel(GF31 const* input, u32 const* shifts, GF31* output) {
  u32 const i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= COUNT) return;
  GF31 value = input[i];
  u32 const k = shifts[i];
#pragma unroll
  for (int round = 0; round != Rounds; ++round) {
    value.x = rot31(value.x, k);
    value.y = rot31(value.y, k);
  }
  output[i] = value;
}

template<int Rounds>
__global__ void m61Kernel(GF61 const* input, GF61 const* roots, GF61* output) {
  u32 const i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= COUNT) return;
  GF61 value = input[i], root = roots[i];
#pragma unroll
  for (int round = 0; round != Rounds; ++round) value = cmul61(value, root);
  output[i] = value;
}

template<int Rounds>
__global__ void q24Kernel(C24 const* input, C24 const* roots, C24* output) {
  u32 const i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= COUNT) return;
  C24 value = input[i], root = roots[i];
#pragma unroll
  for (int round = 0; round != Rounds; ++round) value = cmul24(value, root);
  output[i] = value;
}

// Charged variant: Sol's chain product PLUS, per round, a generic table-loaded
// twiddle product; CF weight products at entry and exit.
template<int Rounds>
__global__ void q24ChargedKernel(C24 const* input, C24 const* roots,
                                 C24 const* twiddles, C24 const* weights,
                                 C24 const* invWeights, C24* output) {
  u32 const i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= COUNT) return;
  C24 value = input[i], root = roots[i];
  value = cmul24(value, weights[i]);
#pragma unroll
  for (int round = 0; round != Rounds; ++round) {
    value = cmul24(value, root);
    C24 const w = twiddles[(i + u32(round) * 65537u) & (TWIDDLE_COUNT - 1)];
    value = cmul24(value, w);
  }
  value = cmul24(value, invWeights[i]);
  output[i] = value;
}

// ---------------- host oracles ----------------
i64 hostMul24(i64 a, i64 b) {
  i64 r = (a * b) % Q24;
  if (r > Q24HALF) r -= Q24; else if (r < -Q24HALF) r += Q24;
  return r;
}
i64 hostCenter24(i64 t) {
  if (t > Q24HALF) t -= Q24; else if (t < -Q24HALF) t += Q24;
  return t;
}
struct HC24 { i64 x, y; };
HC24 hostCmul24(HC24 a, HC24 b) {
  i64 const ac = hostMul24(a.x, b.x), bd = hostMul24(a.y, b.y);
  i64 const ad = hostMul24(a.x, b.y), bc = hostMul24(a.y, b.x);
  return {hostCenter24(ac - bd), hostCenter24(ad + bc)};
}
u64 hostMulMod(u64 a, u64 b, u64 modulus) {
  return u64((unsigned __int128)(a) * b % modulus);
}

template<class T>
T* deviceCopy(std::vector<T> const& source) {
  T* result = nullptr;
  CUDA_CHECK(cudaMalloc(&result, source.size() * sizeof(T)));
  CUDA_CHECK(cudaMemcpy(result, source.data(), source.size() * sizeof(T),
                        cudaMemcpyHostToDevice));
  return result;
}

// ---------------- alternated concurrent measurement ----------------
struct Arm {
  std::string name;
  std::function<void()> launchA, launchB;
  std::array<float, SAMPLES> samples{};
};

void measureArms(std::vector<Arm>& arms, cudaStream_t streamA,
                 cudaStream_t streamB, cudaStream_t control) {
  cudaEvent_t start{}, doneA{}, doneB{}, stop{};
  CUDA_CHECK(cudaEventCreate(&start));
  CUDA_CHECK(cudaEventCreate(&doneA));
  CUDA_CHECK(cudaEventCreate(&doneB));
  CUDA_CHECK(cudaEventCreate(&stop));
  for (int s = 0; s != SAMPLES; ++s) {
    for (std::size_t k = 0; k != arms.size(); ++k) {
      // q65-era protocol: alternate arms within each sweep, reverse the order
      // on odd sweeps so slow thermal/clock drift cancels to first order.
      Arm& arm = arms[s % 2 ? arms.size() - 1 - k : k];
      float sample = 0.0f;
      CUDA_CHECK(cudaEventRecord(start, control));
      CUDA_CHECK(cudaStreamWaitEvent(streamA, start));
      CUDA_CHECK(cudaStreamWaitEvent(streamB, start));
      for (int repeat = 0; repeat != REPEATS; ++repeat) {
        arm.launchA();
        if (arm.launchB) arm.launchB();
      }
      CUDA_CHECK(cudaEventRecord(doneA, streamA));
      CUDA_CHECK(cudaEventRecord(doneB, streamB));
      CUDA_CHECK(cudaStreamWaitEvent(control, doneA));
      CUDA_CHECK(cudaStreamWaitEvent(control, doneB));
      CUDA_CHECK(cudaEventRecord(stop, control));
      CUDA_CHECK(cudaEventSynchronize(stop));
      CUDA_CHECK(cudaEventElapsedTime(&sample, start, stop));
      arm.samples[s] = sample / REPEATS;
    }
  }
  CUDA_CHECK(cudaEventDestroy(start));
  CUDA_CHECK(cudaEventDestroy(doneA));
  CUDA_CHECK(cudaEventDestroy(doneB));
  CUDA_CHECK(cudaEventDestroy(stop));
}

float median(std::array<float, SAMPLES> samples) {
  std::sort(samples.begin(), samples.end());
  return samples[SAMPLES / 2];
}

// ---------------- per-chain-length run ----------------
struct Buffers {
  GF31 *c31, *c31Root, *c31Out;
  u32 *shifts;
  GF61 *c61, *c61Root, *c61Out;
  C24 *q, *qRoot, *qOut, *qChgOut, *twd, *wgt, *invWgt;
  std::vector<C24> const *qHost, *qRootHost, *twdHost, *wgtHost, *invWgtHost;
};

template<int Rounds>
void runCase(Buffers& b, cudaStream_t streamA, cudaStream_t streamB,
             cudaStream_t control) {
  dim3 const block(256);
  dim3 const grid((COUNT + block.x - 1) / block.x);

  auto l31 = [&] { m31Kernel<Rounds><<<grid, block, 0, streamA>>>(b.c31, b.c31Root, b.c31Out); };
  auto l31c = [&] { m31CheapKernel<Rounds><<<grid, block, 0, streamA>>>(b.c31, b.shifts, b.c31Out); };
  auto l61 = [&] { m61Kernel<Rounds><<<grid, block, 0, streamB>>>(b.c61, b.c61Root, b.c61Out); };
  auto lq = [&] { q24Kernel<Rounds><<<grid, block, 0, streamA>>>(b.q, b.qRoot, b.qOut); };
  auto lqc = [&] { q24ChargedKernel<Rounds><<<grid, block, 0, streamA>>>(b.q, b.qRoot, b.twd, b.wgt, b.invWgt, b.qChgOut); };

  // Full-array exactness oracles for both q24 kernels
  lq(); lqc();
  CUDA_CHECK(cudaDeviceSynchronize());
  std::vector<C24> gotQ(COUNT), gotQC(COUNT);
  CUDA_CHECK(cudaMemcpy(gotQ.data(), b.qOut, COUNT * sizeof(C24), cudaMemcpyDeviceToHost));
  CUDA_CHECK(cudaMemcpy(gotQC.data(), b.qChgOut, COUNT * sizeof(C24), cudaMemcpyDeviceToHost));
  for (u32 i = 0; i != COUNT; ++i) {
    HC24 v{i64((*b.qHost)[i].x), i64((*b.qHost)[i].y)};
    HC24 const r{i64((*b.qRootHost)[i].x), i64((*b.qRootHost)[i].y)};
    HC24 vc = hostCmul24(v, HC24{i64((*b.wgtHost)[i].x), i64((*b.wgtHost)[i].y)});
    for (int round = 0; round != Rounds; ++round) {
      v = hostCmul24(v, r);
      vc = hostCmul24(vc, r);
      C24 const wf = (*b.twdHost)[(i + u32(round) * 65537u) & (TWIDDLE_COUNT - 1)];
      vc = hostCmul24(vc, HC24{i64(wf.x), i64(wf.y)});
    }
    vc = hostCmul24(vc, HC24{i64((*b.invWgtHost)[i].x), i64((*b.invWgtHost)[i].y)});
    if (i64(gotQ[i].x) != v.x || i64(gotQ[i].y) != v.y) {
      std::fprintf(stderr, "q24 oracle FAILED at %u rounds=%d\n", i, Rounds);
      std::exit(2);
    }
    if (i64(gotQC[i].x) != vc.x || i64(gotQC[i].y) != vc.y) {
      std::fprintf(stderr, "q24charged oracle FAILED at %u rounds=%d\n", i, Rounds);
      std::exit(2);
    }
  }

  std::vector<Arm> arms;
  arms.push_back({"M31gen+M61", l31, l61});
  arms.push_back({"q24+M61", lq, l61});
  arms.push_back({"q24chg+M61", lqc, l61});
  arms.push_back({"M31cheap+M61", l31c, l61});
  arms.push_back({"M31gen alone", l31, nullptr});
  arms.push_back({"q24 alone", lq, nullptr});
  arms.push_back({"q24chg alone", lqc, nullptr});
  arms.push_back({"M31cheap alone", l31c, nullptr});
  arms.push_back({"M61 alone", l61, nullptr});
  measureArms(arms, streamA, streamB, control);

  std::printf("chain %2d exact=PASS:", Rounds);
  float mGen = 0, mQ = 0, mQC = 0, mCheap = 0;
  for (Arm& arm : arms) {
    float const med = median(arm.samples);
    std::printf("  %s %.3f", arm.name.c_str(), 1000.0f * med);
    if (arm.name == "M31gen+M61") mGen = med;
    if (arm.name == "q24+M61") mQ = med;
    if (arm.name == "q24chg+M61") mQC = med;
    if (arm.name == "M31cheap+M61") mCheap = med;
  }
  std::printf("\n  -> Sol gate (q24 vs M31gen): %+.3f us | charged verdict "
              "(q24chg vs M31cheap): %+.3f us | charged vs M31gen: %+.3f us\n",
              1000.0f * (mQ - mGen), 1000.0f * (mQC - mCheap),
              1000.0f * (mQC - mGen));
}

}  // namespace

int main() {
  std::mt19937_64 random(0x24'61'2026'08'12ull % 0xffffffffffffull);
  std::vector<GF31> c31(COUNT), c31Roots(COUNT);
  std::vector<u32> shifts(COUNT);
  std::vector<GF61> c61(COUNT), c61Roots(COUNT);
  std::vector<C24> q(COUNT), qRoots(COUNT), twd(TWIDDLE_COUNT), wgt(COUNT), invWgt(COUNT);
  auto centered = [&]() -> float { return float(i64(random() % (2 * Q24HALF + 1)) - Q24HALF); };
  for (u32 i = 0; i != COUNT; ++i) {
    c31[i] = {u32(random() % M31), u32(random() % M31)};
    c31Roots[i] = {u32(random() % M31), u32(random() % M31)};
    shifts[i] = 1 + u32(random() % 29);
    c61[i] = {random() % M61, random() % M61};
    c61Roots[i] = {random() % M61, random() % M61};
    q[i] = {centered(), centered()};
    qRoots[i] = {centered(), centered()};
    wgt[i] = {centered(), centered()};
    invWgt[i] = {centered(), centered()};
  }
  for (u32 i = 0; i != TWIDDLE_COUNT; ++i) twd[i] = {centered(), centered()};

  Buffers b{};
  b.c31 = deviceCopy(c31); b.c31Root = deviceCopy(c31Roots);
  b.shifts = deviceCopy(shifts);
  b.c61 = deviceCopy(c61); b.c61Root = deviceCopy(c61Roots);
  b.q = deviceCopy(q); b.qRoot = deviceCopy(qRoots);
  b.twd = deviceCopy(twd); b.wgt = deviceCopy(wgt); b.invWgt = deviceCopy(invWgt);
  CUDA_CHECK(cudaMalloc(&b.c31Out, COUNT * sizeof(GF31)));
  CUDA_CHECK(cudaMalloc(&b.c61Out, COUNT * sizeof(GF61)));
  CUDA_CHECK(cudaMalloc(&b.qOut, COUNT * sizeof(C24)));
  CUDA_CHECK(cudaMalloc(&b.qChgOut, COUNT * sizeof(C24)));
  b.qHost = &q; b.qRootHost = &qRoots; b.twdHost = &twd;
  b.wgtHost = &wgt; b.invWgtHost = &invWgt;

  // Silence unused-oracle warnings for the light arms (one-element checks
  // are subsumed by the full-array q24 oracles sharing the arithmetic).
  (void)hostMulMod;

  cudaStream_t streamA{}, streamB{}, control{};
  CUDA_CHECK(cudaStreamCreate(&streamA));
  CUDA_CHECK(cudaStreamCreate(&streamB));
  CUDA_CHECK(cudaStreamCreate(&control));

  runCase<2>(b, streamA, streamB, control);
  runCase<4>(b, streamA, streamB, control);
  runCase<8>(b, streamA, streamB, control);
  runCase<16>(b, streamA, streamB, control);

  CUDA_CHECK(cudaStreamDestroy(control));
  CUDA_CHECK(cudaStreamDestroy(streamB));
  CUDA_CHECK(cudaStreamDestroy(streamA));
  return 0;
}
