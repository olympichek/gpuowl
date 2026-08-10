// Exact resource-overlap gate for a 3M two-field architecture:
// GF((2^61-1)^2) x GF(q^2), q=2^61-72,351,745.
//
// The control is the incumbent 4M M31/M61 population.  The candidate prime
// is the smallest-complement prime found with q+1 divisible by 2^20, a scalar
// cube root, and an exact 3M Crandall--Fagin root-of-two weight.  This gate
// includes the 25% shorter candidate population and concurrent field streams.

#include <cuda_runtime.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <vector>

namespace {

using u32 = std::uint32_t;
using u64 = std::uint64_t;
using u128 = unsigned __int128;

constexpr u32 M31 = 0x7fffffffu;
constexpr u64 M61 = (u64{1} << 61) - 1;
constexpr u64 Q61 = 2'305'843'009'141'342'207ull;
constexpr u64 Q61_NEG_INV = 0x501298fffbb00001ull;
constexpr u64 Q61_R = 578'813'960ull;
constexpr u64 Q61_R_INV = 721'230'287'205'085'132ull;
constexpr u32 CONTROL_COUNT = 1u << 21;
constexpr u32 CANDIDATE_COUNT = 3u << 19;

static_assert(Q61 == (u64{1} << 61) - 72'351'745);
static_assert((Q61 + 1) % (u64{1} << 20) == 0);
static_assert((Q61 - 1) % 3 == 0);
static_assert(u64(u128(Q61) * Q61_NEG_INV) == u64(-1));

struct GF31 { u32 x, y; };
struct GF61 { u64 x, y; };

void check(cudaError_t result, char const* expression, int line) {
  if (result == cudaSuccess) return;
  std::fprintf(stderr, "%s:%d: %s\n", expression, line,
               cudaGetErrorString(result));
  std::exit(1);
}
#define CUDA_CHECK(expression) check((expression), #expression, __LINE__)

__device__ __forceinline__ u32 add31(u32 a, u32 b) {
  u32 result = a + b;
  result = (result & M31) + (result >> 31);
  return result == M31 ? 0 : result;
}
__device__ __forceinline__ u32 sub31(u32 a, u32 b) {
  return a >= b ? a - b : a + M31 - b;
}
__device__ __forceinline__ u32 mul31(u32 a, u32 b) {
  u64 const product = u64(a) * b;
  u64 result = (product & M31) + (product >> 31);
  result = (result & M31) + (result >> 31);
  return result == M31 ? 0 : u32(result);
}
__device__ __forceinline__ GF31 cmul31(GF31 a, GF31 b) {
  u32 const k1 = mul31(b.x, add31(a.x, a.y));
  u32 const k2 = mul31(a.x, sub31(b.y, b.x));
  u32 const k3 = mul31(a.y, add31(b.y, b.x));
  return {sub31(k1, k3), add31(k1, k2)};
}

__device__ __forceinline__ u64 add61(u64 a, u64 b) {
  u64 result = a + b;
  return result >= M61 ? result - M61 : result;
}
__device__ __forceinline__ u64 sub61(u64 a, u64 b) {
  return a >= b ? a - b : M61 - (b - a);
}
__device__ __forceinline__ u64 mul61(u64 a, u64 b) {
  u64 const lo = a * b;
  u64 const hi = __umul64hi(a, b);
  u64 result = (lo & M61) + (lo >> 61) + (hi << 3);
  if (result >= M61) result -= M61;
  if (result >= M61) result -= M61;
  return result;
}
__device__ __forceinline__ GF61 cmul61(GF61 a, GF61 b) {
  u64 const k1 = mul61(b.x, add61(a.x, a.y));
  u64 const k2 = mul61(a.x, sub61(b.y, b.x));
  u64 const k3 = mul61(a.y, add61(b.y, b.x));
  return {sub61(k1, k3), add61(k1, k2)};
}

__device__ __forceinline__ u64 addQ(u64 a, u64 b) {
  u64 result = a + b;
  return result >= Q61 ? result - Q61 : result;
}
__device__ __forceinline__ u64 subQ(u64 a, u64 b) {
  return a >= b ? a - b : Q61 - (b - a);
}
// R=2^64 Montgomery REDC.  The low word of multiplier*q cancels the
// product low word, so its carry into the high word is productLo!=0.
__device__ __forceinline__ u64 mulQ(u64 a, u64 b) {
  u64 const productLo = a * b;
  u64 const productHi = __umul64hi(a, b);
  u64 const multiplier = productLo * Q61_NEG_INV;
  u64 result = productHi + __umul64hi(multiplier, Q61) +
               (productLo != 0);
  return result >= Q61 ? result - Q61 : result;
}
__device__ __forceinline__ GF61 cmulQ(GF61 a, GF61 b) {
  u64 const k1 = mulQ(b.x, addQ(a.x, a.y));
  u64 const k2 = mulQ(a.x, subQ(b.y, b.x));
  u64 const k3 = mulQ(a.y, addQ(b.y, b.x));
  return {subQ(k1, k3), addQ(k1, k2)};
}

template<int Rounds>
__global__ void m31Kernel(GF31 const* input, GF31 const* roots, GF31* output,
                          u32 count) {
  u32 const i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= count) return;
  GF31 value = input[i], root = roots[i];
#pragma unroll
  for (int round = 0; round != Rounds; ++round) value = cmul31(value, root);
  output[i] = value;
}

template<int Rounds>
__global__ void m61Kernel(GF61 const* input, GF61 const* roots, GF61* output,
                          u32 count) {
  u32 const i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= count) return;
  GF61 value = input[i], root = roots[i];
#pragma unroll
  for (int round = 0; round != Rounds; ++round) value = cmul61(value, root);
  output[i] = value;
}

template<int Rounds>
__global__ void q61Kernel(GF61 const* input, GF61 const* roots, GF61* output) {
  u32 const i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= CANDIDATE_COUNT) return;
  GF61 value = input[i], root = roots[i];
#pragma unroll
  for (int round = 0; round != Rounds; ++round) value = cmulQ(value, root);
  output[i] = value;
}

template<class T>
T* deviceCopy(std::vector<T> const& source) {
  T* result = nullptr;
  CUDA_CHECK(cudaMalloc(&result, source.size() * sizeof(T)));
  CUDA_CHECK(cudaMemcpy(result, source.data(), source.size() * sizeof(T),
                        cudaMemcpyHostToDevice));
  return result;
}

template<class LaunchA, class LaunchB>
float concurrentMedian(LaunchA launchA, LaunchB launchB,
                       cudaStream_t streamA, cudaStream_t streamB,
                       cudaStream_t control) {
  constexpr int Samples = 31, Repeats = 20;
  cudaEvent_t start{}, doneA{}, doneB{}, stop{};
  CUDA_CHECK(cudaEventCreate(&start));
  CUDA_CHECK(cudaEventCreate(&doneA));
  CUDA_CHECK(cudaEventCreate(&doneB));
  CUDA_CHECK(cudaEventCreate(&stop));
  std::array<float, Samples> samples{};
  for (float& sample : samples) {
    CUDA_CHECK(cudaEventRecord(start, control));
    CUDA_CHECK(cudaStreamWaitEvent(streamA, start));
    CUDA_CHECK(cudaStreamWaitEvent(streamB, start));
    for (int repeat = 0; repeat != Repeats; ++repeat) {
      launchA();
      launchB();
    }
    CUDA_CHECK(cudaEventRecord(doneA, streamA));
    CUDA_CHECK(cudaEventRecord(doneB, streamB));
    CUDA_CHECK(cudaStreamWaitEvent(control, doneA));
    CUDA_CHECK(cudaStreamWaitEvent(control, doneB));
    CUDA_CHECK(cudaEventRecord(stop, control));
    CUDA_CHECK(cudaEventSynchronize(stop));
    CUDA_CHECK(cudaEventElapsedTime(&sample, start, stop));
    sample /= Repeats;
  }
  std::sort(samples.begin(), samples.end());
  CUDA_CHECK(cudaEventDestroy(start));
  CUDA_CHECK(cudaEventDestroy(doneA));
  CUDA_CHECK(cudaEventDestroy(doneB));
  CUDA_CHECK(cudaEventDestroy(stop));
  return samples[Samples / 2];
}

u64 hostMul(u64 a, u64 b, u64 modulus) {
  return u64(u128(a) * b % modulus);
}
GF61 hostCmul(GF61 a, GF61 b, u64 modulus) {
  u64 const ac = hostMul(a.x, b.x, modulus);
  u64 const bd = hostMul(a.y, b.y, modulus);
  u64 const real = ac >= bd ? ac - bd : modulus - (bd - ac);
  u64 const imag = (hostMul(a.x, b.y, modulus) +
                    hostMul(a.y, b.x, modulus)) % modulus;
  return {real, imag};
}

template<int Rounds>
void runCase(GF31 const* control31, GF31 const* controlRoot31,
             GF31* controlOut31, GF61 const* control61,
             GF61 const* controlRoot61, GF61* controlOut61,
             GF61 const* candidate61, GF61 const* candidateRoot61,
             GF61* candidateOut61, GF61 const* candidateQ,
             GF61 const* candidateRootQ, GF61* candidateOutQ,
             GF61 candidateQ0, GF61 candidateRootQ0,
             cudaStream_t streamA, cudaStream_t streamB,
             cudaStream_t control) {
  dim3 const block(256);
  dim3 const controlGrid((CONTROL_COUNT + block.x - 1) / block.x);
  dim3 const candidateGrid((CANDIDATE_COUNT + block.x - 1) / block.x);
  auto launchControl31 = [&] {
    m31Kernel<Rounds><<<controlGrid, block, 0, streamA>>>(
        control31, controlRoot31, controlOut31, CONTROL_COUNT);
  };
  auto launchControl61 = [&] {
    m61Kernel<Rounds><<<controlGrid, block, 0, streamB>>>(
        control61, controlRoot61, controlOut61, CONTROL_COUNT);
  };
  auto launchCandidate61 = [&] {
    m61Kernel<Rounds><<<candidateGrid, block, 0, streamA>>>(
        candidate61, candidateRoot61, candidateOut61, CANDIDATE_COUNT);
  };
  auto launchCandidateQ = [&] {
    q61Kernel<Rounds><<<candidateGrid, block, 0, streamB>>>(
        candidateQ, candidateRootQ, candidateOutQ);
  };

  launchCandidateQ();
  CUDA_CHECK(cudaStreamSynchronize(streamB));
  GF61 got{};
  CUDA_CHECK(cudaMemcpy(&got, candidateOutQ, sizeof(got),
                        cudaMemcpyDeviceToHost));
  GF61 expected = candidateQ0;
  for (int round = 0; round != Rounds; ++round)
    expected = hostCmul(expected, candidateRootQ0, Q61);
  // Inputs are Montgomery values, so each host product removes one R.
  for (int round = 0; round != Rounds; ++round) {
    expected.x = hostMul(expected.x, Q61_R_INV, Q61);
    expected.y = hostMul(expected.y, Q61_R_INV, Q61);
  }
  if (got.x != expected.x || got.y != expected.y) {
    std::fprintf(stderr, "q61 exact validation failed for %d rounds\n", Rounds);
    std::exit(2);
  }

  float const incumbent = concurrentMedian(
      launchControl31, launchControl61, streamA, streamB, control);
  float const candidate = concurrentMedian(
      launchCandidate61, launchCandidateQ, streamA, streamB, control);
  std::printf("%2d products exact=PASS: 4M M31+M61 %.3f us, "
              "3M M61+q61 %.3f us, candidate/control %.4f\n",
              Rounds, 1000.0f * incumbent, 1000.0f * candidate,
              candidate / incumbent);
}

}  // namespace

int main() {
  // The two algebraic conditions omitted by the compile-time congruences.
  if (hostMul(Q61, Q61_NEG_INV, u64{1} << 63) != ((u64{1} << 63) - 1) ||
      hostMul(1, 1, Q61) != 1) return 3;

  std::mt19937_64 random(0x3'61'61'03'00ull);
  std::vector<GF31> c31(CONTROL_COUNT), c31Roots(CONTROL_COUNT);
  std::vector<GF61> c61(CONTROL_COUNT), c61Roots(CONTROL_COUNT);
  std::vector<GF61> n61(CANDIDATE_COUNT), n61Roots(CANDIDATE_COUNT);
  std::vector<GF61> q(CANDIDATE_COUNT), qRoots(CANDIDATE_COUNT);
  for (u32 i = 0; i != CONTROL_COUNT; ++i) {
    c31[i] = {u32(random() % M31), u32(random() % M31)};
    c31Roots[i] = {u32(random() % M31), u32(random() % M31)};
    c61[i] = {random() % M61, random() % M61};
    c61Roots[i] = {random() % M61, random() % M61};
  }
  for (u32 i = 0; i != CANDIDATE_COUNT; ++i) {
    n61[i] = {random() % M61, random() % M61};
    n61Roots[i] = {random() % M61, random() % M61};
    GF61 const normal{random() % Q61, random() % Q61};
    GF61 const normalRoot{random() % Q61, random() % Q61};
    q[i] = {hostMul(normal.x, Q61_R, Q61),
            hostMul(normal.y, Q61_R, Q61)};
    qRoots[i] = {hostMul(normalRoot.x, Q61_R, Q61),
                 hostMul(normalRoot.y, Q61_R, Q61)};
  }

  GF31* c31D = deviceCopy(c31), *c31RootD = deviceCopy(c31Roots), *c31OutD{};
  GF61* c61D = deviceCopy(c61), *c61RootD = deviceCopy(c61Roots), *c61OutD{};
  GF61* n61D = deviceCopy(n61), *n61RootD = deviceCopy(n61Roots), *n61OutD{};
  GF61* qD = deviceCopy(q), *qRootD = deviceCopy(qRoots), *qOutD{};
  CUDA_CHECK(cudaMalloc(&c31OutD, c31.size() * sizeof(GF31)));
  CUDA_CHECK(cudaMalloc(&c61OutD, c61.size() * sizeof(GF61)));
  CUDA_CHECK(cudaMalloc(&n61OutD, n61.size() * sizeof(GF61)));
  CUDA_CHECK(cudaMalloc(&qOutD, q.size() * sizeof(GF61)));
  cudaStream_t streamA{}, streamB{}, control{};
  CUDA_CHECK(cudaStreamCreate(&streamA));
  CUDA_CHECK(cudaStreamCreate(&streamB));
  CUDA_CHECK(cudaStreamCreate(&control));

#define RUN_CASE(rounds)                                                        \
  runCase<rounds>(c31D, c31RootD, c31OutD, c61D, c61RootD, c61OutD,           \
                  n61D, n61RootD, n61OutD, qD, qRootD, qOutD, q[0], qRoots[0], \
                  streamA, streamB, control)
  RUN_CASE(2);
  RUN_CASE(4);
  RUN_CASE(8);
  RUN_CASE(16);
#undef RUN_CASE

  CUDA_CHECK(cudaStreamDestroy(control));
  CUDA_CHECK(cudaStreamDestroy(streamB));
  CUDA_CHECK(cudaStreamDestroy(streamA));
  CUDA_CHECK(cudaFree(qOutD)); CUDA_CHECK(cudaFree(qRootD)); CUDA_CHECK(cudaFree(qD));
  CUDA_CHECK(cudaFree(n61OutD)); CUDA_CHECK(cudaFree(n61RootD)); CUDA_CHECK(cudaFree(n61D));
  CUDA_CHECK(cudaFree(c61OutD)); CUDA_CHECK(cudaFree(c61RootD)); CUDA_CHECK(cudaFree(c61D));
  CUDA_CHECK(cudaFree(c31OutD)); CUDA_CHECK(cudaFree(c31RootD)); CUDA_CHECK(cudaFree(c31D));
  return 0;
}
