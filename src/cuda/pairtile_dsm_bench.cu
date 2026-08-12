// Approach-B gate (2026-08-12): can an 8-block cluster with cooperative
// coalesced loads + DSM scatter beat the 16x-amplified direct reads that
// closed FUSED31?  Isolates the load/scatter physics of a pair-resident
// fused kernel on the production INPLACE M31 layout (16-width interleave,
// XOR swizzle).  Arms:
//   A direct:  block gathers 2 far-apart width columns straight from the
//              interleaved plane (FUSED31's pattern, ~16x amplification)
//   B cluster: 8-block cluster covers 16 consecutive widths; 2048 threads
//              sweep the segment range linearly (coalesced) and DSM-scatter
//              each element to the owning block's tile
//   C linear:  pure coalesced sweep, no scatter (lower bound)
// Verdict number: B vs A.  If B ~ A or worse, approach B closes.

#include <cuda_runtime.h>
#include <cooperative_groups.h>
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <vector>

namespace cg = cooperative_groups;
using u32 = std::uint32_t;
using u64 = std::uint64_t;

// production GF31 plane constants (WIDTH=512, MIDDLE=8, SMALL_HEIGHT=512)
constexpr u32 WIDTH = 512, MIDDLE = 8, SH = 512;
constexpr u32 SIZEBLK = SH;                    // 512
constexpr u32 SIZEW = 16 * SIZEBLK + 16;       // 8208
constexpr u32 SIZEM = WIDTH / 16 * SIZEW + 16; // 262672
constexpr u32 PLANE = MIDDLE * SIZEM;          // GF31 units
struct GF31 { u32 x, y; };
__host__ __device__ inline u32 addrOf(u32 w, u32 m, u32 h) {
  return (w / 16) * SIZEW + m * SIZEM + (w % 16) * SIZEBLK + (((h / 16) ^ (w % 16)) * 16) + h % 16;
}

void check(cudaError_t r, int line) {
  if (r != cudaSuccess) { std::fprintf(stderr, "%d: %s\n", line, cudaGetErrorString(r)); std::exit(1); }
}
#define CK(x) check((x), __LINE__)

// Arm A: direct gather of widths {b, WIDTH-1-b} (far apart like Hermitian pairs)
__global__ void directKernel(GF31 const* plane, u64* out) {
  __shared__ GF31 tile[2 * MIDDLE * 64];   // per-m slice staging (looped over h-chunks)
  u32 const b = blockIdx.x;
  u32 const w1 = b, w2 = WIDTH - 1 - b;
  u64 acc = 0;
  for (u32 slot = 0; slot < 2; ++slot) {
    u32 const w = slot ? w2 : w1;
    for (u32 hBase = 0; hBase < SH; hBase += 64) {   // chunked to bound shared
      for (u32 t = threadIdx.x; t < MIDDLE * 64; t += blockDim.x) {
        u32 const m = t / 64, h = hBase + t % 64;
        tile[slot * MIDDLE * 64 + t] = plane[addrOf(w, m, h)];
      }
      __syncthreads();
      for (u32 t = threadIdx.x; t < MIDDLE * 64; t += blockDim.x)
        acc ^= (u64(tile[slot * MIDDLE * 64 + t].x) << 32) | tile[slot * MIDDLE * 64 + t].y;
      __syncthreads();
    }
  }
  if (threadIdx.x == 0) out[b] = acc;
}

// Arm B: 8-block cluster, coalesced sweep + DSM scatter
__global__ void __cluster_dims__(8, 1, 1) clusterKernel(GF31 const* plane, u64* out) {
  __shared__ GF31 tile[2 * MIDDLE * 64];   // this block's 2 widths, chunked by h
  cg::cluster_group cluster = cg::this_cluster();
  u32 const rank = cluster.block_rank();
  u32 const c = blockIdx.x / 8;            // cluster id: widths [16c, 16c+16)
  for (u32 hBase = 0; hBase < SH; hBase += 64) {
    // cooperative coalesced sweep of the cluster's 16-width slab for this
    // h-chunk: elements (w in 16c..16c+15, m, h in hBase..hBase+63).
    // Linear index t enumerates them in ADDRESS order for coalescing:
    // addr = base + m*SIZEM + lane16*SIZEBLK + swiz16*16 + low, so sweep
    // (m, lane16, hi4=h/16 slice, low) with 16*16 contiguous inner span.
    for (u32 t = threadIdx.x + rank * blockDim.x; t < MIDDLE * 16 * 4 * 16; t += 8 * blockDim.x) {
      u32 const low = t % 16;              // h % 16
      u32 const sw = (t / 16) % 4;         // swizzled h/16 slice within chunk
      u32 const lane = (t / 64) % 16;      // w % 16
      u32 const m = t / 1024;
      u32 const hHi = sw ^ (lane & 3) ^ (lane >> 2); // partial unswizzle mix
      u32 const h = hBase + ((sw) ^ 0) * 16 + low;   // NOTE: swizzle handled via addr recompute below
      u32 const w = 16 * c + lane;
      // recompute the true h for this swizzled slot: slot sw holds h/16 = sw ^ (w%16) restricted to chunk
      u32 const hh = hBase + (((sw + hBase / 16) ^ (w % 16)) % 4) * 16 + low;  // approximate inverse within 4-slice chunk
      (void) h; (void) hHi;
      GF31 const v = plane[(w / 16) * SIZEW + m * SIZEM + (w % 16) * SIZEBLK + (((hh / 16) ^ (w % 16)) * 16) + low];
      u32 const owner = lane / 2;
      GF31* remote = (GF31*) cluster.map_shared_rank(tile, owner);
      remote[(lane % 2) * MIDDLE * 64 + m * 64 + (hh - hBase)] = v;
    }
    cluster.sync();
    u64 acc = 0;
    for (u32 t = threadIdx.x; t < 2 * MIDDLE * 64; t += blockDim.x)
      acc ^= (u64(tile[t].x) << 32) | tile[t].y;
    if (threadIdx.x == 0) out[blockIdx.x] ^= acc;
    cluster.sync();
  }
}

// Arm C: pure linear sweep (coalesced lower bound)
__global__ void linearKernel(GF31 const* plane, u64* out) {
  u64 acc = 0;
  for (u32 t = blockIdx.x * blockDim.x + threadIdx.x; t < PLANE; t += gridDim.x * blockDim.x) {
    acc ^= (u64(plane[t].x) << 32) | plane[t].y;
  }
  if ((threadIdx.x & 31) == 0) out[blockIdx.x * (blockDim.x / 32) + threadIdx.x / 32] = acc;
}

int main() {
  std::vector<GF31> host(PLANE);
  for (u32 i = 0; i < PLANE; ++i) host[i] = {i * 2654435761u, i ^ 0x9e3779b9u};
  GF31* plane; u64* out;
  CK(cudaMalloc(&plane, PLANE * sizeof(GF31)));
  CK(cudaMemcpy(plane, host.data(), PLANE * sizeof(GF31), cudaMemcpyHostToDevice));
  CK(cudaMalloc(&out, 65536 * sizeof(u64)));
  CK(cudaMemset(out, 0, 65536 * sizeof(u64)));

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

  float const a = median([&] { directKernel<<<256, 256>>>(plane, out); });
  float const b = median([&] { clusterKernel<<<256, 256>>>(plane, out); });
  float const c = median([&] { linearKernel<<<512, 256>>>(plane, out); });
  CK(cudaGetLastError());
  std::printf("direct(16x-amp) %.3f us | cluster-DSM %.3f us | linear bound %.3f us | cluster/direct %.3f\n",
              1000 * a, 1000 * b, 1000 * c, b / a);
  return 0;
}
