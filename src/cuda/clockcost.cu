// clockcost.cu — E1 "clock-cost table" harness (power-aware codegen track).
//
// Measures per-instruction-class power on GB202 (RTX PRO 6000 Blackwell SE):
// each kernel is an unrolled inner loop of ONE SASS opcode class (verified via
// cuobjdump), at full occupancy.  LO- vs HI-operand-toggle variants use the
// SAME kernel binary and differ ONLY in runtime seed data, so any power delta
// is pure operand-toggle activity (arXiv:2409.18324 effect).
//
// Design notes (learned from SASS v1):
//  - chain inits must be THREAD-VARYING or ptxas promotes the whole chain to
//    the uniform datapath (ULOP3/UIMAD/UFFMA) — wrong ALU;
//  - lo/hi must differ only in runtime data: asm volatile does not survive
//    into PTX, and ptxas folds add-0 / mul-1 / mov chains;
//  - shared-mem chains must load from a *different* address than they store
//    (neighbor slot) or ptxas store-to-load-forwards the LDS away.
//
// Observables while a kernel streams for N seconds:
//   - sustained SM clock (NVML, 20 Hz)               -> J/cycle sensor at a binding power cap
//   - nvmlDeviceGetTotalEnergyConsumption delta      -> integrated watts
//   - instantaneous power samples (cross-check), temperature
//   - ops/s from launch count                        -> pJ/op
//
// Protocols (set outside via nvidia-smi):
//   P1: -lgc <f>,<f>  fixed clock  -> watts table (J/cycle = W/f)
//   P2: -pl  <cap>    binding cap  -> sustained-MHz table
//
// Usage: clockcost <run|list|all|header> [seconds=75]

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <vector>
#include <string>
#include <thread>
#include <atomic>
#include <algorithm>
#include <chrono>
#include <cuda_runtime.h>
#include <nvml.h>

typedef unsigned u32;
typedef unsigned long long u64;

#define CUDA_CHECK(x) do { cudaError_t e = (x); if (e != cudaSuccess) { \
  fprintf(stderr, "CUDA error %s at %s:%d\n", cudaGetErrorString(e), __FILE__, __LINE__); exit(1); } } while (0)
#define NVML_CHECK(x) do { nvmlReturn_t e = (x); if (e != NVML_SUCCESS) { \
  fprintf(stderr, "NVML error %s at %s:%d\n", nvmlErrorString(e), __FILE__, __LINE__); exit(1); } } while (0)

// 16x repetition; with 4 chains inside = 64 class-ops per loop iteration.
#define R16(X) X X X X X X X X X X X X X X X X

#define TV_INIT \
  u32 tv = (threadIdx.x + blockIdx.x * blockDim.x) * 0x9E3779B9u ^ 0x243F6A88u; \
  (void)tv;

#define TAILGUARD(expr) if ((expr) == 0xdeadbeefu && threadIdx.x == 1023u) out[blockIdx.x] = (expr)

// ---------------------------------------------------------------- kernels ----

// IADD3: a += b.  lo: seedA=0 -> b=0 (bits static).  hi: seedA=~0 -> b=0x55555555 (carry ripple).
__global__ void k_iadd3(u64 iters, u32 seedA, u32 seedB, u32* out) {
  TV_INIT;
  u32 a0 = tv, a1 = tv ^ 0x33333333u, a2 = ~tv, a3 = tv + 0x77777777u;
  u32 b = seedA & 0x55555555u;
  for (u64 i = 0; i < iters; ++i) {
    R16(asm volatile("add.u32 %0,%0,%4;\n\tadd.u32 %1,%1,%4;\n\tadd.u32 %2,%2,%4;\n\tadd.u32 %3,%3,%4;"
                     : "+r"(a0), "+r"(a1), "+r"(a2), "+r"(a3) : "r"(b));)
  }
  TAILGUARD(a0 + a1 + a2 + a3 + seedB);
}

// LOP3 (XOR3): a ^= b ^ c.  lo: b=c=0 (static).  hi: b^c=0xFFFFFFFF (every bit flips per op).
__global__ void k_lop3(u64 iters, u32 seedA, u32 seedB, u32* out) {
  TV_INIT;
  u32 a0 = tv, a1 = tv ^ 0x33333333u, a2 = ~tv, a3 = tv + 0x77777777u;
  u32 b = seedA & 0xAAAAAAAAu, c = seedA & 0x55555555u;
  for (u64 i = 0; i < iters; ++i) {
    R16(asm volatile("lop3.b32 %0,%0,%4,%5,0x96;\n\tlop3.b32 %1,%1,%4,%5,0x96;\n\t"
                     "lop3.b32 %2,%2,%4,%5,0x96;\n\tlop3.b32 %3,%3,%4,%5,0x96;"
                     : "+r"(a0), "+r"(a1), "+r"(a2), "+r"(a3) : "r"(b), "r"(c));)
  }
  TAILGUARD(a0 + a1 + a2 + a3 + seedB);
}

// SHF rotate-by-1.  lo: seedB=0 -> rotating zeros (static).  hi: seedB=~0 -> rotating tv (toggle).
__global__ void k_shf(u64 iters, u32 seedA, u32 seedB, u32* out) {
  TV_INIT;
  u32 a0 = tv & seedB, a1 = (tv ^ 0x33333333u) & seedB, a2 = ~tv & seedB, a3 = (tv + 0x77777777u) & seedB;
  u32 s = (seedA & 0x1Fu) | 1u;
  for (u64 i = 0; i < iters; ++i) {
    R16(asm volatile("shf.l.wrap.b32 %0,%0,%0,%4;\n\tshf.l.wrap.b32 %1,%1,%1,%4;\n\t"
                     "shf.l.wrap.b32 %2,%2,%2,%4;\n\tshf.l.wrap.b32 %3,%3,%3,%4;"
                     : "+r"(a0), "+r"(a1), "+r"(a2), "+r"(a3) : "r"(s));)
  }
  TAILGUARD(a0 + a1 + a2 + a3);
}

// IMAD 32-bit: a = a*b + c.  lo: b=1, c=0 (static).  hi: b=0x9E3779B1, c=0x7F4A7C15 (avalanche).
__global__ void k_imad32(u64 iters, u32 seedA, u32 seedB, u32* out) {
  TV_INIT;
  u32 a0 = tv, a1 = tv ^ 0x33333333u, a2 = ~tv, a3 = tv + 0x77777777u;
  u32 b = (seedA & 0xFFFFFFFEu) | 1u, c = seedB;
  for (u64 i = 0; i < iters; ++i) {
    R16(asm volatile("mad.lo.u32 %0,%0,%4,%5;\n\tmad.lo.u32 %1,%1,%4,%5;\n\t"
                     "mad.lo.u32 %2,%2,%4,%5;\n\tmad.lo.u32 %3,%3,%4,%5;"
                     : "+r"(a0), "+r"(a1), "+r"(a2), "+r"(a3) : "r"(b), "r"(c));)
  }
  TAILGUARD(a0 + a1 + a2 + a3);
}

// IMAD.WIDE static: A += x*y with runtime-constant x,y (x=seedA=0 -> A static; the lo floor).
__global__ void k_imadwide_st(u64 iters, u32 seedA, u32 seedB, u32* out) {
  TV_INIT;
  u64 A0 = tv, A1 = tv ^ 0x33333333u, A2 = ~(u64)tv, A3 = tv + 0x77777777u;
  u32 x = seedA, y = seedB | 1u;
  for (u64 i = 0; i < iters; ++i) {
    R16(asm volatile("mad.wide.u32 %0,%4,%5,%0;\n\tmad.wide.u32 %1,%4,%5,%1;\n\t"
                     "mad.wide.u32 %2,%4,%5,%2;\n\tmad.wide.u32 %3,%4,%5,%3;"
                     : "+l"(A0), "+l"(A1), "+l"(A2), "+l"(A3) : "r"(x), "r"(y));)
  }
  TAILGUARD((u32)(A0 + A1 + A2 + A3));
}

// IMAD.WIDE feedback: x = lo32(A) & seedA; A += x*y.
// lo: seedA=0 -> x=0 (multiplier operand static).  hi: seedA=~0 -> avalanche through the multiplier.
__global__ void k_imadwide_fb(u64 iters, u32 seedA, u32 seedB, u32* out) {
  TV_INIT;
  u64 A0 = tv, A1 = tv ^ 0x33333333u, A2 = ~(u64)tv, A3 = tv + 0x77777777u;
  u32 y = seedB | 1u;
  u32 x0, x1, x2, x3;
  for (u64 i = 0; i < iters; ++i) {
    R16(asm volatile("cvt.u32.u64 %4,%0;\n\tand.b32 %4,%4,%9;\n\tmad.wide.u32 %0,%4,%8,%0;\n\t"
                     "cvt.u32.u64 %5,%1;\n\tand.b32 %5,%5,%9;\n\tmad.wide.u32 %1,%5,%8,%1;\n\t"
                     "cvt.u32.u64 %6,%2;\n\tand.b32 %6,%6,%9;\n\tmad.wide.u32 %2,%6,%8,%2;\n\t"
                     "cvt.u32.u64 %7,%3;\n\tand.b32 %7,%7,%9;\n\tmad.wide.u32 %3,%7,%8,%3;"
                     : "+l"(A0), "+l"(A1), "+l"(A2), "+l"(A3),
                       "=r"(x0), "=r"(x1), "=r"(x2), "=r"(x3) : "r"(y), "r"(seedA));)
  }
  TAILGUARD((u32)(A0 + A1 + A2 + A3) + x0 + x1 + x2 + x3);
}

// 64-bit add (IADD3 + IADD3.X pair): A += b, b = seedA | seedB<<32.
__global__ void k_add64(u64 iters, u32 seedA, u32 seedB, u32* out) {
  TV_INIT;
  u64 A0 = tv, A1 = tv ^ 0x33333333u, A2 = ~(u64)tv, A3 = tv + 0x77777777u;
  u64 b = (u64)seedA | ((u64)seedB << 32);
  for (u64 i = 0; i < iters; ++i) {
    R16(asm volatile("add.u64 %0,%0,%4;\n\tadd.u64 %1,%1,%4;\n\tadd.u64 %2,%2,%4;\n\tadd.u64 %3,%3,%4;"
                     : "+l"(A0), "+l"(A1), "+l"(A2), "+l"(A3) : "l"(b));)
  }
  TAILGUARD((u32)(A0 + A1 + A2 + A3));
}

// FFMA: a = a*b + c (b,c bit-cast from seeds).
// lo: b=1.0, c=0.0 (static).  hi: b=-1.0, c=1.0 -> a <- 1-a oscillates forever (mantissa+sign toggle).
__global__ void k_ffma(u64 iters, u32 seedA, u32 seedB, u32* out) {
  TV_INIT;
  float a0 = __uint_as_float((tv & 0x007FFFFFu) | 0x3F000000u);          // [0.5, 1)
  float a1 = __uint_as_float(((tv >> 3) & 0x007FFFFFu) | 0x3E800000u);
  float a2 = __uint_as_float(((tv >> 5) & 0x007FFFFFu) | 0x3E000000u);
  float a3 = __uint_as_float(((tv >> 7) & 0x007FFFFFu) | 0x3D800000u);
  float b = __uint_as_float(seedA), c = __uint_as_float(seedB);
  for (u64 i = 0; i < iters; ++i) {
    R16(asm volatile("fma.rn.f32 %0,%0,%4,%5;\n\tfma.rn.f32 %1,%1,%4,%5;\n\t"
                     "fma.rn.f32 %2,%2,%4,%5;\n\tfma.rn.f32 %3,%3,%4,%5;"
                     : "+f"(a0), "+f"(a1), "+f"(a2), "+f"(a3) : "f"(b), "f"(c));)
  }
  TAILGUARD((u32)(a0 + a1 + a2 + a3));
}

// LDS/STS: ld from neighbor slot (s^4: provably != s, defeats store-to-load fwd),
// xor with seedB, st back.  lo: seedB=0 (values settle static).  hi: seedB=~0 (slots alternate v/~v).
// Class = 2 shared ops + 1 LOP3 per step (documented dilution, identical lo/hi).
__global__ void k_lds(u64 iters, u32 seedA, u32 seedB, u32* out) {
  __shared__ u32 sh[1032];   // 8 pad slots: the [+4] neighbor read of tid 255 chain 3 lands at byte 4096
  TV_INIT;
  for (u32 j = threadIdx.x; j < 1032; j += blockDim.x) sh[j] = j * 0x9E3779B9u ^ seedA;
  __syncthreads();
  u32 base = (u32)__cvta_generic_to_shared(&sh[threadIdx.x]);
  u32 s0 = base, s1 = base + 1024, s2 = base + 2048, s3 = base + 3072;
  u32 v0 = tv, v1 = ~tv, v2 = tv ^ 0x55555555u, v3 = tv ^ 0xAAAAAAAAu;
  for (u64 i = 0; i < iters; ++i) {
    R16(asm volatile("ld.shared.u32 %0,[%4+4];\n\txor.b32 %0,%0,%8;\n\tst.shared.u32 [%4],%0;\n\t"
                     "ld.shared.u32 %1,[%5+4];\n\txor.b32 %1,%1,%8;\n\tst.shared.u32 [%5],%1;\n\t"
                     "ld.shared.u32 %2,[%6+4];\n\txor.b32 %2,%2,%8;\n\tst.shared.u32 [%6],%2;\n\t"
                     "ld.shared.u32 %3,[%7+4];\n\txor.b32 %3,%3,%8;\n\tst.shared.u32 [%7],%3;"
                     : "+r"(v0), "+r"(v1), "+r"(v2), "+r"(v3)
                     : "r"(s0), "r"(s1), "r"(s2), "r"(s3), "r"(seedB));)
  }
  TAILGUARD(v0 + v1 + v2 + v3);
}

// M61 production-mix: t = lo32(A)*y + A; A = (t & M61) + (t >> 61).
// IMAD.WIDE + LOP3 + SHF + IADD3 on 64-bit — the M61 modmul skeleton, naturally toggling.
__global__ void k_mix61(u64 iters, u32 seedA, u32 seedB, u32* out) {
  TV_INIT;
  const u64 M61 = (1ull << 61) - 1;
  u64 A0 = tv | 1u, A1 = (u64)(tv ^ 0x33333333u) | 1u, A2 = ~(u64)tv | 1u, A3 = ((u64)tv << 3) | 1u;
  u32 y = (seedA & 0xFFFFFFFEu) | 0x9E3779B1u;
  for (u64 i = 0; i < iters; ++i) {
    #define M61STEP(A) { u64 t; u32 x; \
      asm volatile("cvt.u32.u64 %0,%1;" : "=r"(x) : "l"(A)); \
      asm volatile("mad.wide.u32 %0,%2,%3,%1;" : "=l"(t) : "l"(A), "r"(x), "r"(y)); \
      asm volatile("and.b64 %0,%1,%2;" : "=l"(A) : "l"(t), "l"(M61)); \
      asm volatile("shr.u64 %1,%1,61;\n\tadd.u64 %0,%0,%1;" : "+l"(A), "+l"(t)); }
    R16(M61STEP(A0) M61STEP(A1) M61STEP(A2) M61STEP(A3))
    #undef M61STEP
  }
  TAILGUARD((u32)(A0 + A1 + A2 + A3) + seedB);
}

// DRAM stream (memory-system power reference; lo/hi = buffer contents).
__global__ void k_dram_stream(u64 iters, u32 seed, u32* out, const uint4* buf, u32 n4) {
  u32 acc = seed;
  u32 stride = gridDim.x * blockDim.x;
  for (u64 i = 0; i < iters; ++i) {
    for (u32 j = blockIdx.x * blockDim.x + threadIdx.x; j < n4; j += stride) {
      uint4 v = __ldcv(&buf[j]);   // .cv: always hit DRAM/L2, no L1 reuse
      acc ^= v.x ^ v.y ^ v.z ^ v.w;
    }
  }
  TAILGUARD(acc);
}

// ---------------------------------------------------------------- host ------

typedef void (*KFn)(u64, u32, u32, u32*);

struct RunDesc {
  const char* name;
  KFn fn;             // null => dram
  u32 seedA, seedB;
  int opsPerIter;     // class-ops per thread per inner loop iteration
};

// opsPerIter = SASS class-ops per thread per source iteration (cuobjdump-verified):
// iadd3: ptxas fuses paired adds into 3-input IADD3 -> 32; imadwide_fb: AND+IMAD.WIDE -> 128;
// lds: LDS+STS+LOP3 per step -> 192; mix61: 5 core ops per M61 step -> 320 (+64 MOV overhead).
// imadwide_st is NOT run: ptxas strength-reduces loop-invariant x*y to a pure IADD.64 chain.
static RunDesc runs[] = {
  {"iadd3_lo",       k_iadd3,       0,           0,           32},
  {"iadd3_hi",       k_iadd3,       0xFFFFFFFFu, 0,           32},
  {"lop3_lo",        k_lop3,        0,           0,           64},
  {"lop3_hi",        k_lop3,        0xFFFFFFFFu, 0,           64},
  {"shf_lo",         k_shf,         0,           0,           64},
  {"shf_hi",         k_shf,         0,           0xFFFFFFFFu, 64},
  {"imad32_lo",      k_imad32,      0,           0,           64},
  {"imad32_hi",      k_imad32,      0x9E3779B1u, 0x7F4A7C15u, 64},
  {"imadwide_fb_lo", k_imadwide_fb, 0,           0x9E3779B0u, 128},
  {"imadwide_fb_hi", k_imadwide_fb, 0xFFFFFFFFu, 0x9E3779B0u, 128},
  {"add64_lo",       k_add64,       0,           0,           64},
  {"add64_hi",       k_add64,       0x55555555u, 0x55555555u, 64},
  {"ffma_lo",        k_ffma,        0x3F800000u, 0,           64},  // b=1.0 c=0.0
  {"ffma_hi",        k_ffma,        0xBF800000u, 0x3F800000u, 64},  // b=-1.0 c=1.0
  {"lds_lo",         k_lds,         0,           0,           192},
  {"lds_hi",         k_lds,         0,           0xFFFFFFFFu, 192},
  {"mix61",          k_mix61,       0,           0,           320},
  {"dram_lo",        nullptr,       0,           0,           0},
  {"dram_hi",        nullptr,       1,           0,           0},
};

struct Samples {
  std::vector<u32> mhz;
  std::vector<u32> mw;
  std::vector<u32> temp;
};

int main(int argc, char** argv) {
  if (argc < 2) { fprintf(stderr, "usage: %s <run|list|all|header> [seconds=75]\n", argv[0]); return 1; }
  std::string const what = argv[1];
  double seconds = argc > 2 ? atof(argv[2]) : 75.0;

  if (what == "list") { for (auto& r : runs) printf("%s\n", r.name); return 0; }
  if (what == "header") {
    printf("run,seconds,launches,mean_mhz,median_mhz,min_mhz,max_mhz,energy_J,watts_energy,watts_samples,mean_temp,gops,pj_per_op\n");
    return 0;
  }

  CUDA_CHECK(cudaSetDevice(0));
  int smCount = 0;
  CUDA_CHECK(cudaDeviceGetAttribute(&smCount, cudaDevAttrMultiProcessorCount, 0));

  NVML_CHECK(nvmlInit_v2());
  nvmlDevice_t dev;
  NVML_CHECK(nvmlDeviceGetHandleByIndex_v2(0, &dev));

  u32* dOut = nullptr;
  CUDA_CHECK(cudaMalloc(&dOut, 1 << 20));
  uint4* dBuf = nullptr;
  const u32 BUF_BYTES = 256u << 20;

  std::vector<std::string> names;
  if (what == "all") { for (auto& r : runs) names.push_back(r.name); }
  else names.push_back(what);

  for (auto& name : names) {
    RunDesc* rd = nullptr;
    for (auto& r : runs) if (name == r.name) { rd = &r; break; }
    if (!rd) { fprintf(stderr, "unknown run '%s'\n", name.c_str()); return 1; }
    bool isDram = rd->fn == nullptr;

    int block = 256;
    int grid;

    if (isDram) {
      if (!dBuf) CUDA_CHECK(cudaMalloc(&dBuf, BUF_BYTES));
      if (rd->seedA == 0) {
        CUDA_CHECK(cudaMemset(dBuf, 0, BUF_BYTES));
      } else {  // xorshift pattern: high bus toggle
        std::vector<u32> h(BUF_BYTES / 4);
        u32 s = 0x243F6A88u;
        for (auto& v : h) { s ^= s << 13; s ^= s >> 17; s ^= s << 5; v = s; }
        CUDA_CHECK(cudaMemcpy(dBuf, h.data(), BUF_BYTES, cudaMemcpyHostToDevice));
      }
      grid = smCount * 6;
    } else {
      int maxBlocks = 0;
      CUDA_CHECK(cudaOccupancyMaxActiveBlocksPerMultiprocessor(&maxBlocks, rd->fn, block, 0));
      grid = smCount * maxBlocks;
    }

    // Calibrate: aim ~250 ms per launch.
    u64 iters = isDram ? 4 : (1u << 18);
    double launchSec = 0;
    for (int tries = 0; tries < 12; ++tries) {
      auto t0 = std::chrono::steady_clock::now();
      if (isDram) k_dram_stream<<<grid, block>>>(iters, rd->seedA, dOut, dBuf, BUF_BYTES / 16);
      else rd->fn<<<grid, block>>>(iters, rd->seedA, rd->seedB, dOut);
      CUDA_CHECK(cudaDeviceSynchronize());
      launchSec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
      if (launchSec > 0.15 && launchSec < 0.4) break;
      double scale = 0.25 / std::max(launchSec, 1e-4);
      scale = std::min(std::max(scale, 0.05), 32.0);
      iters = std::max((u64)1, (u64)(iters * scale));
    }

    // Measure.
    std::atomic<bool> stop{false};
    Samples smp;
    std::thread sampler([&] {
      while (!stop.load(std::memory_order_relaxed)) {
        u32 mhz = 0, mw = 0, t = 0;
        if (nvmlDeviceGetClockInfo(dev, NVML_CLOCK_SM, &mhz) == NVML_SUCCESS) smp.mhz.push_back(mhz);
        if (nvmlDeviceGetPowerUsage(dev, &mw) == NVML_SUCCESS) smp.mw.push_back(mw);
        if (nvmlDeviceGetTemperature(dev, NVML_TEMPERATURE_GPU, &t) == NVML_SUCCESS) smp.temp.push_back(t);
        std::this_thread::sleep_for(std::chrono::milliseconds(50));
      }
    });

    unsigned long long e0 = 0, e1 = 0;
    bool haveEnergy = nvmlDeviceGetTotalEnergyConsumption(dev, &e0) == NVML_SUCCESS;

    auto t0 = std::chrono::steady_clock::now();
    u64 launches = 0;
    double elapsed = 0;
    if (isDram) k_dram_stream<<<grid, block>>>(iters, rd->seedA, dOut, dBuf, BUF_BYTES / 16);
    else rd->fn<<<grid, block>>>(iters, rd->seedA, rd->seedB, dOut);
    ++launches;
    while (true) {
      if (isDram) k_dram_stream<<<grid, block>>>(iters, rd->seedA, dOut, dBuf, BUF_BYTES / 16);
      else rd->fn<<<grid, block>>>(iters, rd->seedA, rd->seedB, dOut);
      ++launches;
      CUDA_CHECK(cudaDeviceSynchronize());
      elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
      if (elapsed >= seconds) break;
    }
    if (haveEnergy) haveEnergy = nvmlDeviceGetTotalEnergyConsumption(dev, &e1) == NVML_SUCCESS;
    stop.store(true);
    sampler.join();
    CUDA_CHECK(cudaGetLastError());

    auto mean = [](std::vector<u32>& v) -> double {
      if (v.empty()) return 0;
      double s = 0; for (u32 x : v) s += x; return s / v.size();
    };
    std::vector<u32> mz = smp.mhz;
    std::sort(mz.begin(), mz.end());
    u32 med = mz.empty() ? 0 : mz[mz.size() / 2];
    u32 mn = mz.empty() ? 0 : mz.front(), mx = mz.empty() ? 0 : mz.back();

    double energyJ = haveEnergy ? (e1 - e0) / 1000.0 : 0;
    double wattsE = haveEnergy ? energyJ / elapsed : 0;
    double wattsS = mean(smp.mw) / 1000.0;

    double totalOps, gops, pj;
    if (isDram) {
      totalOps = (double)launches * iters * BUF_BYTES;                 // bytes read
      gops = totalOps / elapsed / 1e9;                                 // GB/s
      pj = wattsE > 0 ? wattsE / (totalOps / elapsed) * 1e12 : 0;      // pJ/byte
    } else {
      totalOps = (double)launches * iters * rd->opsPerIter * (double)grid * block;
      gops = totalOps / elapsed / 1e9;
      pj = wattsE > 0 ? wattsE / (totalOps / elapsed) * 1e12 : 0;
    }

    printf("%s,%.1f,%llu,%.0f,%u,%u,%u,%.1f,%.1f,%.1f,%.0f,%.1f,%.4f\n",
           name.c_str(), elapsed, (unsigned long long)launches,
           mean(smp.mhz), med, mn, mx, energyJ, wattsE, wattsS,
           mean(smp.temp), gops, pj);
    fflush(stdout);
  }

  nvmlShutdown();
  return 0;
}
