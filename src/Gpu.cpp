// Copyright (C) Mihai Preda and George Woltman.

#include "Gpu.h"
#include "Proof.h"
#include "TimeInfo.h"
#include "Trig.h"
#include "state.h"
#include "Args.h"
#include "Signal.h"
#include "FFTConfig.h"
#include "Queue.h"
#include "Task.h"
#include "KernelCompiler.h"
#include "Saver.h"
#include "timeutil.h"
#include "TrigBufCache.h"
#include "fs.h"
#include "Sha3Hash.h"

#include <algorithm>
#include <limits>
#include <iomanip>
#include <array>
#include <cinttypes>
#include <cstring>

#define _USE_MATH_DEFINES
#include <cmath>
#include <ranges>
#include <utility>

#ifndef M_PIl
#define M_PIl 3.141592653589793238462643383279502884L
#endif

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884
#endif

#ifndef M_LN2l
#define M_LN2l 0.69314718055994530941723212145818L
#endif

#ifndef M_LN2
#define M_LN2 0.69314718055994530941723212145818
#endif

enum {
DEFAULT_CARRY_LEN = 8
};

namespace {

u32 carryLength(Args const& args, FFTConfig const& fft) {
  return fft.shape.fft_type == FFT31R2 ?
         args.value("RIESEL_CARRY_LEN", DEFAULT_CARRY_LEN) :
         DEFAULT_CARRY_LEN;
}

u32 kAt(u32 H, u32 line, u32 col) { return (line + col * H) * 2; }

double weight(u32 N, u64 E, u32 H, u32 line, u32 col, u32 rep) {
  return exp2l((long double)(extra(N, E, kAt(H, line, col) + rep)) / N);
}

double invWeight(u32 N, u64 E, u32 H, u32 line, u32 col, u32 rep) {
  return exp2l(-(long double)(extra(N, E, kAt(H, line, col) + rep)) / N);
}

// MSVC does not truly support long double.  Use expm1 rather than exp2 and subtracting one.
double weightM1(u32 N, u64 E, u32 H, u32 line, u32 col, u32 rep) {
  return expm1l(M_LN2l * (long double)(extra(N, E, kAt(H, line, col) + rep)) / N);
}

double invWeightM1(u32 N, u64 E, u32 H, u32 line, u32 col, u32 rep) {
  return expm1l(M_LN2l * - (long double)(extra(N, E, kAt(H, line, col) + rep)) / N);
}

constexpr array<u32, 2> RIESEL_Q_CANONICAL {2090860543u, 2141192191u};
constexpr array<u32, 2> RIESEL_THETA_CANONICAL {872791198u, 898129257u};
// Sub-2^30 fields permit Harvey's redundant [0,2q) Montgomery arithmetic.
// Lane order remains ascending for the simple two-prime Garner reconstruction.
constexpr array<u32, 2> RIESEL_Q_LAZY {1031798783u, 1038090239u};
constexpr array<u32, 2> RIESEL_THETA_LAZY {286716770u, 571191171u};
constexpr u32 M19_Q = (u32(1) << 19) - 1;
constexpr u32 M19_THETA = 256;

auto const& rieselQ(bool lazy) {
  return lazy ? RIESEL_Q_LAZY : RIESEL_Q_CANONICAL;
}

auto const& rieselTheta(bool lazy) {
  return lazy ? RIESEL_THETA_LAZY : RIESEL_THETA_CANONICAL;
}

u32 powMod32(u32 value, u64 exponent, u32 modulus) {
  u32 result = 1;
  while (exponent != 0) {
    if (exponent & 1) result = u32(u64(result) * value % modulus);
    exponent >>= 1;
    if (exponent != 0) value = u32(u64(value) * value % modulus);
  }
  return result;
}

u64 rieselMontPower(u32 exponent, bool inverse, bool lazy) {
  auto const& qValues = rieselQ(lazy);
  auto const& thetaValues = rieselTheta(lazy);
  u64 packed = 0;
  for (u32 lane = 0; lane != 2; ++lane) {
    u32 const q = qValues[lane];
    u32 theta = thetaValues[lane];
    if (inverse) theta = powMod32(theta, q - 2, q);
    u32 value = powMod32(theta, exponent, q);
    value = u32(u64(value) * ((u64(1) << 32) % q) % q);
    packed |= u64(value) << (lane * 32);
  }
  return packed;
}

u64 rieselMontScalar(u64 value, bool lazy) {
  auto const& qValues = rieselQ(lazy);
  u64 packed = 0;
  for (u32 lane = 0; lane != 2; ++lane) {
    u32 const q = qValues[lane];
    u32 normal = u32(value % q);
    normal = u32(u64(normal) * ((u64(1) << 32) % q) % q);
    packed |= u64(normal) << (lane * 32);
  }
  return packed;
}

// The Good--Thomas M19 experiment reuses the q1 storage/scheduling lane but
// stores that lane in the ordinary Mersenne representation.  Keep q0 packed
// exactly as before so all existing layouts and offsets remain unchanged.
u64 rieselFieldPower(u32 exponent, bool inverse, bool lazy, bool m19Field) {
  u64 packed = rieselMontPower(exponent, inverse, lazy);
  if (!m19Field) return packed;
  u32 theta = M19_THETA;
  if (inverse) theta = powMod32(theta, M19_Q - 2, M19_Q);
  u32 const value = powMod32(theta, exponent, M19_Q);
  return (packed & 0xffffffffULL) | (u64(value) << 32);
}

u64 rieselFieldScalar(u64 value, bool lazy, bool m19Field) {
  u64 packed = rieselMontScalar(value, lazy);
  if (!m19Field) return packed;
  return (packed & 0xffffffffULL) | (u64(value % M19_Q) << 32);
}

constexpr u64 GOLD_HOST_Q = 0xffffffff00000001ull;
constexpr u64 GOLD_HOST_THETA = 0x983d4271ce1d45bbull;

u64 goldHostMul(u64 a, u64 b) {
  return u64((unsigned __int128)a * b % GOLD_HOST_Q);
}

u64 goldHostPower(u64 value, u64 exponent) {
  u64 result = 1;
  while (exponent != 0) {
    if (exponent & 1) result = goldHostMul(result, value);
    exponent >>= 1;
    if (exponent != 0) value = goldHostMul(value, value);
  }
  return result;
}

u64 goldWeightPower(u32 exponent, bool inverse) {
  u64 theta = GOLD_HOST_THETA;
  if (inverse) theta = goldHostPower(theta, GOLD_HOST_Q - 2);
  return goldHostPower(theta, exponent);
}

double boundUnderOne(double x) { return std::min(x, nexttoward(1, 0)); }

float weight32(u32 N, u64 E, u32 H, u32 line, u32 col, u32 rep) {
  return float(exp2((double)(extra(N, E, kAt(H, line, col) + rep)) / N));
}

float invWeight32(u32 N, u64 E, u32 H, u32 line, u32 col, u32 rep) {
  return float(exp2(-(double)(extra(N, E, kAt(H, line, col) + rep)) / N));
}

float weightM132(u32 N, u64 E, u32 H, u32 line, u32 col, u32 rep) {
  return float(expm1(M_LN2 * (double)(extra(N, E, kAt(H, line, col) + rep)) / N));
}

float invWeightM132(u32 N, u64 E, u32 H, u32 line, u32 col, u32 rep) {
  return float(expm1(M_LN2 * - (double)(extra(N, E, kAt(H, line, col) + rep)) / N));
}

Weights genWeights(FFTConfig fft, u64 E, u32 W, u32 H, u32 nW,
                   bool nvidiaGpu, bool lazyRiesel, bool goldPair,
                   bool m19Field) {
  u32 const N = 2u * W * H;
  u32 const groupWidth = W / nW;

  vector<double> weightsConstIF;
  vector<double> weightsIF;

  if (fft.FFT_FP64) {
    // Inverse + Forward
    for (u32 thread = 0; thread < groupWidth; ++thread) {
      auto iw = invWeight(N, E, H, 0, thread, 0);
      auto w = weight(N, E, H, 0, thread, 0);
      weightsIF.push_back(2 * boundUnderOne(iw));
      weightsIF.push_back(2 * w);
    }

    // the group order matches CarryA/M (not fftP/CarryFused).
    for (u32 gy = 0; gy < H; ++gy) {
      weightsIF.push_back(invWeightM1(N, E, H, gy, 0, 0));
      weightsIF.push_back(weightM1(N, E, H, gy, 0, 0));
    }

    // nVidia GPUs have a fast constant cache that only works on buffer sizes less than 64KB.  Create two smaller buffers
    // that can be used to create the large group order buffer with a single multiply.
    if (nvidiaGpu) {
      for (u32 gy = 0; gy < 64; ++gy) {
        weightsConstIF.push_back(invWeightM1(N, E, H, gy, 0, 0));
        weightsConstIF.push_back(weightM1(N, E, H, gy, 0, 0));
      }
      for (u32 gy = 0; gy < H; gy += 64) {
        weightsConstIF.push_back(invWeightM1(N, E, H, gy, 0, 0));
        weightsConstIF.push_back(weightM1(N, E, H, gy, 0, 0));
      }
    }
  }

  else if (fft.FFT_FP32) {
    vector<float> weightsConstIF32;
    vector<float> weightsIF32;
    // Inverse + Forward
    for (u32 thread = 0; thread < groupWidth; ++thread) {
      auto iw = invWeight32(N, E, H, 0, thread, 0) ;
      auto w = weight32(N, E, H, 0, thread, 0) ;
      // Weights are scaled by 2^-24 and 2^48 so that multiplicaton by 1/epsilon does not generate infinty results (width and height variant 2).
      iw = iw * 281474976710656.0f;
      w = w * 0.000000059604644775390625f;
      weightsIF32.push_back(iw);
      weightsIF32.push_back(w);
    }

    // the group order matches CarryA/M (not fftP/CarryFused).
    for (u32 gy = 0; gy < H; ++gy) {
      weightsIF32.push_back(invWeightM132(N, E, H, gy, 0, 0));
      weightsIF32.push_back(weightM132(N, E, H, gy, 0, 0));
    }

    // nVidia GPUs have a fast constant cache that only works on buffer sizes less than 64KB.  Create two smaller buffers
    // that can be used to create the large group order buffer with a single multiply.
    if (nvidiaGpu) {
      for (u32 gy = 0; gy < 64; ++gy) {
        weightsConstIF32.push_back(invWeightM132(N, E, H, gy, 0, 0));
        weightsConstIF32.push_back(weightM132(N, E, H, gy, 0, 0));
      }
      for (u32 gy = 0; gy < H; gy += 64) {
        weightsConstIF32.push_back(invWeightM132(N, E, H, gy, 0, 0));
        weightsConstIF32.push_back(weightM132(N, E, H, gy, 0, 0));
      }
    }

    // Copy the float vectors to the double vectors
    weightsIF.resize(weightsIF32.size() / 2);
    memcpy((double *) weightsIF.data(), weightsIF32.data(), weightsIF32.size() * sizeof(float));
    weightsConstIF.resize(weightsConstIF32.size() / 2);
    memcpy((double *) weightsConstIF.data(), weightsConstIF32.data(), weightsConstIF32.size() * sizeof(float));
  }

  else if (fft.shape.fft_type == FFT3161 && goldPair) {
    vector<u64> goldWeights;
    goldWeights.reserve(2 * (W + H));
    auto append = [&](u32 exponent) {
      goldWeights.push_back(goldWeightPower(exponent, true));
      goldWeights.push_back(goldWeightPower(exponent, false));
    };
    for (u32 x = 0; x != W; ++x) append(extra(N, E, 2 * x * H));
    for (u32 line = 0; line != H; ++line) append(extra(N, E, 2 * line));
    weightsIF.resize(goldWeights.size());
    memcpy(weightsIF.data(), goldWeights.data(),
           goldWeights.size() * sizeof(u64));
  }

  else if (fft.shape.fft_type == FFT31R2) {
    if (fft.NTT_GF61) {
      // q0/q1 are vectorized into the low/high halves of each u64.  They are
      // still reduced independently; packing only changes scheduling and the
      // in-register/global-memory representation.
      vector<u64> packedWeights;
      packedWeights.reserve(2 * (W + H));
      auto append = [&](u32 exponent) {
        packedWeights.push_back(rieselFieldPower(exponent, true, lazyRiesel,
                                                 m19Field));
        packedWeights.push_back(rieselFieldPower(exponent, false, lazyRiesel,
                                                 m19Field));
      };
      for (u32 x = 0; x != W; ++x) append(extra(N, E, 2 * x * H));
      for (u32 line = 0; line != H; ++line) append(extra(N, E, 2 * line));
      weightsIF.resize(packedWeights.size());
      memcpy(weightsIF.data(), packedWeights.data(), packedWeights.size() * sizeof(u64));
    } else {
      // Independent representation: one uint plane per Riesel field.
      vector<u32> fieldWeights;
      fieldWeights.reserve(4 * (W + H));
      for (u32 lane = 0; lane != 2; ++lane) {
        auto append = [&](u32 exponent) {
          fieldWeights.push_back(u32(rieselFieldPower(exponent, true, lazyRiesel,
                                                      m19Field) >> (32 * lane)));
          fieldWeights.push_back(u32(rieselFieldPower(exponent, false, lazyRiesel,
                                                      m19Field) >> (32 * lane)));
        };
        for (u32 x = 0; x != W; ++x) append(extra(N, E, 2 * x * H));
        for (u32 line = 0; line != H; ++line) append(extra(N, E, 2 * line));
      }
      weightsIF.resize((fieldWeights.size() + 1) / 2);
      memcpy(weightsIF.data(), fieldWeights.data(), fieldWeights.size() * sizeof(u32));
    }
  }

  return Weights{.weightsConstIF=weightsConstIF, .weightsIF=weightsIF};
}

string toLiteral(i32 value) { return to_string(value); }
string toLiteral(u32 value) { return to_string(value) + 'u'; }
[[maybe_unused]] string toLiteral(long value) { return to_string(value) + "ll"; }
[[maybe_unused]] string toLiteral(unsigned long value) { return to_string(value) + "ull"; }
[[maybe_unused]] string toLiteral(long long value) { return to_string(value) + "ll"; }
[[maybe_unused]] string toLiteral(unsigned long long value) { return to_string(value) + "ull"; }

template<typename F>
string toLiteral(F value) {
  std::ostringstream ss;
  ss << std::setprecision(numeric_limits<F>::max_digits10) << value;
  if (sizeof(F) == 4) ss << "f";
  string s = std::move(ss).str();

  // verify exact roundtrip
  [[maybe_unused]] F back = 0;
  sscanf(s.c_str(), (sizeof(F) == 4) ? "%f" : "%lf", &back);
  assert(back == value);
  
  return s;
}

template<typename T>
string toLiteral(const vector<T>& v) {
  assert(!v.empty());
  string s = "{";
  for (auto x : v) {
    s += toLiteral(x) + ",";
  }
  s += "}";
  return s;
}

template<typename T, size_t N>
string toLiteral(const std::array<T, N>& v) {
  string s = "{";
  for (T x : v) {
    s += toLiteral(x) + ",";
  }
  s += "}";
  return s;
}

string toLiteral(const string& s) { return s; }

[[maybe_unused]] string toLiteral(float2 cs) { return "U2("s + toLiteral(cs.first) + ',' + toLiteral(cs.second) + ')'; }
[[maybe_unused]] string toLiteral(double2 cs) { return "U2("s + toLiteral(cs.first) + ',' + toLiteral(cs.second) + ')'; }
[[maybe_unused]] string toLiteral(int2 cs) { return "U2("s + toLiteral(cs.first) + ',' + toLiteral(cs.second) + ')'; }
[[maybe_unused]] string toLiteral(uint2 cs) { return "U2("s + toLiteral(cs.first) + ',' + toLiteral(cs.second) + ')'; }
[[maybe_unused]] string toLiteral(ulong2 cs) { return "U2("s + toLiteral(cs.first) + ',' + toLiteral(cs.second) + ')'; }

template<typename T>
string toDefine(const string& k, const T& v) { return " -D"s + k + '=' + toLiteral(v); }

template<typename T>
string toDefine(const T& vect) {
  string s;
  for (const auto& [k, v] : vect) { s += toDefine(k, v); }
  return s;
}

constexpr bool isInList(const string& s, initializer_list<string> list) {
  for (const string& e : list) { if (e == s) { return true; }}
  return false;
}

string clDefines(Args& args, cl_device_id id, FFTConfig fft, const vector<KeyVal>& extraConf, u64 E, bool doLog,
                 bool &tail_single_wide, bool &tail_single_kernel, u32 &in_place, u32 &pad_size, u32 &wmul) {
  map<string, string> config;

  // Highest priority is the requested "extra" conf
  config.insert(extraConf.begin(), extraConf.end());

  // Next, args config
  config.insert(args.flags.begin(), args.flags.end());

  // Lowest priority: the per-FFT config if any
  if (auto it = args.perFftConfig.find(fft.shape.spec()); it != args.perFftConfig.end()) {
    // log("Found %s\n", fft.shape.spec().c_str());
    config.insert(it->second.begin(), it->second.end());
  }

  // Default value for -use options that must also be parsed in C++ code
  tail_single_wide = false, tail_single_kernel = true;         // Default tailSquare is double-wide in one kernel
  in_place = 0;                                         // Default is not in-place
  wmul = 2;                                             // Default is carryFused processes two lines at a time
  pad_size = isAmdGpu(id) ? 256 : 0;                    // Default is 256 bytes for AMD, 0 for others

  // Validate -use options
  for (const auto& [k, v] : config) {
    bool const isValid = isInList(k, {
                              "FAST_BARRIER",
                              "STATS",
                              "SPIN_STATS",
                              "ROE_COUNT",
                              "PARITY_CORRECT",
                              "PARITY_PREPARED",        // 1: overlap; 2: direct scatter; 3: physical natural parity
                              "PARITY_PACKED",          // Ballot-packed natural/prepared parity planes (CUDA only)
                              "PARITY_LAZY",            // Load expected parity only for residual-risk candidates
                              "FOLD_SYNDROME",          // Sparse-error folded M31 sidecar; 2 enables diagnostic validation
                              "FOLD_TRANSFORM",         // Run the 512K cyclic M31 square on the folded sidecar
                              "FOLD_FACTOR",            // 4: lifted 1M fold; 8: parity-free 512K; 16: parity-assisted 256K
                              "FOLD_SHAPE",             // 0: 256:4:256, 1: 256:2:512, 2: 512:2:256
                              "FOLD_ZERO_WIDTH",        // Timing lower bound: omit both side width edges after priming
                              "QUOTIENT_STATS",
                              "FP32_WIDE_QUOTIENT",      // Sign-magnitude rounding for the full FP32/M61 quotient range
                              "FP32_CMUL64",
                              "FP32_CFMA64",
                              "COEFF_RANGE_STATS",
                              "IN_SIZEX",
                              "IN_WG",
                              "OUT_SIZEX",
                              "OUT_WG",
                              "UNROLL_H",
                              "UNROLL_W",
                              "WIDTH_NW",              // Experimental per-thread width ownership (2, 4, or 8)
                              "ZEROHACK_H",
                              "ZEROHACK_W",
                              "NO_ASM",
                              "DEBUG",
                              "CARRY64",
                              "BIGLIT",                 // Deprecated
                              "NONTEMPORAL",            // Deprecated
                              "INPLACE",
                              "PAD",
                              "MIDDLE_IN_LDS_TRANSPOSE",
                              "MIDDLE_OUT_LDS_TRANSPOSE",
                              "TAIL_KERNELS",
                              "TAIL_TRIGS",
                              "TAIL_TRIGS31",
                              "TAIL_TRIGS32",
                              "TAIL_TRIGS61",
                              "TABMUL_CHAIN",
                              "TABMUL_CHAIN31",
                              "TABMUL_CHAIN32",
                              "TABMUL_CHAIN61",
                              "ENABLE_BETTER_ONEPAIRSQ", // Alternate M61 pair-square identity in tailsquare.cl
                              "ENABLE_PARALLEL_ONEPAIRSQ", // Nine-product pair square with independent square chains
                              "M61_CMUL4",               // Four independent products instead of Karatsuba GF61 multiply
                              "FULL_MIDDLE_ROOTS61",     // Two coalesced 32-MiB M61 middleMul2 root planes
                              "DETACHED_M31_EDGE",       // FFT3161: exact external M31 width kernels; FFT323161: timing-only scaffold
                              "MODM31",
                              "LOADS","STORES",
                              "NOREG",                  // CUDA - experimental
                              "WMUL",
                              "MULTI_Q",
                              "RIESEL_FUSED",           // Experimental fused carry for the three-prime M31R2 path
                              "RIESEL_PACKED",          // Vectorize the two Riesel fields into one transform stream
                              "RIESEL_FOUR",            // Add M61 to M31R2 for an experimental 154-bit, 2M path
                              "GOOD_THOMAS3",           // 3M NTT fields as three packed 1M Good-Thomas channels
                              "GOOD_THOMAS7",           // 3.5M fields as seven packed 512K Good-Thomas channels
                              "GOOD_THOMAS9",           // 4.5M FP32+M61 as nine packed 512K channels
                              "RIESEL_CARRY_LEN",       // Experimental long-carry block length for M31R2
                              "RIESEL_PTX_MONT",        // Experimental inline-PTX Montgomery multiply for Riesel fields
                              "RIESEL_LAZY",            // Use sub-2^30 Riesel fields and Harvey redundant ranges
                              "M19_FIELD",              // Replace q1 with M19 in the 3M Good-Thomas experiment
                              "GOLD_PAIR",              // Experimental full scalar Goldilocks plane in even/odd pairs
                              "MIDCARRY_FUSED",         // Two-phase middle/width/carry fusion for CUDA FFT3161
                              "GRAPHS",
                              "L1CUDA",
                              "L2_PERSIST",             // CUDA persisting-L2 window for the largest trig table
                              "ASYNC_MID61",            // cp.async tile-looped fftMiddleOutGF61: 1=full staging, 2=half
                              "ASYNC_TILES",            // Tiles per group for ASYNC_MID61 (default 4)
                              "FUSED31",                // Fused middle-in + tail for the M31 plane (pair-resident tile)
                              "CARRY_NOWAIT",           // Timing scaffold: skip the carry shuttle (WRONG results)
                              "CARRY_EARLY",            // Hoist carryFused shuttle wait+load before weights/shuffle
                              "CARRY_ACQREL"            // Release/acquire shuttle handshake instead of device fences
                            });
    if (!isValid) {
      log("Warning: unrecognized -use key '%s'\n", k.c_str());
    }

    // Some -use options are needed in both OpenCL code and C++ initialization code
    if (k == "TAIL_KERNELS") {
      if (atoi(v.c_str()) == 0) tail_single_wide = true, tail_single_kernel = true;
      if (atoi(v.c_str()) == 1) tail_single_wide = true, tail_single_kernel = false;
      if (atoi(v.c_str()) == 2) tail_single_wide = false, tail_single_kernel = true;
      if (atoi(v.c_str()) == 3) tail_single_wide = false, tail_single_kernel = false;
    }
    if (k == "INPLACE") in_place = atoi(v.c_str());
    if (k == "WMUL") wmul = atoi(v.c_str());
    if (k == "PAD") pad_size = atoi(v.c_str());
  }

  // The padded shuffle layouts only have specialized address formulas for
  // radix 4 and radix 8.  Radix 2 uses the compact, generic shuffle layout.
  if (config.contains("WIDTH_NW") &&
      stoul(config.at("WIDTH_NW")) == 2) {
    config["LDSPAD_W"] = "0";
  }

  // Maximum WMUL is 32KB / (WIDTH * SHUFL_BYTES_W).  If using the 32KB maximum, LDS padding must be disabled.
  // Furthermore, I've seen the CUDA compiler refuse to create a kernel with 1024 threads.  Thus, we limit WMUL to 2 for a 1K width and to 1 for a 4K width.
  {
    u32 const shufl_bytes_w = args.value("SHUFL_BYTES_W", 8);
    u32 max_wmul = 32768 / (fft.shape.width * shufl_bytes_w);
    if (max_wmul > 2 && fft.shape.width >= 1024) max_wmul = 2;
    if (max_wmul > 1 && fft.shape.width >= 4096) max_wmul = 1;
    if (wmul > max_wmul) {
      wmul = max_wmul;
      config["WMUL"] = to_string(wmul);
      log("WMUL setting too large for this FFT width.  Changing to WMUL=%d\n", wmul);
    }
    if (fft.shape.width * shufl_bytes_w * wmul >= 32768) {
      log("Local shared memory limit of 32KB exceeded.  Changing to LDSPAD_W=0\n");
      config["LDSPAD_W"] = to_string(0);
    }
  }

  // L2_STRIPING is not allowed if INPLACE=0.  Maximum L2_STRIPING is WIDTH/64 if MULTI_Q=0 and WIDTH/128 if MULTI_Q=1.
  // Technically, L2_STRIPING of WIDTH/32, MULTI_Q=0 could be allowed but that is just a more complicated way to implement L2_STRIPING=0.
  // Also, WIDTH/64, MULTI_Q=1 could be allowed with some marker/sync code changes but that is very similar to L2_STRIPING=0.
  {
    u32 l2_striping = args.value("L2_STRIPING", 0);
    u32 multi_q = args.value("MULTI_Q", 0);
    if (l2_striping && !in_place) {
      config["L2_STRIPING"] = to_string(0);
      args.flags["L2_STRIPING"] = to_string(0);
      log("L2_STRIPING is only allowed if INPLACE=1.  Changing to L2_STRIPING=0.\n");
    }
    else if (multi_q == 0 && l2_striping > fft.shape.width/64) {
      config["L2_STRIPING"] = to_string(fft.shape.width/64);
      args.flags["L2_STRIPING"] = to_string(fft.shape.width/64);
      log("Max L2_STRIPING when MULTI_Q=0 exceeded.  Changing to L2_STRIPING=%u.\n", fft.shape.width/64);
    }
    else if (multi_q > 0 && l2_striping > fft.shape.width/128) {
      config["L2_STRIPING"] = to_string(fft.shape.width/128);
      args.flags["L2_STRIPING"] = to_string(fft.shape.width/128);
      log("Max L2_STRIPING when MULTI_Q=1 exceeded.  Changing to L2_STRIPING=%u.\n", fft.shape.width/128);
    }
  }

  string defines = toDefine(config);
  if (doLog) { log("config: %s\n", defines.c_str()); }

  defines += toDefine("EXP", E);
  u32 const widthNW = config.contains("WIDTH_NW") ?
    u32(stoul(config.at("WIDTH_NW"))) : fft.shape.nW();
  if ((widthNW != 2 && widthNW != 4 && widthNW != 8) ||
      fft.shape.width % widthNW != 0) {
    throw runtime_error("WIDTH_NW must be 2, 4, or 8 and divide the FFT width");
  }

  defines += toDefine(initializer_list<pair<string, u32>>{
                    {"WIDTH", fft.shape.width},
                    {"SMALL_HEIGHT", fft.shape.height},
                    {"MIDDLE", fft.shape.middle},
                    {"CARRY_LEN", carryLength(args, fft)},
                    {"NW", widthNW},
                    {"NH", fft.shape.nH()}
                  });

  if (isAmdGpu(id)) { defines += toDefine("AMDGPU", 1); }
  if (isNvidiaGpu(id)) { defines += toDefine("NVIDIAGPU", 1); }
  if (isNvidiaGpu(id)) { defines += toDefine("CC", getNvidiaComputeCapability(id)); }

  if ((fft.carry == CARRY_AUTO && fft.shape.needsLargeCarry(E)) || (fft.carry == CARRY_64)) {
    if (doLog) { log("Using CARRY64\n"); }
    defines += toDefine("CARRY64", 1);
  }

  u32 const N = fft.shape.size();
  defines += toDefine("FFT_VARIANT", fft.variant);
  bool const rieselPair = fft.shape.fft_type == FFT31R2 && fft.NTT_GF61 &&
                          !fft.NTT_RIESEL;
  bool const rieselFour = fft.shape.fft_type == FFT31R2 && fft.NTT_GF61 &&
                          fft.NTT_RIESEL;
  bool const rieselLazy = fft.shape.fft_type == FFT31R2 &&
                          args.value("RIESEL_LAZY", 0);
  bool const goldPair = fft.shape.fft_type == FFT3161 &&
                        args.value("GOLD_PAIR", 0);
  bool const goodThomas3 =
    (fft.shape.fft_type == FFT31R2 || fft.shape.fft_type == FFT323161) &&
    args.value("GOOD_THOMAS3", 0);
  bool const goodThomas7 = fft.shape.fft_type == FFT31R2 &&
                           args.value("GOOD_THOMAS7", 0);
  bool const goodThomas9 = fft.shape.fft_type == FFT3261 &&
                           args.value("GOOD_THOMAS9", 0);
  bool const m19Field = fft.shape.fft_type == FFT31R2 &&
                        (goodThomas3 || goodThomas7) &&
                        args.value("M19_FIELD", 0);
  defines += toDefine("RIESEL_PAIR", (int)rieselPair);
  defines += toDefine("RIESEL_FOUR", (int)rieselFour);
  defines += toDefine("RIESEL_LAZY", (int)rieselLazy);
  defines += toDefine("GOLD_PAIR", (int)goldPair);
  defines += toDefine("GOOD_THOMAS3", (int)goodThomas3);
  defines += toDefine("GOOD_THOMAS7", (int)goodThomas7);
  defines += toDefine("GOOD_THOMAS9", (int)goodThomas9);
  defines += toDefine("M19_FIELD", (int)m19Field);
  defines += toDefine("MAXBPW", (u32)(fft.maxBpw() * 100.0f));

  if (fft.FFT_FP64 || fft.FFT_FP32) {
    defines += toDefine("WEIGHT_STEP", weightM1(N, E, fft.shape.height * fft.shape.middle, 0, 0, 1));
    defines += toDefine("IWEIGHT_STEP", invWeightM1(N, E, fft.shape.height * fft.shape.middle, 0, 0, 1));
    if (fft.FFT_FP64) defines += toDefine("TAILT", root1Fancy(fft.shape.height * 2, 1));
    else defines += toDefine("TAILT", root1FancyFP32(fft.shape.height * 2, 1));

    TrigCoefs const coefs = trigCoefs(fft.shape.size() / 4);
    defines += toDefine("TRIG_SCALE", int(coefs.scale));
    defines += toDefine("TRIG_SIN",  coefs.sinCoefs);
    defines += toDefine("TRIG_COS",  coefs.cosCoefs);
  }
  if (fft.NTT_GF31) {
    defines += toDefine("TAILTGF31", root1GF31(fft.shape.height * 2, 1));
  }
  if (fft.NTT_GF61) {
    defines += toDefine("TAILTGF61", goldPair ?
                        root1GoldPair(fft.shape.height * 2, 1) :
                        (rieselPair ?
                         root1Riesel(fft.shape.height * 2, 1, rieselLazy) :
                         root1GF61(fft.shape.height * 2, 1)));
  }
  if (goldPair) {
    u32 const deltaOne = step(N, E);
    u32 const deltaX = u32(u64(deltaOne) *
      (2u * (fft.shape.width / fft.shape.nW()) *
       fft.shape.height * fft.shape.middle) % N);
    u64 const inverseNd = goldHostPower(N / 2, GOLD_HOST_Q - 2);
    assert(goldHostPower(GOLD_HOST_THETA, N) == 2);
    defines += toDefine("GOLD_DELTA_ONE", deltaOne);
    defines += toDefine("GOLD_DELTA_X", deltaX);
    defines += toDefine("GOLD_FWD_ONE", goldWeightPower(deltaOne, false));
    defines += toDefine("GOLD_INV_ONE", goldWeightPower(deltaOne, true));
    defines += toDefine("GOLD_FWD_X", goldWeightPower(deltaX, false));
    defines += toDefine("GOLD_INV_X", goldWeightPower(deltaX, true));
    defines += toDefine("GOLD_INV_ND", inverseNd);
  }
  if (fft.shape.fft_type == FFT31R2) {
    auto qValues = rieselQ(rieselLazy);
    if (m19Field) qValues[1] = M19_Q;
    u32 const deltaOne = step(N, E);
    u32 const deltaX = u32(u64(step(N, E)) * (2u * (fft.shape.width / fft.shape.nW()) *
                                               fft.shape.height * fft.shape.middle) % N);
    // The unnormalised inverse radix-3 leaves a factor of three for the joint
    // CRT to divide exactly.  Each field therefore removes only the 2B scale
    // of its packed B-word power-of-two DGT here.
    u64 const scale = u64(2) * (goodThomas3 ? N / 3 : N);
    u64 const scaleInverse0 = powMod32(u32(scale % qValues[0]), qValues[0] - 2, qValues[0]);
    u64 const scaleInverse1 = powMod32(u32(scale % qValues[1]), qValues[1] - 2, qValues[1]);
    ulong2 tail = root1Riesel(fft.shape.height * 2, 1, rieselLazy);
    if (m19Field) {
      uint2 const tail19 = root1M19(fft.shape.height * 2, 1);
      tail.first = (tail.first & 0xffffffffULL) | (u64(tail19.first) << 32);
      tail.second = (tail.second & 0xffffffffULL) | (u64(tail19.second) << 32);
    }
    auto lane = [](u64 packed, u32 index) { return u32(packed >> (32 * index)); };
    u64 const fwdOne = rieselFieldPower(deltaOne, false, rieselLazy, m19Field);
    u64 const invOne = rieselFieldPower(deltaOne, true, rieselLazy, m19Field);
    u64 const fwdX = rieselFieldPower(deltaX, false, rieselLazy, m19Field);
    u64 const invX = rieselFieldPower(deltaX, true, rieselLazy, m19Field);
    u64 const montOne = rieselFieldPower(0, false, rieselLazy, m19Field);
    u64 const invScale = (rieselFieldScalar(scaleInverse0, rieselLazy, m19Field) & 0xffffffffULL) |
                         (rieselFieldScalar(scaleInverse1, rieselLazy, m19Field) & 0xffffffff00000000ULL);
    defines += toDefine("RIESEL_DELTA_ONE", deltaOne);
    defines += toDefine("RIESEL_DELTA_X", deltaX);
    defines += toDefine("RIESEL_WEIGHT_STRIDE", 2u * (fft.shape.width + fft.shape.height * fft.shape.middle));
    defines += toDefine("RIESEL_FWD_ONE", fwdOne);
    defines += toDefine("RIESEL_INV_ONE", invOne);
    defines += toDefine("RIESEL_FWD_X", fwdX);
    defines += toDefine("RIESEL_INV_X", invX);
    defines += toDefine("RIESEL_MONT_ONE", montOne);
    defines += toDefine("RIESEL_INV_SCALE", invScale);
    for (u32 i = 0; i != 2; ++i) {
      string const prefix = "RIESEL" + to_string(i) + '_';
      defines += toDefine(prefix + "FWD_ONE", lane(fwdOne, i));
      defines += toDefine(prefix + "INV_ONE", lane(invOne, i));
      defines += toDefine(prefix + "FWD_X", lane(fwdX, i));
      defines += toDefine(prefix + "INV_X", lane(invX, i));
      defines += toDefine(prefix + "MONT_ONE", lane(montOne, i));
      defines += toDefine(prefix + "INV_SCALE", lane(invScale, i));
    }
    defines += toDefine("TAILTR0", uint2{lane(tail.first, 0), lane(tail.second, 0)});
    defines += toDefine("TAILTR1", uint2{lane(tail.first, 1), lane(tail.second, 1)});
  }

  // Send the FFT/NTT type and booleans that enable/disable code for each possible FP and NTT
  defines += toDefine("FFT_TYPE", (int) fft.shape.fft_type);
  defines += toDefine("WordSize", fft.WordSize);

  // When using multiple NTT primes or hybrid FFT/NTT, each FFT/NTT prime's data buffer and trig values are combined into one buffer.
  // The openCL code needs to know the offset to the data and trig values.  Distances are in "number of double2 values".
  bool const widthTrigUsesCombo = fft.shape.width == fft.shape.height &&
                                  widthNW == fft.shape.nH();
  u32 const wTrigMiddle = widthTrigUsesCombo ? fft.shape.middle : 0;
  u32 const wTrigHeight = widthTrigUsesCombo ? fft.shape.height : 0;
  u32 const wTrigNH = widthTrigUsesCombo ? fft.shape.nH() : 0;
  if (fft.NTT_RIESEL) {
    u32 const data31Plane = GF31_DATA_SIZE(fft.shape.width, fft.shape.middle, fft.shape.height, in_place, pad_size) / 2;
    u32 const data61Plane = fft.NTT_GF61 ?
      GF61_DATA_SIZE(fft.shape.width, fft.shape.middle, fft.shape.height, in_place, pad_size) / 2 : 0;
    u32 const w31Plane = SMALLTRIG_GF31_DIST(fft.shape.width, wTrigMiddle, wTrigHeight, wTrigNH);
    u32 const w61Plane = fft.NTT_GF61 ?
      SMALLTRIG_GF61_DIST(fft.shape.width, wTrigMiddle, wTrigHeight, wTrigNH) : 0;
    u32 const m31Plane = MIDDLETRIG_GF31_DIST(fft.shape.width, fft.shape.middle, fft.shape.height);
    u32 const m61Plane = fft.NTT_GF61 ?
      MIDDLETRIG_GF61_DIST(fft.shape.width, fft.shape.middle, fft.shape.height) : 0;
    u32 const h31Plane = SMALLTRIGCOMBO_GF31_DIST(fft.shape.width, fft.shape.middle, fft.shape.height, fft.shape.nH());
    u32 const h61Plane = fft.NTT_GF61 ?
      SMALLTRIGCOMBO_GF61_DIST(fft.shape.width, fft.shape.middle, fft.shape.height, fft.shape.nH()) : 0;
    defines += toDefine("DISTGF31", 0);
    defines += toDefine("DISTWTRIGGF31", 0);
    defines += toDefine("DISTMTRIGGF31", 0);
    defines += toDefine("DISTHTRIGGF31", 0);
    defines += toDefine("DISTGF61", data31Plane);
    defines += toDefine("DISTWTRIGGF61", w31Plane);
    defines += toDefine("DISTMTRIGGF61", m31Plane);
    defines += toDefine("DISTHTRIGGF61", h31Plane);
    defines += toDefine("DISTR0", data31Plane + data61Plane);
    defines += toDefine("DISTR1", 2 * data31Plane + data61Plane);
    defines += toDefine("DISTWTR0", w31Plane + w61Plane);
    defines += toDefine("DISTWTR1", 2 * w31Plane + w61Plane);
    defines += toDefine("DISTMTR0", m31Plane + m61Plane);
    defines += toDefine("DISTMTR1", 2 * m31Plane + m61Plane);
    defines += toDefine("DISTHTR0", h31Plane + h61Plane);
    defines += toDefine("DISTHTR1", 2 * h31Plane + h61Plane);
  }
  else if (fft.FFT_FP64 && fft.NTT_GF31) {
    // GF31 data is located after the FP64 data.  Compute size of the FP64 data and trigs.
    defines += toDefine("DISTGF31",      FP64_DATA_SIZE(fft.shape.width, fft.shape.middle, fft.shape.height, in_place, pad_size) / 2);
    defines += toDefine("DISTWTRIGGF31", SMALLTRIG_FP64_DIST(fft.shape.width, wTrigMiddle, wTrigHeight, wTrigNH));
    defines += toDefine("DISTMTRIGGF31", MIDDLETRIG_FP64_DIST(fft.shape.width, fft.shape.middle, fft.shape.height));
    defines += toDefine("DISTHTRIGGF31", SMALLTRIGCOMBO_FP64_DIST(fft.shape.width, fft.shape.middle, fft.shape.height, fft.shape.nH()));
  }
  else if (fft.FFT_FP32 && fft.NTT_GF31 && fft.NTT_GF61) {
    // GF31 and GF61 data is located after the FP32 data.  Compute size of the FP32 data and trigs.
    u32 sz1, sz2, sz3, sz4;
    defines += toDefine("DISTGF31",      sz1 = FP32_DATA_SIZE(fft.shape.width, fft.shape.middle, fft.shape.height, in_place, pad_size) / 2);
    defines += toDefine("DISTWTRIGGF31", sz2 = SMALLTRIG_FP32_DIST(fft.shape.width, wTrigMiddle, wTrigHeight, wTrigNH));
    defines += toDefine("DISTMTRIGGF31", sz3 = MIDDLETRIG_FP32_DIST(fft.shape.width, fft.shape.middle, fft.shape.height));
    defines += toDefine("DISTHTRIGGF31", sz4 = SMALLTRIGCOMBO_FP32_DIST(fft.shape.width, fft.shape.middle, fft.shape.height, fft.shape.nH()));
    defines += toDefine("DISTGF61",      sz1 + GF31_DATA_SIZE(fft.shape.width, fft.shape.middle, fft.shape.height, in_place, pad_size) / 2);
    defines += toDefine("DISTWTRIGGF61", sz2 + SMALLTRIG_GF31_DIST(fft.shape.width, wTrigMiddle, wTrigHeight, wTrigNH));
    defines += toDefine("DISTMTRIGGF61", sz3 + MIDDLETRIG_GF31_DIST(fft.shape.width, fft.shape.middle, fft.shape.height));
    defines += toDefine("DISTHTRIGGF61", sz4 + SMALLTRIGCOMBO_GF31_DIST(fft.shape.width, fft.shape.middle, fft.shape.height, fft.shape.nH()));
  }
  else if (fft.FFT_FP32 && fft.NTT_GF31) {
    // GF31 data is located after the FP32 data.  Compute size of the FP32 data and trigs.
    defines += toDefine("DISTGF31",      FP32_DATA_SIZE(fft.shape.width, fft.shape.middle, fft.shape.height, in_place, pad_size) / 2);
    defines += toDefine("DISTWTRIGGF31", SMALLTRIG_FP32_DIST(fft.shape.width, wTrigMiddle, wTrigHeight, wTrigNH));
    defines += toDefine("DISTMTRIGGF31", MIDDLETRIG_FP32_DIST(fft.shape.width, fft.shape.middle, fft.shape.height));
    defines += toDefine("DISTHTRIGGF31", SMALLTRIGCOMBO_FP32_DIST(fft.shape.width, fft.shape.middle, fft.shape.height, fft.shape.nH()));
  }
  else if (fft.FFT_FP32 && fft.NTT_GF61) {
    // GF61 data is located after the FP32 data.  Compute size of the FP32 data and trigs.
    defines += toDefine("DISTGF61",      FP32_DATA_SIZE(fft.shape.width, fft.shape.middle, fft.shape.height, in_place, pad_size) / 2);
    defines += toDefine("DISTWTRIGGF61", SMALLTRIG_FP32_DIST(fft.shape.width, wTrigMiddle, wTrigHeight, wTrigNH));
    defines += toDefine("DISTMTRIGGF61", MIDDLETRIG_FP32_DIST(fft.shape.width, fft.shape.middle, fft.shape.height));
    defines += toDefine("DISTHTRIGGF61", SMALLTRIGCOMBO_FP32_DIST(fft.shape.width, fft.shape.middle, fft.shape.height, fft.shape.nH()));
  }
  else if (fft.NTT_GF31 && fft.NTT_GF61) {
    defines += toDefine("DISTGF31",      0);
    defines += toDefine("DISTWTRIGGF31", 0);
    defines += toDefine("DISTMTRIGGF31", 0);
    defines += toDefine("DISTHTRIGGF31", 0);
    // GF61 data is located after the GF31 data.  Compute size of the GF31 data and trigs.
    defines += toDefine("DISTGF61",      GF31_DATA_SIZE(fft.shape.width, fft.shape.middle, fft.shape.height, in_place, pad_size) / 2);
    defines += toDefine("DISTWTRIGGF61", SMALLTRIG_GF31_DIST(fft.shape.width, wTrigMiddle, wTrigHeight, wTrigNH));
    defines += toDefine("DISTMTRIGGF61", MIDDLETRIG_GF31_DIST(fft.shape.width, fft.shape.middle, fft.shape.height));
    defines += toDefine("DISTHTRIGGF61", SMALLTRIGCOMBO_GF31_DIST(fft.shape.width, fft.shape.middle, fft.shape.height, fft.shape.nH()));
  }
  else if (fft.NTT_GF31) {
    defines += toDefine("DISTGF31",      0);
    defines += toDefine("DISTWTRIGGF31", 0);
    defines += toDefine("DISTMTRIGGF31", 0);
    defines += toDefine("DISTHTRIGGF31", 0);
  }
  else if (fft.NTT_GF61) {
    defines += toDefine("DISTGF61",      0);
    defines += toDefine("DISTWTRIGGF61", 0);
    defines += toDefine("DISTMTRIGGF61", 0);
    defines += toDefine("DISTHTRIGGF61", 0);
  }

  // Calculate fractional bits-per-word = (E % N) / N * 2^64
  u32 const bpw_hi = (u64(E % N) << 32) / N;
  u32 const bpw_lo = (((u64(E % N) << 32) % N) << 32) / N;
  u64 bpw = (u64(bpw_hi) << 32) + bpw_lo;
  bpw--; // bpw must not be an exact value -- it must be less than exact value to get last biglit value right
  defines += toDefine("FRAC_BPW_HI", (u32) (bpw >> 32));
  defines += toDefine("FRAC_BPW_LO", (u32) bpw);

  return defines;
}

template<typename T>
pair<vector<T>, vector<T>> split(const vector<T>& v, const vector<u32>& select) {
  vector<T> a;
  vector<T> b;
  auto selIt = select.begin();
  u32 selNext = selIt == select.end() ? u32(-1) : *selIt;
  for (u32 i = 0; i < v.size(); ++i) {
    if (i == selNext) {
      b.push_back(v[i]);
      ++selIt;
      selNext = selIt == select.end() ? u32(-1) : *selIt;
    } else {
      a.push_back(v[i]);
    }
  }
  return {a, b};
}

RoeInfo roeStat(const vector<float>& roe) {
  double sumRoe = 0;
  double sum2Roe = 0;
  double maxRoe = 0;

  for (auto xf : roe) {
    double const x = xf;
    assert(x >= 0);
    maxRoe = max(x, maxRoe);
    sumRoe  += x;
    sum2Roe += x * x;
  }
  u32 const n = u32(roe.size());

  double const sdRoe = sqrt(n * sum2Roe - sumRoe * sumRoe) / n;
  double const meanRoe = sumRoe / n;

  return {n, maxRoe, meanRoe, sdRoe};
}

class IterationTimer {
  Timer timer;
  u64 kStart;

public:
  explicit IterationTimer(u64 kStart) : kStart(kStart) { }

  float reset(u64 k) {
    float const secs = float(timer.reset());

    u64 const its = max(u64(1), k - kStart);
    kStart = k;
    return secs / its;
  }
};

u32 baseCheckStep(u32 blockSize) {
  switch (blockSize) {
    case 200:  return 40'000;
    case 400:  return 160'000;
    case 500:  return 200'000;
    case 1000: return 1'000'000;
    default:
      assert(false);
      return 0;
  }
}

u32 checkStepForErrors(u32 blockSize, u32 nErrors) {
  u32 const step = baseCheckStep(blockSize);
  return nErrors ? step / 2 : step;
}

string toHex(u32 x) {
  char buf[16];
  snprintf(buf, sizeof(buf), "%08x", x);
  return buf;
}

string toHex(const vector<u32>& v) {
  string s;
  for (auto it = v.rbegin(), end = v.rend(); it != end; ++it) {
    s += toHex(*it);
  }
  return s;
}

string formatSecsPerIter(float secsPerIter) {
  char buf[64];
  float const usecsPerIter = secsPerIter * 1.0e6f;    // Convert to micro-seconds
  if (usecsPerIter > 1000.0f) {
    snprintf(buf, sizeof(buf), "%4.0f", usecsPerIter);
  } else {
    snprintf(buf, sizeof(buf), "%6.1f", usecsPerIter);
  }
  return string(buf);  
}

} // namespace

// --------

unique_ptr<Gpu> Gpu::make(u64 E, GpuCommon shared, FFTConfig fftConfig, const vector<KeyVal>& extraConf, bool logFftSize) {
  bool const lazyRiesel = fftConfig.NTT_RIESEL &&
                          shared.args->value("RIESEL_LAZY", 0);
  bool const goodThomas3 = shared.args->value("GOOD_THOMAS3", 0);
  bool const goodThomas7 = shared.args->value("GOOD_THOMAS7", 0);
  bool const goodThomas9 = shared.args->value("GOOD_THOMAS9", 0);
  bool const m19Field = shared.args->value("M19_FIELD", 0);
  if (goodThomas3 &&
      ((fftConfig.shape.fft_type != FFT31R2 &&
        fftConfig.shape.fft_type != FFT323161) ||
       (fftConfig.shape.middle != 3 && fftConfig.shape.middle != 6) ||
       fftConfig.shape.size() != 3u * 1024u * 1024u)) {
    throw "GOOD_THOMAS3 requires FFT31R2 or FFT323161 at a 3M middle-3/6 geometry";
  }
  if (goodThomas7 &&
      (fftConfig.shape.fft_type != FFT31R2 ||
       fftConfig.shape.middle != 7 ||
       fftConfig.shape.size() != 7u * 512u * 1024u)) {
    throw "GOOD_THOMAS7 requires FFT31R2 at a 3.5M middle-7 geometry";
  }
  if (goodThomas9 &&
      (fftConfig.shape.fft_type != FFT3261 ||
       fftConfig.shape.middle != 9 ||
       fftConfig.shape.size() != 9u * 512u * 1024u)) {
    throw "GOOD_THOMAS9 requires FFT3261 at a 4.5M middle-9 geometry";
  }
  if ((int(goodThomas3) + int(goodThomas7) + int(goodThomas9)) > 1) {
    throw "Good-Thomas experimental modes are mutually exclusive";
  }
  if ((goodThomas3 || goodThomas7) && lazyRiesel) {
    throw "Good-Thomas modes require an admissible ordinary field, not RIESEL_LAZY";
  }
  if (m19Field && (!(goodThomas3 || goodThomas7) ||
                   fftConfig.shape.fft_type != FFT31R2)) {
    throw "M19_FIELD requires FFT31R2 with a Good-Thomas mode";
  }
  if (lazyRiesel && (shared.args->value("RIESEL_PACKED", 0) ||
                     shared.args->value("RIESEL_FOUR", 0) ||
                     shared.args->value("RIESEL_FUSED", 0))) {
    throw "RIESEL_LAZY currently requires independent q fields and split carry";
  }
  // Experimental M31R2 representation: keep q0/q1 as independent residues,
  // but vectorize them into the low/high halves of each 64-bit storage word.
  // This is a two-stream M31 + (q0,q1) architecture, not M61 arithmetic.
  if (fftConfig.NTT_RIESEL && shared.args->value("RIESEL_PACKED", 0)) {
    fftConfig.NTT_RIESEL = false;
    fftConfig.NTT_GF61 = true;
  }
  // Four independent fields: retain the M31 and q0/q1 planes and add the
  // production M61 plane.  This is intentionally opt-in until the 154-bit
  // carry and checkpoint representation are complete.
  if (fftConfig.NTT_RIESEL &&
      (shared.args->value("RIESEL_FOUR", 0) || goodThomas3 || goodThomas7)) {
    fftConfig.NTT_GF61 = true;
  }
  return make_unique<Gpu>(shared, fftConfig, E, extraConf, logFftSize);
}

Gpu::~Gpu() {
  // Background tasks may have captured *this*, so wait until those are complete before destruction
  background->waitEmpty();
}

// Part of GPU initialization is to compute the default number of registers each kernel should target during compilation.
// Kernel register usage is critical for maximizing GPU occupancy.  The default values can be overrriden with command line arguments.
// This feature currently only works for the CUDA compiler.
// Most kernels have occupancy limited by register usage.  For reference, the following guidelines dictate where an "uptick" in occupancy occurs.
// If kernel threads=256, register crossovers are at 128, 80, 64, 48, 40
// If kernel threads=128, register crossovers are at 128, 96, 80, 72, 64, 56, 48, 40
// If kernel threads=64,  register crossovers are at 128, 112, 96, 88, 80, 72, 64, 56, 48, 40
string Gpu::numCudaRegisters([[maybe_unused]] enum WHICH_KERNEL which_kernel) {
#if CUDA_BACKEND
  int regs = 0;
  const char *use_override = "";
  // Allow command line to prefer the CUDA compiler's default number of registers
  if (args.value("NOREG", 0)) return string("");
  // Determine a kernel specific default maximum number of GPU registers (values set to -1 have not been tuned for best default value)
  switch (which_kernel) {
  case CARRYFUSED:         // Register usage depends on NW, the FFT/NTT type, and perhaps the long carry setting
    switch (fft.shape.fft_type) {
    case FFT64:
      regs = nW == 8 ? 80 : 64;
      use_override = "REGCF64";
      break;
    case FFT3161:
      regs = nW == 8 ? 96 : 64;         // Tested on 4090, nW=8, CUDA 13.0 (88 regs is possible without spilling but is slower)
      use_override = "REGCF3161";
      break;
    case FFT31R2:
      regs = nW == 8 ? 128 : 80;
      use_override = "REGCF31R2";
      break;
    case FFT3261:
      regs = nW == 8 ? 96 : 64;
      use_override = "REGCF3261";
      break;
    case FFT61:
      regs = nW == 8 ? 80 : 64;
      use_override = "REGCF61";
      break;
    case FFT323161:
      regs = nW == 8 ? 128 : 80;
      use_override = "REGCF323161";
      break;
    case FFT3231:
      regs = -1;
      use_override = "REGCF3231";
      break;
    case FFT6431:
      regs = -1;
      use_override = "REGCF6431";
      break;
    case FFT31:
      regs = -1;
      use_override = "REGCF31";
      break;
    case FFT32:
      regs = -1;
      use_override = "REGCF32";
      break;
    }
    break;
  case MIDIN:              // Register usage depends on MIDDLE and the FP32/FP64
    if (fft.FFT_FP64) {
      if (fft.shape.middle >= 16) regs = 96;
      else if (fft.shape.middle >= 14) regs = 88;
      else if (fft.shape.middle >= 13) regs = 80;
      else if (fft.shape.middle >= 10) regs = 72;
      else if (fft.shape.middle >= 7) regs = 64;
      else if (fft.shape.middle >= 5) regs = 56;
      else if (fft.shape.middle >= 4) regs = 48;
      else regs = -1;
      use_override = "REGMI64";
    } else {
      if (fft.shape.middle == 16) regs = 56;
      else if (fft.shape.middle == 8) regs = 40;
      else if (fft.shape.middle == 4) regs = 32;
      else regs = -1;
      use_override = "REGMI32";
    }
    break;
  case MIDIN31:            // Register usage depends on MIDDLE
    if (fft.shape.middle == 16) regs = 56;
    else if (fft.shape.middle == 8) regs = 48;         // Tested on 4090, CUDA 13.0 (40 regs is possible without spilling but is not measurably faster)
    else if (fft.shape.middle == 4) regs = 32;         // Tested on 5070Ti, CUDA 13.2.
    else regs = -1;  
    use_override = "REGMI31";
    break;
  case MIDIN61:            // Register usage depends on MIDDLE
    if (fft.shape.middle == 16) regs = 96;
    else if (fft.shape.middle == 8) regs = 64;         // Tested on 4090, CUDA 13.0
    else if (fft.shape.middle == 4) regs = -1;         // Tested on 5070Ti, CUDA 13.2 (48 regs is possible without spilling but is slower), best is -1.
    else regs = -1;  
    use_override = "REGMI61";
    break;
  case TAIL:               // Register usage depends on NH and the FP32/FP64 (assumes double-wide kernel)
    if (fft.FFT_FP64) {
      regs = nH == 8 ? 88 : 64;
      use_override = "REGTS64";
    } else {
      regs = nH == 8 ? 64 : 48;
      use_override = "REGTS32";
    }
    break;
  case TAIL31:             // Register usage depends on NH (assumes double-wide kernel)
    regs = nH == 8 ? -1 : 48;         // Tested on 4090, NH=8, CUDA 13.0.  Occupancy is limited by LDS memory use, register usage of 48 is possible, best is 64.
                                      // Tested on 5070Ti, NH=8, CUDA 13.2.  Occupancy is limited by LDS memory use, register usage of 56 is possible, best is -1.
    use_override = "REGTS31";
    break;
  case TAIL61:             // Register usage depends on NH (assumes double-wide kernel)
    regs = nH == 8 ? 96 : 64;         // Tested on 4090, nH=8, CUDA 13.0 (80 regs is possible without spilling but is slower)
    use_override = "REGTS61";
    break;
  case MIDOUT:             // Register usage depends on MIDDLE and the FFT/NTT type
    if (fft.FFT_FP64) {
      if (fft.shape.middle >= 15) regs = 96;
      else if (fft.shape.middle >= 14) regs = 88;
      else if (fft.shape.middle >= 11) regs = 80;
      else if (fft.shape.middle >= 10) regs = 72;
      else if (fft.shape.middle >= 7) regs = 64;
      else if (fft.shape.middle >= 5) regs = 56;
      else if (fft.shape.middle >= 4) regs = 48;
      else regs = -1;
      use_override = "REGMO64";
    } else {
      if (fft.shape.middle == 16) regs = 56;
      else if (fft.shape.middle == 8) regs = 40;
      else if (fft.shape.middle == 4) regs = 32;
      else regs = -1;
      use_override = "REGMO32";
    }
    break;
  case MIDOUT31:           // Register usage depends on MIDDLE
    if (fft.shape.middle == 16) regs = 48;
    else if (fft.shape.middle == 8) regs = 40;         // Tested on 4090, CUDA 13.0
    else if (fft.shape.middle == 4) regs = -1;         // Tested on 5070Ti, CUDA 13.2 (32 regs is possible without spilling but is slower), best is -1.
    else regs = -1;
    use_override = "REGMO31";
    break;
  case MIDOUT61:           // Register usage depends on MIDDLE
    if (fft.shape.middle == 16) regs = 96;
    else if (fft.shape.middle == 8) regs = 64;         // Tested on 4090, CUDA 13.0, best is 64.
                                                       // Tested on 5070Ti, CUDA 13.2, best is 72.
    else if (fft.shape.middle == 4) regs = 64;         // Tested on 5070Ti, CUDA 13.2 (48 regs is possible without spilling but is slower), best is 64.
    else regs = -1;
    use_override = "REGMO61";
    break;
  }
  // Get the optional override register count
  int const override_regs = args.value(use_override, 0);
  // If a specified override is small, use the count as a CUDA launch_bounds rather than a maximum register count
  if (override_regs && (override_regs > 0 && override_regs <= 16)) return string("-DCUDA_MIN_BLOCKS=") + to_string(override_regs) + " ";
  // If specified, override the default maximum register count
  if (override_regs) regs = override_regs;
  // Sometimes the results using CUDA compiler's default launch_bounds without setting an explicit launch bounds or maxrrregcount can't be beat
  if (regs == -1) return string("");
  // Format an explicit register count setting
  return string("--maxrregcount=") + to_string(regs) + " ";
#else
  return string("");
#endif
}

// Kernels are compiled one at a time, but OpenCL source files contain multiple kernels.  This routine set the #defines necessary so that only one kernel is compiled.
// While not strictly necessary, startup speed will be a bit faster if we do less compilations.
string Gpu::kernelDefines(enum WHICH_KERNEL_TYPE which_kernel) {
  string defines;
  // Determine the kernel specific #defines
  switch (which_kernel) {
  case KFP:         // FP64 or FP32 kernel
    defines += toDefine("FFT_FP64", (int) fft.FFT_FP64);
    defines += toDefine("FFT_FP32", (int) fft.FFT_FP32);
    defines += toDefine("NTT_GF31", 0);
    defines += toDefine("NTT_GF61", 0);
    break;
  case K31:         // GF1 kernel
    defines += toDefine("FFT_FP64", 0);
    defines += toDefine("FFT_FP32", 0);
    defines += toDefine("NTT_GF31", (int) fft.NTT_GF31);
    defines += toDefine("NTT_GF61", 0);
    break;
  case K61:         // GF61 kernel
    defines += toDefine("FFT_FP64", 0);
    defines += toDefine("FFT_FP32", 0);
    defines += toDefine("NTT_GF31", 0);
    defines += toDefine("NTT_GF61", (int) fft.NTT_GF61);
    break;
  case KALL:        // Kernels, like carryFused, that need all defines set properly
    defines += toDefine("FFT_FP64", (int) fft.FFT_FP64);
    defines += toDefine("FFT_FP32", (int) fft.FFT_FP32);
    defines += toDefine("NTT_GF31", (int) fft.NTT_GF31);
    defines += toDefine("NTT_GF61", (int) fft.NTT_GF61);
    break;
  }
  return defines + " ";
}

enum {
ROE_SIZE = 100000,
CARRY_SIZE = 100000
};

// A deliberately small, independently scheduled M31 cyclic square for the
// folded error syndrome.  The input is already full-transform weighted, so its
// first width edge does not apply Crandall--Fagin weights again.
class FoldTransform {
  u32 factor;
  u32 sideN;
  u32 hn;
  FFTConfig fft;
  u32 width;
  u32 middle;
  u32 height;
  u32 bigHeight;
  u32 nw;
  u32 nh;
  bool tailSingleWide{};
  bool tailSingleKernel{};
  u32 inPlace{};
  u32 padSize{};
  u32 wmul{};
  Queue queue;
  KernelCompiler compiler;
  Kernel kCompact;
  Kernel kP;
  Kernel kMidIn;
  Kernel kTail;
  Kernel kMidOut;
  Kernel kW;
  TrigPtr trigH;
  TrigPtr trigM;
  TrigPtr trigW;
  Buffer<u32> compact;
  Buffer<double> data;
  EventHolder inputRead;
  bool pending{};
  bool zeroWidth{};
  bool primed{};

  static u32 foldFactor(Args* args) {
    u32 const factor = args->value("FOLD_FACTOR", 8);
    if (factor != 4 && factor != 8 && factor != 16) {
      throw runtime_error("FOLD_FACTOR must be 4, 8, or 16");
    }
    return factor;
  }

  static vector<KeyVal> config(u32 factor) {
    return {{"INPLACE", "1"}, {"PAD", "0"}, {"TAIL_KERNELS", "2"},
            {"MULTI_Q", "0"}, {"GRAPHS", "0"},
            {"PARITY_CORRECT", "0"}, {"PARITY_PREPARED", "0"},
            {"PARITY_PACKED", "0"}, {"PARITY_LAZY", "0"}, {"FOLD_SYNDROME", "0"},
            {"FOLD_TRANSFORM", "1"}, {"FOLD_FACTOR", to_string(factor)}};
  }

  static string fftSpec(Args* args) {
    if (foldFactor(args) == 4) {
      if (args->value("FOLD_SHAPE", 0) != 0) {
        throw runtime_error("FOLD_FACTOR=4 currently requires FOLD_SHAPE=0");
      }
      return "52:512:2:512:202";
    }
    if (foldFactor(args) == 16) {
      if (args->value("FOLD_SHAPE", 0) != 0) {
        throw runtime_error("FOLD_FACTOR=16 currently requires FOLD_SHAPE=0");
      }
      return "52:256:2:256:202";
    }
    switch (args->value("FOLD_SHAPE", 0)) {
    case 0: return "52:256:4:256:202";
    case 1: return "52:256:2:512:202";
    case 2: return "52:512:2:256:202";
    default: throw runtime_error("FOLD_SHAPE must be 0, 1, or 2");
    }
  }

public:
  FoldTransform(GpuCommon shared, Profile& profile, u32 mainN) :
    factor{foldFactor(shared.args)},
    sideN{mainN / factor},
    hn{sideN / 2},
    fft{fftSpec(shared.args)},
    width{fft.shape.width},
    middle{fft.shape.middle},
    height{fft.shape.height},
    bigHeight{middle * height},
    nw{fft.shape.nW()},
    nh{fft.shape.nH()},
    queue{*shared.context, shared.args->profile, true},
    compiler{*shared.args, shared.context,
             clDefines(*shared.args, shared.context->deviceId(), fft, config(factor),
                       u64(sideN) * 4, false, tailSingleWide,
                       tailSingleKernel, inPlace, padSize, wmul)},
    kCompact{"foldCompact", &compiler, profile.make("foldCompact"), &queue,
       "foldp.cl", "foldCompact", hn,
       "-DFFT_FP64=0 -DFFT_FP32=0 -DNTT_GF31=1 -DNTT_GF61=0 "},
    kP{"foldP", &compiler, profile.make("foldP"), &queue,
       "foldp.cl", "foldP", hn / nw,
       "-DFFT_FP64=0 -DFFT_FP32=0 -DNTT_GF31=1 -DNTT_GF61=0 "},
    kMidIn{"foldMidIn", &compiler, profile.make("foldMidIn"), &queue,
       "fftmiddlein.cl", "fftMiddleInGF31", hn / (bigHeight / height),
       "-DFFT_FP64=0 -DFFT_FP32=0 -DNTT_GF31=1 -DNTT_GF61=0 --maxrregcount=32 "},
    kTail{"foldTail", &compiler, profile.make("foldTail"), &queue,
       "tailsquare.cl", "tailSquareGF31", hn / nh,
       "-DFFT_FP64=0 -DFFT_FP32=0 -DNTT_GF31=1 -DNTT_GF61=0 --maxrregcount=48 "},
    kMidOut{"foldMidOut", &compiler, profile.make("foldMidOut"), &queue,
       "fftmiddleout.cl", "fftMiddleOutGF31", hn / (bigHeight / height),
       "-DFFT_FP64=0 -DFFT_FP32=0 -DNTT_GF31=1 -DNTT_GF61=0 "},
    kW{"foldW", &compiler, profile.make("foldW"), &queue,
       "fftw.cl", "fftWGF31", hn / nw,
       "-DFFT_FP64=0 -DFFT_FP32=0 -DNTT_GF31=1 -DNTT_GF61=0 "},
    trigH{shared.bufCache->smallTrigCombo(shared.args, fft, width, middle,
                                          height, nh, tailSingleWide)},
    trigM{shared.bufCache->middleTrig(shared.args, fft, height, middle, width)},
    trigW{shared.bufCache->smallTrig(shared.args, fft, width, nw, middle,
                                     height, nh, tailSingleWide)},
    compact{profile.make("foldCompactData"), &queue, factor == 16 ? sideN : 1},
    data{profile.make("foldData"), &queue,
         TOTAL_DATA_SIZE(fft, width, middle, height, inPlace, padSize)},
    zeroWidth{shared.args->value("FOLD_ZERO_WIDTH", 0) != 0}
  {
#if CUDA_BACKEND
    // The main M61 path owns the critical integer pipeline.  Let the compact
    // sidecar fill gaps at low priority; only its short input edge is a hard
    // dependency, and protectInput() enforces that boundary explicitly.
    queue.setCudaPriority(-1);
#endif
    kP.setFixedArgs(2, trigW);
    kMidIn.setFixedArgs(3, trigM);
    kTail.setFixedArgs(3, trigH);
    kMidOut.setFixedArgs(3, trigM);
    kW.setFixedArgs(2, trigW);
  }

  void wait(Queue& main) {
    if (!pending) return;
    EventHolder done = queue.createSyncEvent();
    main.waitForSyncEvent(&done);
    pending = false;
  }

  // Carry may overwrite the compact input as soon as foldP has consumed it;
  // it need not wait for the middle, square, and inverse edge.  The side queue
  // remains in order, so its private data buffer can safely pipeline iterations.
  void protectInput(Queue& main) {
    if (inputRead) main.waitForSyncEvent(&inputRead);
  }

  void launch(Queue& main, Buffer<u32>& input) {
    EventHolder ready = main.createSyncEvent();
    queue.waitForSyncEvent(&ready);
    if (factor == 16) {
      kCompact(compact, input);
      inputRead = queue.createSyncEvent();
      kP(data, compact);
    } else {
      // Prime the private buffer once so the lower-bound mode still exercises
      // dense, valid residues.  Subsequent iterations deliberately grant both
      // width transforms zero cost; this is not a correct side recurrence.
      if (!zeroWidth || !primed) kP(data, input);
      inputRead = queue.createSyncEvent();
    }
    kMidIn(data, data, 0);
    kTail(data, data, 0);
    kMidOut(data, data, 0);
    if (!zeroWidth) kW(data, data);
    primed = true;
    pending = true;
  }

  Buffer<double>& output() { return data; }
};

Gpu::Gpu(GpuCommon s, FFTConfig fft, u64 E, const vector<KeyVal>& extraConf, bool logFftSize) :
  shared(s),
  background{shared.background},
  args{*shared.args},
  E(E),
  N(fft.shape.size()),
  fft(fft),
  WIDTH(fft.shape.width),
  SMALL_H(fft.shape.height),
  BIG_H(SMALL_H * fft.shape.middle),
  hN(N / 2),
  nW(args.value("WIDTH_NW", fft.shape.nW())),
  nH(fft.shape.nH()),
  carryLen{carryLength(args, fft)},
  useLongCarry{args.carry == CARRY_64},
  useMiddleCarry{args.value("MIDCARRY_FUSED", 0) == 1},
  useSplitCarry{args.value("MIDCARRY_FUSED", 0) == 2},
  parityCorrect{args.value("PARITY_CORRECT", 0) &&
                (fft.shape.fft_type == FFT3261 || fft.shape.fft_type == FFT323161)},
  parityPrepared{args.value("PARITY_PREPARED", 0) == 1 && parityCorrect},
  parityDirect{args.value("PARITY_PREPARED", 0) == 2 && parityCorrect},
  parityPhysical{args.value("PARITY_PREPARED", 0) == 3 && parityCorrect},
  parityPacked{args.value("PARITY_PACKED", 0) != 0 && parityCorrect},
  foldSyndrome{args.value("FOLD_SYNDROME", 0) != 0},
  foldValidate{args.value("FOLD_SYNDROME", 0) > 1},
  foldTransformEnabled{args.value("FOLD_TRANSFORM", 0) != 0},
  queue{*shared.context, args.profile},
      
  compiler{args, shared.context, clDefines(args, shared.context->deviceId(), fft, extraConf, E, logFftSize, tail_single_wide, tail_single_kernel, in_place, pad_size, wmul)},

#define K(name, ...) name(#name, &compiler, profile.make(#name), &queue, __VA_ARGS__)

  K(kfftMidIn,             "fftmiddlein.cl",  "fftMiddleIn",  hN / (BIG_H / SMALL_H), kernelDefines(KFP) + numCudaRegisters(MIDIN)),
  K(kfftHin,               "ffthin.cl",  "fftHin",  hN / nH, kernelDefines(KFP)),
  K(ktailSquareZero,       "tailsquare.cl", "tailSquareZero", SMALL_H / nH * 2, kernelDefines(KFP)),
  K(ktailSquare,           "tailsquare.cl", "tailSquare",
                                               !tail_single_wide && !tail_single_kernel ? hN / nH - SMALL_H / nH * 2 : // Double-wide tailSquare with two kernels
                                               !tail_single_wide ? hN / nH :                                           // Double-wide tailSquare with one kernel
                                               !tail_single_kernel ? hN / nH / 2 - SMALL_H / nH :                      // Single-wide tailSquare with two kernels
                                               hN / nH / 2, kernelDefines(KFP) + numCudaRegisters(TAIL)),    // Single-wide tailSquare with one kernel
  K(ktailMul,              "tailmul.cl", "tailMul", hN / nH / 2, kernelDefines(KFP)),
  K(ktailMulLow,           "tailmul.cl", "tailMul", hN / nH / 2, kernelDefines(KFP) + "-DMUL_LOW=1"),
  K(kfftMidOut,            "fftmiddleout.cl", "fftMiddleOut", hN / (BIG_H / SMALL_H), kernelDefines(KFP) + numCudaRegisters(MIDOUT)),
  K(kfftW,                 "fftw.cl", "fftW", hN / nW, kernelDefines(KFP)),

  K(kfftMidInGF31,         "fftmiddlein.cl",  "fftMiddleInGF31",  hN / (BIG_H / SMALL_H), kernelDefines(K31) + numCudaRegisters(MIDIN31)),
  K(kfftHinGF31,           "ffthin.cl",  "fftHinGF31",  hN / nH, kernelDefines(K31)),
  K(ktailSquareZeroGF31,   "tailsquare.cl", "tailSquareZeroGF31", SMALL_H / nH * 2, kernelDefines(K31)),
  K(ktailSquareGF31,       "tailsquare.cl", "tailSquareGF31",
                                               !tail_single_wide && !tail_single_kernel ? hN / nH - SMALL_H / nH * 2 : // Double-wide tailSquare with two kernels
                                               !tail_single_wide ? hN / nH :                                           // Double-wide tailSquare with one kernel
                                               !tail_single_kernel ? hN / nH / 2 - SMALL_H / nH :                      // Single-wide tailSquare with two kernels
                                               hN / nH / 2, kernelDefines(K31) + numCudaRegisters(TAIL31)),  // Single-wide tailSquare with one kernel
  // FUSED31: fused middle-in + tail for the M31 plane; (WIDTH/2+1) pair-resident
  // blocks of G_H*2 threads, 64-KiB dynamic-shared tile, max shared carveout.
  K(kFusedMidTail31,       "fusedmidtail31.cl", "fusedMidTail31",
                                               (WIDTH / 2 + 1) * (SMALL_H / nH * 4),
                                               kernelDefines(K31), 2 * fft.shape.middle * SMALL_H * 8,
                                               args.value("FUSED31", 0) ? 100 : 0),
  K(ktailMulGF31,          "tailmul.cl", "tailMulGF31", hN / nH / 2, kernelDefines(K31)),
  K(ktailMulLowGF31,       "tailmul.cl", "tailMulGF31", hN / nH / 2, kernelDefines(K31) + "-DMUL_LOW=1"),
  K(kfftMidOutGF31,        "fftmiddleout.cl", "fftMiddleOutGF31", hN / (BIG_H / SMALL_H), kernelDefines(K31) + numCudaRegisters(MIDOUT31)),
  K(kfftWGF31,             "fftw.cl", "fftWGF31", hN / nW, kernelDefines(K31)),
  K(kfftWOut31,            "fftw.cl", "fftWOut31", hN / nW, kernelDefines(K31)),

  K(kfftMidInR0,           "fftmiddlein.cl", "fftMiddleInGF31", hN / (BIG_H / SMALL_H), kernelDefines(K31) + "-DRIESEL_FIELD=1 " + numCudaRegisters(MIDIN31)),
  K(kfftHinR0,             "ffthin.cl", "fftHinGF31", hN / nH, kernelDefines(K31) + "-DRIESEL_FIELD=1 "),
  K(ktailSquareZeroR0,     "tailsquare.cl", "tailSquareZeroGF31", SMALL_H / nH * 2, kernelDefines(K31) + "-DRIESEL_FIELD=1 "),
  K(ktailSquareR0,         "tailsquare.cl", "tailSquareGF31",
                                               !tail_single_wide && !tail_single_kernel ? hN / nH - SMALL_H / nH * 2 :
                                               !tail_single_wide ? hN / nH :
                                               !tail_single_kernel ? hN / nH / 2 - SMALL_H / nH :
                                               hN / nH / 2, kernelDefines(K31) + "-DRIESEL_FIELD=1 " + numCudaRegisters(TAIL31)),
  K(ktailMulR0,            "tailmul.cl", "tailMulGF31", hN / nH / 2, kernelDefines(K31) + "-DRIESEL_FIELD=1 "),
  K(ktailMulLowR0,         "tailmul.cl", "tailMulGF31", hN / nH / 2, kernelDefines(K31) + "-DRIESEL_FIELD=1 -DMUL_LOW=1 "),
  K(kfftMidOutR0,          "fftmiddleout.cl", "fftMiddleOutGF31", hN / (BIG_H / SMALL_H), kernelDefines(K31) + "-DRIESEL_FIELD=1 " + numCudaRegisters(MIDOUT31)),
  K(kfftWR0,               "fftw.cl", "fftWGF31", hN / nW, kernelDefines(K31) + "-DRIESEL_FIELD=1 "),
  K(kfftPR0,               "fftp.cl", "fftP", hN / nW, kernelDefines(K31) + "-DRIESEL_FIELD=1 "),

  K(kfftMidInR1,           "fftmiddlein.cl", "fftMiddleInGF31", hN / (BIG_H / SMALL_H), kernelDefines(K31) + "-DRIESEL_FIELD=2 " + numCudaRegisters(MIDIN31)),
  K(kfftHinR1,             "ffthin.cl", "fftHinGF31", hN / nH, kernelDefines(K31) + "-DRIESEL_FIELD=2 "),
  K(ktailSquareZeroR1,     "tailsquare.cl", "tailSquareZeroGF31", SMALL_H / nH * 2, kernelDefines(K31) + "-DRIESEL_FIELD=2 "),
  K(ktailSquareR1,         "tailsquare.cl", "tailSquareGF31",
                                               !tail_single_wide && !tail_single_kernel ? hN / nH - SMALL_H / nH * 2 :
                                               !tail_single_wide ? hN / nH :
                                               !tail_single_kernel ? hN / nH / 2 - SMALL_H / nH :
                                               hN / nH / 2, kernelDefines(K31) + "-DRIESEL_FIELD=2 " + numCudaRegisters(TAIL31)),
  K(ktailMulR1,            "tailmul.cl", "tailMulGF31", hN / nH / 2, kernelDefines(K31) + "-DRIESEL_FIELD=2 "),
  K(ktailMulLowR1,         "tailmul.cl", "tailMulGF31", hN / nH / 2, kernelDefines(K31) + "-DRIESEL_FIELD=2 -DMUL_LOW=1 "),
  K(kfftMidOutR1,          "fftmiddleout.cl", "fftMiddleOutGF31", hN / (BIG_H / SMALL_H), kernelDefines(K31) + "-DRIESEL_FIELD=2 " + numCudaRegisters(MIDOUT31)),
  K(kfftWR1,               "fftw.cl", "fftWGF31", hN / nW, kernelDefines(K31) + "-DRIESEL_FIELD=2 "),
  K(kfftPR1,               "fftp.cl", "fftP", hN / nW, kernelDefines(K31) + "-DRIESEL_FIELD=2 "),
  K(kfftP31R2,             "fftp.cl", "fftP", hN / nW, kernelDefines(K31)),
  K(kfftPRPair,            "fftp.cl", "fftP", hN / nW, kernelDefines(K61)),

  K(kfftMidInGF61,         "fftmiddlein.cl",  "fftMiddleInGF61",  hN / (BIG_H / SMALL_H), kernelDefines(K61) + numCudaRegisters(MIDIN61)),
  K(kfftHinGF61,           "ffthin.cl",  "fftHinGF61",  hN / nH, kernelDefines(K61)),
  K(ktailSquareZeroGF61,   "tailsquare.cl", "tailSquareZeroGF61", SMALL_H / nH * 2, kernelDefines(K61)),
  K(ktailSquareGF61,       "tailsquare.cl", "tailSquareGF61",
                                               !tail_single_wide && !tail_single_kernel ? hN / nH - SMALL_H / nH * 2 : // Double-wide tailSquare with two kernels
                                               !tail_single_wide ? hN / nH :                                           // Double-wide tailSquare with one kernel
                                               !tail_single_kernel ? hN / nH / 2 - SMALL_H / nH :                      // Single-wide tailSquare with two kernels
                                               hN / nH / 2, kernelDefines(K61) + numCudaRegisters(TAIL61)),  // Single-wide tailSquare with one kernel
  K(ktailMulGF61,          "tailmul.cl", "tailMulGF61", hN / nH / 2, kernelDefines(K61)),
  K(ktailMulLowGF61,       "tailmul.cl", "tailMulGF61", hN / nH / 2, kernelDefines(K61) + "-DMUL_LOW=1"),
  K(kfftMidOutGF61,        "fftmiddleout.cl", "fftMiddleOutGF61", hN / (BIG_H / SMALL_H) / (args.value("ASYNC_MID61", 0) ? args.value("ASYNC_TILES", 4) : 1), kernelDefines(K61) + numCudaRegisters(MIDOUT61), 0, args.value("ASYNC_MID61", 0) ? 100 : 0),
  K(kfftWGF61,             "fftw.cl", "fftWGF61", hN / nW, kernelDefines(K61)),

  K(kfftP,                 "fftp.cl", "fftP", hN / nW, kernelDefines(KALL)),
  K(kCarryA,               "carry.cl", "carry", hN / carryLen, kernelDefines(KALL)),
  K(kCarryAROE,            "carry.cl", "carry", hN / carryLen, kernelDefines(KALL) + "-DROE=1"),
  K(kCarryAParity,         "carry.cl", "carry", hN / carryLen, kernelDefines(KALL) + "-DPARITY_SQUARE=1 " + (parityPrepared || parityDirect || parityPacked ? "-DPARITY_PREPARED=1 " : "") + (parityDirect ? "-DPARITY_DIRECT=1 " : "") + (parityPhysical ? "-DPARITY_PHYSICAL=1 " : "") + (parityPacked ? "-DPARITY_PACKED=1 " : "")),
  K(kCarryAROEParity,      "carry.cl", "carry", hN / carryLen, kernelDefines(KALL) + "-DPARITY_SQUARE=1 -DROE=1 " + (parityPrepared || parityDirect || parityPacked ? "-DPARITY_PREPARED=1 " : "") + (parityDirect ? "-DPARITY_DIRECT=1 " : "") + (parityPhysical ? "-DPARITY_PHYSICAL=1 " : "") + (parityPacked ? "-DPARITY_PACKED=1 " : "")),
  K(kCarryM,               "carry.cl", "carry", hN / carryLen, kernelDefines(KALL) + "-DMUL3=1"),
  K(kCarryMROE,            "carry.cl", "carry", hN / carryLen, kernelDefines(KALL) + "-DMUL3=1 -DROE=1"),
  K(kCarryLL,              "carry.cl", "carry", hN / carryLen, kernelDefines(KALL) + "-DLL=1"),
  K(kCarryFused,           "carryfused.cl", "carryFused", WIDTH * (BIG_H + wmul) / nW, kernelDefines(KALL) + numCudaRegisters(CARRYFUSED) + (parityCorrect ? "-DPARITY_SQUARE=1 " : "") + (parityPrepared || parityDirect || parityPacked ? "-DPARITY_PREPARED=1 " : "") + (parityDirect ? "-DPARITY_DIRECT=1 " : "") + (parityPhysical ? "-DPARITY_PHYSICAL=1 " : "") + (parityPacked ? "-DPARITY_PACKED=1 " : "") + (foldValidate ? "-DFOLD_SYNDROME_VALIDATE=1 " : "")),
  K(kCarryFusedROE,        "carryfused.cl", "carryFused", WIDTH * (BIG_H + wmul) / nW, kernelDefines(KALL) + numCudaRegisters(CARRYFUSED) + "-DROE=1 " + (parityCorrect ? "-DPARITY_SQUARE=1 " : "") + (parityPrepared || parityDirect || parityPacked ? "-DPARITY_PREPARED=1 " : "") + (parityDirect ? "-DPARITY_DIRECT=1 " : "") + (parityPhysical ? "-DPARITY_PHYSICAL=1 " : "") + (parityPacked ? "-DPARITY_PACKED=1 " : "") + (foldValidate ? "-DFOLD_SYNDROME_VALIDATE=1 " : "")),
  K(kCarryFusedMul,        "carryfused.cl", "carryFused", WIDTH * (BIG_H + wmul) / nW, kernelDefines(KALL) + numCudaRegisters(CARRYFUSED) + "-DMUL3=1"),
  K(kCarryFusedMulROE,     "carryfused.cl", "carryFused", WIDTH * (BIG_H + wmul) / nW, kernelDefines(KALL) + numCudaRegisters(CARRYFUSED) + "-DMUL3=1 -DROE=1"),
  K(kCarryFusedLL,         "carryfused.cl", "carryFused", WIDTH * (BIG_H + wmul) / nW, kernelDefines(KALL) + numCudaRegisters(CARRYFUSED) + "-DLL=1 " + (parityCorrect ? "-DPARITY_SQUARE=1 " : "") + (parityPrepared || parityDirect || parityPacked ? "-DPARITY_PREPARED=1 " : "") + (parityDirect ? "-DPARITY_DIRECT=1 " : "") + (parityPhysical ? "-DPARITY_PHYSICAL=1 " : "") + (parityPacked ? "-DPARITY_PACKED=1 " : "") + (foldValidate ? "-DFOLD_SYNDROME_VALIDATE=1 " : "")),
  K(kCarryMiddleOut,       "carrymiddle.cl", "carryMiddleOut", SMALL_H * 512, kernelDefines(KALL) + "-DLDSPAD_W=0", 8 * 512 * (8 + 16)),
  K(kCarryMiddleOutROE,    "carrymiddle.cl", "carryMiddleOut", SMALL_H * 512, kernelDefines(KALL) + "-DLDSPAD_W=0 -DROE=1", 8 * 512 * (8 + 16)),
  K(kCarryMiddleIn,        "carrymiddle.cl", "carryMiddleIn", SMALL_H * 512, kernelDefines(KALL) + "-DLDSPAD_W=0", 8 * 512 * (8 + 16)),
  K(kCarrySplitOut,        "carrymiddle.cl", "carryMiddleOut", SMALL_H * 512, kernelDefines(KALL) + "-DLDSPAD_W=0 -DCM_SKIP_MIDDLE=1", 8 * 512 * (8 + 16)),
  K(kCarrySplitOutROE,     "carrymiddle.cl", "carryMiddleOut", SMALL_H * 512, kernelDefines(KALL) + "-DLDSPAD_W=0 -DCM_SKIP_MIDDLE=1 -DROE=1", 8 * 512 * (8 + 16)),
  K(kCarrySplitIn,         "carrymiddle.cl", "carryMiddleIn", SMALL_H * 512, kernelDefines(KALL) + "-DLDSPAD_W=0 -DCM_SKIP_MIDDLE=1", 8 * 512 * (8 + 16)),

  K(carryB,                "carryb.cl", "carryB",   hN / carryLen, kernelDefines(KALL)),

  // 64
  K(transpIn,  "transpose.cl", "transposeIn",  hN / 64),
  K(transpOut, "transpose.cl", "transposeOut", hN / 64),

  K(readResidue, "etc.cl", "readResidue", 32, "-DREADRESIDUE=1"),
  K(parityInit, "etc.cl", "parityInit", hN, kernelDefines(KALL) + "-DPARITY_INIT=1 " + (parityPrepared || parityDirect || parityPacked ? "-DPARITY_PREPARED=1 " : "") + (parityDirect ? "-DPARITY_DIRECT=1 " : "") + (parityPhysical ? "-DPARITY_PHYSICAL=1 " : "") + (parityPacked ? "-DPARITY_PACKED=1 " : "")),
  K(parityPrepare, "etc.cl", "parityPrepare", hN, kernelDefines(KALL) + "-DPARITY_PREPARE=1 " + (parityPrepared || parityPacked ? "-DPARITY_PREPARED=1 " : "") + (parityPacked ? "-DPARITY_PACKED=1 " : "")),
  K(kFoldValidate, "foldsyndrome.cl", "foldValidate", N / 16, kernelDefines(KALL) + (foldValidate ? "-DFOLD_SYNDROME_VALIDATE=1 " : "")),

  // 256
  K(kernIsEqual, "etc.cl", "isEqual", 256 * 256, "-DISEQUAL=1"),
  K(sum64,       "etc.cl", "sum64",   256 * 256, "-DSUM64=1"),

  K(testTrig, "selftest.cl", "testTrig", 256 * 256),
  K(testFFT4, "selftest.cl", "testFFT4", 256),
  K(testFFT14, "selftest.cl", "testFFT14", 256),
  K(testFFT15, "selftest.cl", "testFFT15", 256),
  K(testFFT, "selftest.cl", "testFFT", 256),
  K(testTime, "selftest.cl", "testTime", 4096 * 64),

#undef K

  bufTrigH{shared.bufCache->smallTrigCombo(shared.args, fft, WIDTH, fft.shape.middle, SMALL_H, nH, tail_single_wide)},
  bufTrigM{shared.bufCache->middleTrig(shared.args, fft, SMALL_H, BIG_H / SMALL_H, WIDTH)},
  bufTrigW{shared.bufCache->smallTrig(shared.args, fft, WIDTH, nW, fft.shape.middle, SMALL_H, nH, tail_single_wide)},

  weights{genWeights(fft, E, WIDTH, BIG_H, nW,
                     isNvidiaGpu(shared.context->deviceId()),
                     args.value("RIESEL_LAZY", 0),
                     args.value("GOLD_PAIR", 0),
                     args.value("M19_FIELD", 0))},
  bufConstWeights{shared.context, std::move(weights.weightsConstIF)},
  bufWeights{shared.context,      std::move(weights.weightsIF)},

#define BUF(name, ...) name{profile.make(#name), &queue, __VA_ARGS__}

  // GPU Buffers containing integer data.  Since this buffer is type i64, if fft.WordSize < 8 then we need less memory allocated.
  BUF(bufData, N * fft.WordSize / sizeof(Word)),
  BUF(bufAux, N * fft.WordSize / sizeof(Word)),
  BUF(bufCheck, N * fft.WordSize / sizeof(Word)),
  // Every double-word (i.e. N/2) produces one carry. In addition we may have one extra group thus WIDTH more carries.
  // FFT31R2 fused carry uses a small private spill area for workgroup zero so
  // it cannot overwrite line zero before the duplicate wraparound workgroup
  // has read the original transform data.
  BUF(bufCarry, N / 2 + WIDTH + (fft.NTT_RIESEL ? 2 * WIDTH * wmul : 0)),
  BUF(bufReady, (N / 2 + WIDTH) / 32), // Every wavefront (32 or 64 lanes) needs to signal "carry is ready"
  BUF(bufParityA, parityCorrect ? (parityPacked ? hN / 16 : hN) : 1),
  BUF(bufParityB, parityCorrect ? (parityPacked ? hN / 16 : hN) : 1),
  BUF(bufParityExpected, parityPrepared || parityPacked ? (parityPacked ? hN / 32 : hN) : 1),
  parityIn{&bufParityA},
  parityOut{&bufParityB},
  // Factor four emits twice as many packed GF31 values as the original
  // factor-eight fold.  Factor sixteen still uses the N/8 partial fold and
  // compacts it on the side queue.
  BUF(bufFolded, foldSyndrome ?
      N / (args.value("FOLD_FACTOR", 8) == 4 ? 4 : 8) : 1),
  BUF(bufFoldWords, foldValidate ? N : 1),
  BUF(bufFoldMismatches, 1),

  BUF(bufSmallOut, 256),
  BUF(bufSumOut,     1),
  BUF(bufTrue,       1),
  BUF(bufROE, ROE_SIZE + 2),
  BUF(bufStatsCarry, CARRY_SIZE + 2),

  BUF(buf1, TOTAL_DATA_SIZE(fft, WIDTH, fft.shape.middle, SMALL_H, in_place, pad_size)),
  BUF(buf2, TOTAL_DATA_SIZE(fft, WIDTH, fft.shape.middle, SMALL_H, in_place, pad_size)),
  BUF(buf3, TOTAL_DATA_SIZE(fft, WIDTH, fft.shape.middle, SMALL_H, in_place, pad_size)),
  // One flat GF31 value per coefficient pair for the detached M31 width edge.
  BUF(bufDetach31, (fft.shape.fft_type == FFT3161 && args.value("DETACHED_M31_EDGE", 0))
                   ? size_t(WIDTH) * BIG_H : 16),
#undef BUF

  statsBits{u32(args.value("STATS", 0))},
  timeBufVect{profile.make("proofBufVect")},

  recorded_kernels{},
  recorded_kernel_args{},

  use_graphs{},
  graph_square{}
{    
  float const bitsPerWord = E / float(N);
  if (logFftSize) {
    log("FFT: %s %s (%.2f bpw)\n", numberK(N).c_str(), fft.spec().c_str(), bitsPerWord);

    // Sometimes we do want to run a FFT beyond a reasonable BPW (e.g. during -ztune), and these situations
    // coincide with logFftSize == false
    if (fft.maxExp() < E) {
      log("Warning: %s (max %" PRIu64 ") may be too small for %" PRIu64 "\n", fft.spec().c_str(), fft.maxExp(), E);
    }
  }

  if (bitsPerWord < fft.minBpw()) {
    log("FFT size too large for exponent (%.2f bits/word < %.2f bits/word).\n", bitsPerWord, fft.minBpw());
    throw "FFT size too large";
  }

  // The independent field queues are currently substantially faster for M31R2.
  // Keep the fused three-prime carry available for experiments without making
  // its present performance cost the default.
  useLongCarry = useLongCarry ||
                 (fft.shape.fft_type == FFT31R2 && args.carry == CARRY_AUTO &&
                  !args.value("RIESEL_FUSED", 0)) ||
                 args.value("GOOD_THOMAS9", 0) ||
                 (bitsPerWord < 10.0);

  if (useLongCarry) { log("Using long carry!\n"); }

  if (foldSyndrome) {
    if (fft.shape.fft_type != FFT3261 || WIDTH != 512 ||
        fft.shape.middle != 8 || SMALL_H != 512 || !in_place ||
        useLongCarry) {
      throw std::runtime_error(
        "FOLD_SYNDROME requires in-place FFT3261 512:8:512 and short carry");
    }
    log(foldValidate ? "Using folded M31 sidecar with standalone validation\n" :
                       "Using folded M31 sidecar\n");
  }
  if (foldTransformEnabled && !foldSyndrome) {
    throw std::runtime_error("FOLD_TRANSFORM requires FOLD_SYNDROME=1");
  }
  if (args.value("FOLD_FACTOR", 8) == 16 && (!foldTransformEnabled || !parityPacked)) {
    throw std::runtime_error("FOLD_FACTOR=16 requires FOLD_TRANSFORM=1 and PARITY_PACKED=1");
  }
  if (args.value("FOLD_FACTOR", 8) != 8 && foldValidate) {
    throw std::runtime_error("FOLD_SYNDROME=2 validation currently supports only FOLD_FACTOR=8");
  }
  if (args.value("PARITY_PREPARED", 0) && !parityCorrect) {
    throw std::runtime_error("PARITY_PREPARED requires the FFT3261/FFT323161 PARITY_CORRECT path");
  }
  if (args.value("PARITY_PREPARED", 0) < 0 || args.value("PARITY_PREPARED", 0) > 3) {
    throw std::runtime_error("PARITY_PREPARED must be 0, 1 (overlap), 2 (direct scatter), or 3 (physical natural)");
  }
  if (args.value("PARITY_PACKED", 0) && !parityCorrect) {
    throw std::runtime_error("PARITY_PACKED requires the FFT3261/FFT323161 PARITY_CORRECT path");
  }
  if (args.value("PARITY_PACKED", 0) && args.value("PARITY_PREPARED", 0)) {
    throw std::runtime_error("PARITY_PACKED is a separate prepared-parity mode; PARITY_PREPARED must be 0");
  }
#if !CUDA_BACKEND
  if (parityPacked) {
    throw std::runtime_error("PARITY_PACKED currently requires the CUDA backend");
  }
#endif
  if (args.value("PARITY_LAZY", 0) &&
      (fft.shape.fft_type != FFT3261 || (!parityPrepared && !parityDirect && !parityPhysical))) {
    throw std::runtime_error("PARITY_LAZY requires FFT3261 and PARITY_PREPARED=1, 2, or 3");
  }

  int const middleCarryMode = args.value("MIDCARRY_FUSED", 0);
  if (middleCarryMode < 0 || middleCarryMode > 2) {
    throw std::runtime_error("MIDCARRY_FUSED must be 0, 1 (middle fusion), or 2 (split-carry diagnostic)");
  }
  if (useMiddleCarry || useSplitCarry) {
#if !CUDA_BACKEND
    throw std::runtime_error("MIDCARRY_FUSED is CUDA-only");
#endif
    if (fft.shape.fft_type != FFT3161 || WIDTH != 512 ||
        fft.shape.middle != 8 || SMALL_H != 512 || !in_place ||
        useLongCarry || parityCorrect) {
      throw std::runtime_error(
        "MIDCARRY_FUSED requires in-place FFT3161 512:8:512, short carry, and no parity mode");
    }
    log(useMiddleCarry ? "Using two-phase middle/carry fusion\n" :
                         "Using two-phase split-carry diagnostic\n");
  }

  // Exact detached M31 width edge: the M31 inverse width runs as fftWGF31
  // before carryFused and the forward width runs as fftWOut31 after it, both
  // scheduled inside the M31 bottom half so they hide under the longer M61
  // chain.  The FFT323161 spelling of DETACHED_M31_EDGE remains the recorded
  // correctness-disabled timing scaffold; only FFT3161 has the exact producer
  // and consumer.
  detachedM31 = fft.shape.fft_type == FFT3161 ? args.value("DETACHED_M31_EDGE", 0) : 0;
  if (detachedM31) {
    if (detachedM31 != 1 && detachedM31 != 2) {
      throw std::runtime_error("DETACHED_M31_EDGE on FFT3161 must be 1 (both widths) or 2 (inverse only)");
    }
    if (!in_place || useLongCarry || useMiddleCarry || useSplitCarry ||
        args.value("GOLD_PAIR", 0) || args.value("L2_STRIPING", 0) ||
        args.value("GRAPHS", 0)) {
      throw std::runtime_error(
        "DETACHED_M31_EDGE on FFT3161 requires in-place, short fused carry, "
        "no GOLD_PAIR, no L2_STRIPING, and no GRAPHS");
    }
    log(detachedM31 == 1 ? "Using exact detached M31 width edge\n"
                         : "Using exact detached M31 inverse width\n");
  }

  if (fft.FFT_FP64 || fft.FFT_FP32) {
    kfftMidIn.setFixedArgs(3, bufTrigM);
    kfftHin.setFixedArgs(3, bufTrigH);
    ktailSquareZero.setFixedArgs(2, bufTrigH);
    ktailSquare.setFixedArgs(3, bufTrigH);
    ktailMulLow.setFixedArgs(4, bufTrigH);
    ktailMul.setFixedArgs(4, bufTrigH);
    kfftMidOut.setFixedArgs(3, bufTrigM);
    kfftW.setFixedArgs(2, bufTrigW);
  }

  fused31 = args.value("FUSED31", 0);

  if (fft.NTT_GF31) {
    kFusedMidTail31.setFixedArgs(3, bufTrigM, bufTrigH);
    kfftMidInGF31.setFixedArgs(3, bufTrigM);
    kfftHinGF31.setFixedArgs(3, bufTrigH);
    ktailSquareZeroGF31.setFixedArgs(2, bufTrigH);
    ktailSquareGF31.setFixedArgs(3, bufTrigH);
    ktailMulLowGF31.setFixedArgs(4, bufTrigH);
    ktailMulGF31.setFixedArgs(4, bufTrigH);
    kfftMidOutGF31.setFixedArgs(3, bufTrigM);
    kfftWGF31.setFixedArgs(2, bufTrigW);
    if (detachedM31) kfftWOut31.setFixedArgs(2, bufTrigW);
  }

  if (fft.NTT_GF61) {
    kfftMidInGF61.setFixedArgs(3, bufTrigM);
    kfftHinGF61.setFixedArgs(3, bufTrigH);
    ktailSquareZeroGF61.setFixedArgs(2, bufTrigH);
    ktailSquareGF61.setFixedArgs(3, bufTrigH);
    ktailMulLowGF61.setFixedArgs(4, bufTrigH);
    ktailMulGF61.setFixedArgs(4, bufTrigH);
    kfftMidOutGF61.setFixedArgs(3, bufTrigM);
    kfftWGF61.setFixedArgs(2, bufTrigW);
  }

  if (fft.NTT_RIESEL) {
    for (Kernel* k : {&kfftMidInR0, &kfftMidInR1}) k->setFixedArgs(3, bufTrigM);
    for (Kernel* k : {&kfftHinR0, &kfftHinR1}) k->setFixedArgs(3, bufTrigH);
    for (Kernel* k : {&ktailSquareZeroR0, &ktailSquareZeroR1}) k->setFixedArgs(2, bufTrigH);
    for (Kernel* k : {&ktailSquareR0, &ktailSquareR1}) k->setFixedArgs(3, bufTrigH);
    for (Kernel* k : {&ktailMulLowR0, &ktailMulLowR1, &ktailMulR0, &ktailMulR1}) k->setFixedArgs(4, bufTrigH);
    for (Kernel* k : {&kfftMidOutR0, &kfftMidOutR1}) k->setFixedArgs(3, bufTrigM);
    for (Kernel* k : {&kfftWR0, &kfftWR1}) k->setFixedArgs(2, bufTrigW);
    kfftPR0.setFixedArgs(2, bufTrigW, bufWeights);
    kfftPR1.setFixedArgs(2, bufTrigW, bufWeights);
  }

  if (fft.shape.fft_type == FFT31R2) {
    // M31's premultiplication is independent of the Riesel transform stream(s).
    kfftP.setFixedArgs(2, bufTrigW);
    kfftP31R2.setFixedArgs(2, bufTrigW);
    if (fft.NTT_RIESEL && fft.NTT_GF61) kfftPRPair.setFixedArgs(2, bufTrigW);
    else kfftPRPair.setFixedArgs(2, bufTrigW, bufWeights);
    for (Kernel* k : {&kCarryA, &kCarryAROE, &kCarryM, &kCarryMROE, &kCarryLL}) {
      k->setFixedArgs(3, bufCarry, bufWeights);
    }
    for (Kernel* k : {&kCarryA, &kCarryM, &kCarryLL}) { k->setFixedArgs(5, bufStatsCarry); }
    for (Kernel* k : {&kCarryAROE, &kCarryMROE})      { k->setFixedArgs(5, bufROE); }
    for (Kernel* k : {&kCarryFused, &kCarryFusedROE, &kCarryFusedMul, &kCarryFusedMulROE, &kCarryFusedLL}) {
      k->setFixedArgs(3, bufCarry, bufReady, bufTrigW, bufWeights);
    }
    for (Kernel* k : {&kCarryFusedROE, &kCarryFusedMulROE})           { k->setFixedArgs(7, bufROE); }
    for (Kernel* k : {&kCarryFused, &kCarryFusedMul, &kCarryFusedLL}) { k->setFixedArgs(7, bufStatsCarry); }
  } else if (fft.FFT_FP64 || fft.FFT_FP32) {                  // The FP versions take bufWeight arguments
    kfftP.setFixedArgs(2, bufTrigW, bufWeights);
    for (Kernel* k : {&kCarryA, &kCarryAROE, &kCarryM, &kCarryMROE, &kCarryLL}) { k->setFixedArgs(3, bufCarry, bufWeights); }
    for (Kernel* k : {&kCarryAParity, &kCarryAROEParity}) { k->setFixedArgs(3, bufCarry, bufWeights); }
    for (Kernel* k : {&kCarryA, &kCarryM, &kCarryLL}) { k->setFixedArgs(5, bufStatsCarry); }
    for (Kernel* k : {&kCarryAROE, &kCarryMROE})      { k->setFixedArgs(5, bufROE); }
    kCarryAParity.setFixedArgs(5, bufStatsCarry);
    kCarryAROEParity.setFixedArgs(5, bufROE);
    for (Kernel* k : {&kCarryFused, &kCarryFusedROE, &kCarryFusedMul, &kCarryFusedMulROE, &kCarryFusedLL}) {
      k->setFixedArgs(3, bufCarry, bufReady, bufTrigW, bufConstWeights, bufWeights);
    }
    for (Kernel* k : {&kCarryFusedROE, &kCarryFusedMulROE})           { k->setFixedArgs(8, bufROE); }
    for (Kernel* k : {&kCarryFused, &kCarryFusedMul, &kCarryFusedLL}) { k->setFixedArgs(8, bufStatsCarry); }
  } else {
    bool const useGoldWeights = fft.shape.fft_type == FFT3161 &&
                                args.value("GOLD_PAIR", 0);
    if (useGoldWeights) kfftP.setFixedArgs(2, bufTrigW, bufWeights);
    else kfftP.setFixedArgs(2, bufTrigW);
    for (Kernel* k : {&kCarryA, &kCarryAROE, &kCarryM, &kCarryMROE, &kCarryLL}) {
      if (useGoldWeights) k->setFixedArgs(3, bufCarry, bufWeights);
      else k->setFixedArgs(3, bufCarry);
    }
    for (Kernel* k : {&kCarryAParity, &kCarryAROEParity}) { k->setFixedArgs(3, bufCarry); }
    u32 const carryStatsArg = useGoldWeights ? 5 : 4;
    for (Kernel* k : {&kCarryA, &kCarryM, &kCarryLL}) { k->setFixedArgs(carryStatsArg, bufStatsCarry); }
    for (Kernel* k : {&kCarryAROE, &kCarryMROE})      { k->setFixedArgs(carryStatsArg, bufROE); }
    kCarryAParity.setFixedArgs(4, bufStatsCarry);
    kCarryAROEParity.setFixedArgs(4, bufROE);
    for (Kernel* k : {&kCarryFused, &kCarryFusedROE, &kCarryFusedMul, &kCarryFusedMulROE, &kCarryFusedLL}) {
      if (useGoldWeights) k->setFixedArgs(3, bufCarry, bufReady, bufTrigW, bufWeights);
      else k->setFixedArgs(3, bufCarry, bufReady, bufTrigW);
    }
    u32 const statsArg = useGoldWeights ? 7 : 6;
    for (Kernel* k : {&kCarryFusedROE, &kCarryFusedMulROE}) { k->setFixedArgs(statsArg, bufROE); }
    for (Kernel* k : {&kCarryFused, &kCarryFusedMul, &kCarryFusedLL}) { k->setFixedArgs(statsArg, bufStatsCarry); }
    if (detachedM31) {
      for (Kernel* k : {&kCarryFused, &kCarryFusedROE, &kCarryFusedMul, &kCarryFusedMulROE, &kCarryFusedLL}) {
        k->setFixedArgs(statsArg + 1, bufDetach31);
      }
    }
  }

  if (foldSyndrome) {
    // The folded sidecar follows the optional parity input/output arguments.
    for (Kernel* k : {&kCarryFused, &kCarryFusedROE, &kCarryFusedLL}) {
      u32 const foldArg = parityCorrect ? 11 : 9;
      k->setFixedArgs(foldArg, bufFolded);
      if (foldValidate) k->setFixedArgs(foldArg + 1, bufFoldWords);
    }
  }

  kCarryMiddleOut.setFixedArgs(4, bufTrigM, bufTrigW, bufStatsCarry);
  kCarryMiddleOutROE.setFixedArgs(4, bufTrigM, bufTrigW, bufROE);
  kCarryMiddleIn.setFixedArgs(3, bufTrigM, bufTrigW);
  kCarrySplitOut.setFixedArgs(4, bufTrigM, bufTrigW, bufStatsCarry);
  kCarrySplitOutROE.setFixedArgs(4, bufTrigM, bufTrigW, bufROE);
  kCarrySplitIn.setFixedArgs(3, bufTrigM, bufTrigW);

  carryB.setFixedArgs(1, bufCarry);

  kernIsEqual.setFixedArgs(2, bufTrue);

  bufReady.zero();
  bufROE.zero();
  bufStatsCarry.zero();
  bufTrue.write({1});

  if (args.verbose >= 99) {
    selftestTrig();
  }

  // Three-plane Good-Thomas transforms can schedule one independent field per
  // queue.  Other mixed transforms retain the historical two-queue setup.
  if (args.value("MULTI_Q", 0)) {
    auxQueues.push_back(Queue{*shared.context, args.profile, true});
    if (fft.NTT_RIESEL ||
        (args.value("GOOD_THOMAS3", 0) && fft.shape.fft_type == FFT323161))
      auxQueues.push_back(Queue{*shared.context, args.profile, true});
  }
  if (parityPrepared || parityPacked) {
    parityQueue = std::make_unique<Queue>(*shared.context, args.profile, true);
    parityPrepare.setQueue(parityQueue.get());
  }
  if (foldTransformEnabled) {
    foldTransform = std::make_unique<FoldTransform>(shared, profile, N);
  }

#if CUDA_BACKEND
  if (args.value("L2_PERSIST", 0)) {
    std::vector<cl_mem> const readOnlyBuffers{
      bufTrigH->get(), bufTrigM->get(), bufTrigW->get(),
      bufConstWeights.get(), bufWeights.get()
    };
    cudaSetL2Persistent(queue.get(), readOnlyBuffers);
    for (auto& auxQueue : auxQueues)
      cudaSetL2Persistent(auxQueue.get(), readOnlyBuffers);
  }
#endif

  // Set flag indicating we're going to use CUDA graphs
  use_graphs = graph_square[0].isSupported(shared.context->deviceId()) && args.value("GRAPHS", 1);

  // Set L1 cache configuration.  Really we should only do this once rather than once per worker.
  // However, the current way PRPLL is organized would then make this option hard to tune.
#if CUDA_BACKEND
    cudaSetL1Config(args.value("L1CUDA", 0));
#endif

  // Process the queue.  I don't know if this is really needed.
  queue.finish();
}

// Optionallly split some of the MiddleIn/Tail/MiddleOut kernels off of executing on the main queue to run on an auxiliary queue.
// This will increase GPU occupancy but will negatively impact L2 cache coherency.
// If the L2 cache is large enough so that all FFT data fits in the cache, this ought to be a win.
// If the L2 cache is small enough such that L2 cache hits are very low anyway, this might be a win.
void Gpu::splitQueue() {
  // Queue a sync event in the main queue.  Have all auxiliary queues wait on the event.
  EventHolder event = queue.createSyncEvent();
  for (auto & auxQueue : auxQueues) {
    auxQueue.waitForSyncEvent(&event);
  }
}

void Gpu::mergeQueue() {
  // Queue a sync event in each auxiliary queue(s).  Wait on the event(s) in the main queue.
  for (auto & auxQueue : auxQueues) {
    EventHolder event = auxQueue.createSyncEvent();
    queue.waitForSyncEvent(&event);
  }
}

// Convert the natural packed parity of the current integer words into one
// carry-ready bit per output coefficient.  It runs on its own queue so the
// permutation can overlap the much longer FP32/M61 bottom half.
void Gpu::prepareParity() {
  assert((parityPrepared || parityPacked) && parityQueue && !parityExpectedPending);
  EventHolder event = queue.createSyncEvent();
  parityQueue->waitForSyncEvent(&event);
  parityPrepare(*parityIn, bufParityExpected);
  parityExpectedPending = true;
}

void Gpu::waitParity() {
  if (!parityPrepared && !parityPacked) return;
  assert(parityExpectedPending && parityQueue);
  EventHolder event = parityQueue->createSyncEvent();
  queue.waitForSyncEvent(&event);
  parityExpectedPending = false;
}

// We've finished the "bottom half" of a squaring or multiply.  Replay the recorded bottom half kernel calls.
void Gpu::endBottomHalf() {
  replay();
  // Increment the squaring count.  The queue's squarings/multiplies count determines how long to sleep when the queue is full.
  // We only do this for the main command queue.  Auxiliary queues are not allowed to cause a CPU sleep.
  queue.incSquareCount();
}

// Replay the recorded bottom half kernels in a cache friendly order.  We support several options here using multiple openCl command queues.
void Gpu::replay() {
  // If there are no recorded kernels to replay, we're done
  if (recorded_kernels.size() == 0) return;

  // FUSED31 substitutes middle-in + tailSquare; mul bottom halves (tailMul)
  // must run the unfused path.
  replay_square_pass = false;
  {
    bool hasMidIn = false, hasTailSquare = false;
    for (auto kern : recorded_kernels) {
      if (kern == KMIDIN) hasMidIn = true;
      if (kern == KTAILSQUARE) hasTailSquare = true;
    }
    replay_square_pass = hasMidIn && hasTailSquare;   // fuse only complete midIn+tail replays
  }
  if (fused31) {
    static int logged = 0;
    if (logged < 6) {
      ++logged;
      std::string seq;
      for (auto kern : recorded_kernels) seq += std::to_string((int) kern) + " ";
      log("replay seq (square=%d): %s\n", (int) replay_square_pass, seq.c_str());
    }
  }

  // Get MULTI_Q and L2_STRIPING settings
  bool multi_q = args.value("MULTI_Q", 0);
  int l2_striping = args.value("L2_STRIPING", 0);
  auto inactiveCacheGroup = [&](int cache_group) {
    if (cache_group == 1) return !(fft.FFT_FP64 || fft.FFT_FP32);
    if (cache_group == 2) return !fft.NTT_GF31;
    if (cache_group == 3) return !fft.NTT_GF61;
    if (cache_group == 4)
      return !fft.NTT_RIESEL || args.value("GOOD_THOMAS3", 0) ||
             args.value("GOOD_THOMAS7", 0);
    if (cache_group == 5) return !fft.NTT_RIESEL;
    return true;
  };

  // In the simplest case, we use one command queue and process one data type at a time.  By processing one data type at a time, we reduce maximum L2 cache used.
  // For example, a 4M GF61+GF31 NTT needs just 32MB L2 cache during GF61 processing of fftMiddleIn, tailSquare, and fftMiddleOut (and only 16MB duing GF31 processing).
  // Without MULTI_Q, PRPLL needs 32MB + 16MB of L2 cache by processing both fftMiddleIns, then both tailSquares, then both fftMiddleOuts.

  if ((!multi_q || fft.shape.fft_type == FFT64 || fft.shape.fft_type == FFT61 || fft.shape.fft_type == FFT31 || fft.shape.fft_type == FFT32) && !l2_striping) {
    for (int cache_group = 1; cache_group <= NUM_CACHE_GROUPS; ++cache_group) {
      // Check for irrelevant cache group
      if (inactiveCacheGroup(cache_group)) continue;

      // Iterate over the recorded kernels.  Execute each.
      int arg = 0;
      for (auto kern : recorded_kernels) {
        replay_one(kern, cache_group, arg);
        arg = replay_next_arg(kern, arg);
      }
    }
  }

  // The next simple case, we use two command queues and process one data type in each queue.  This works well for large L2 caches where all FFT data fits in the cache.
  // The extra command queue can hide the latency in starting up kernels for each data type.  Also occupancy may benefit as the queue may be executing kernels with different
  // workgroup size, register usage, and local memory usage.  This case requires an FFT using at least two data types.

  else if (multi_q && !l2_striping) {
    // Switch tp using multiple command queues.
    splitQueue();
    for (int cache_group = 1; cache_group <= NUM_CACHE_GROUPS; ++cache_group) {
      // Check for irrelevant cache group
      if (inactiveCacheGroup(cache_group)) continue;

      // To better balance the load on the two command queues, put a 64-bit data type in one queue and two 32-bit data types in the other queue.
      Queue *q;
      if (cache_group == 1) q = &queue;
      if (cache_group == 2)
        q = (args.value("GOOD_THOMAS3", 0) &&
             fft.shape.fft_type == FFT323161) ? &auxQueues[1] :
            ((fft.shape.fft_type == FFT323161 ||
              fft.shape.fft_type == FFT3161 ||
              fft.shape.fft_type == FFT31R2) ? &queue : &auxQueues[0]);
      if (cache_group == 3) q = &auxQueues[0];
      if (cache_group == 4) q = (fft.NTT_RIESEL && fft.NTT_GF61) ? &queue : &auxQueues[0];
      if (cache_group == 5) q = &auxQueues[1];

      // Iterate over the recorded kernels.  Execute each.
      int arg = 0;
      for (auto kern : recorded_kernels) {
        replay_one(kern, cache_group, arg, q);
        arg = replay_next_arg(kern, arg);
      }
    }
    // Using multiple command queues, go back to a single command queue
    mergeQueue();
  }

  // The next case is L2 striping in one command queue.  The hope is two stripe groups plus one stripe are small enough to fit in the L2 caches
  // for the fftMiddleIn, tailSquare, and fftMiddleOut kernels.
  // Sadly, testing thusfar shows the extra overhead of more kernel launches and events/syncs outweighs the benefit of more L2 cache hits.
#ifdef ORIGINAL_VERSION         // Very readable, replaced by version below which merges the teo base_lo and base_hi kernel calls into one combined kernel call
  else if (!multi_q && l2_striping) {
    for (int cache_group = 1; cache_group <= NUM_CACHE_GROUPS; ++cache_group) {
      // Check for irrelevant cache group
      if (inactiveCacheGroup(cache_group)) continue;

      // Allow larger caches to do several L2 stripes in a single "stripe group" at a time (increases occupancy, reduces kernel launch costs).
      u32 stripe_group_size = l2_striping;

//For now, only support INPLACE with its 16x16 transpose      

      // Loop over all the stripes for this data type.  Let's define a stripe as one "column" of data processed by fftMiddleIn producing 16*MIDDLE tailSquare lines.
      // Since tailSquare operates on Hermetian pairs of lines, fftMiddleIn alternates operating on low stripes and high stripes.
      u32 num_stripes = fft.shape.width / 16;
      u32 num_stripe_groups = num_stripes / stripe_group_size;
      for (u32 i = 0; i < num_stripe_groups / 2; ++i) {
        // Base_lo refers to the starting fftMiddleIn column number (x coordinate).  When fftMiddleIn uses a 16x16 transpose, base_lo advances 16 at a time.
        // Base_hi refers to the fftMiddleIn column number that outputs the lines needed for Hermetian matching in tailSquare.
        u32 base_lo = i * stripe_group_size * 16;
        u32 base_hi = (num_stripe_groups - 1 - i) * stripe_group_size * 16;
        bool last_block = (base_hi == fft.shape.width / 2);

        // Iterate over the recorded kernels.  Execute each.
        int arg = 0;
        for (auto kern : recorded_kernels) {

          if (kern == KMIDIN) {
            u32 oneStripeKernelsToExecute = fft.shape.height / 16;
            // If MIDDLE is odd, the first fftMiddleIn call must be preceeded by a fftMiddleIn call to produce the special N/2 tailSquare line from the WIDTH/2 column.
            if (base_lo == 0 && (fft.shape.middle & 1)) {
              replay_one(kern, cache_group, arg, &queue, fft.shape.width / 2, oneStripeKernelsToExecute);
            }

            // Produce one stripe group.  The last base value must take into account that the N/2 stripe has already been done if MIDDLE is odd.
            u32 kernelsToExecute = stripe_group_size * oneStripeKernelsToExecute;
            replay_one(kern, cache_group, arg, &queue, base_lo, kernelsToExecute);

            u32 base = base_hi;                                                                                   // Read full base_hi stripe groups (usually)
            if (last_block && (fft.shape.middle & 1)) base += 16, kernelsToExecute -= oneStripeKernelsToExecute;  // Skip first stripe for last block if MIDDLE is odd
            if (kernelsToExecute) replay_one(kern, cache_group, arg, &queue, base, kernelsToExecute);
          }

          else if (kern == KFFTHIN) {
            // fftMiddleIn produces 16 * MIDDLE lines
            replay_one(kern, cache_group, arg, &queue, base_lo, stripe_group_size * 16 * fft.shape.middle);
            replay_one(kern, cache_group, arg, &queue, base_hi, stripe_group_size * 16 * fft.shape.middle);
          }
 
          else if (kern == KTAILSQUARE || kern == KTAILMUL || kern == KTAILMULLOW) {
            // fftMiddleIn produces 16 * MIDDLE lines in base_lo and base_hi.  tailSquare kernel processes 2 lines linked by Hermetian symmetry.
            // Tail kernels use two dimensions, the X coordinate bumps line number by one, the Y coordinate bumps line number by WIDTH.
            u32 kernelsToExecuteX = stripe_group_size * 16;
            u32 kernelsToExecuteY = (fft.shape.middle + 1) / 2;  // For base_lo, round odd middles up.

            // We can now completely process lines from the lower half of the base_lo stripe group (Hermetian mates are mostly in upper half of base_hi stripe group)
            replay_one(kern, cache_group, arg, &queue, base_lo, kernelsToExecuteX, kernelsToExecuteY);

            // We can process most lines from the lower half of the base_hi stripe group (Hermetian mates are mostly in upper half of base_lo stripe group)
            // The first line in the stripe group is the only line that is not ready for base_hi tail processing.
            u32 base = base_hi + 1;                                           // Skip first line in base_hi (usually)
            if (base_lo == 0) kernelsToExecuteX--;                            // Do one fewer line for the first tail call.
            if (base_hi == fft.shape.width / 2) base--, kernelsToExecuteX++;  // Last tail call does not skip first line
            if (fft.shape.middle & 1) kernelsToExecuteY--;                    // For base_hi, round odd middles down.
            replay_one(kern, cache_group, arg, &queue, base, kernelsToExecuteX, kernelsToExecuteY);
          }

          else if (kern == KMIDOUT) {
            u32 oneStripeKernelsToExecute = fft.shape.height / 16;
            u32 kernelsToExecute = stripe_group_size * oneStripeKernelsToExecute;

            // We've completely processed lines from the base_lo stripe group
            replay_one(kern, cache_group, arg, &queue, base_lo, kernelsToExecute);

            // The first stripe in the base_hi stripe group is not ready for output.
            u32 base = base_hi + 16;                                                                        // Skip first stripe in base_hi (usually)
            if (base_lo == 0) kernelsToExecute -= oneStripeKernelsToExecute;                                // Do one fewer stripe for the first midOut call.
            if (base_hi == fft.shape.width / 2) base -= 16, kernelsToExecute += oneStripeKernelsToExecute;  // Last midOut call does not skip first stripe
            if (kernelsToExecute) replay_one(kern, cache_group, arg, &queue, base, kernelsToExecute);
          }

          // Skip other kernels (KFFTW)
          else;

          // Advance argument index
          arg = replay_next_arg(kern, arg);
        }
      }

      // Iterate over the recorded kernels.  Execute any not already executed (KFFTW).
      // FFTW cannot benefit from L2 striping, it can only benefit from datatype grouping.
      int arg = 0;
      for (auto kern : recorded_kernels) {
        if (kern == KFFTW) replay_one(kern, cache_group, arg, &queue);
        arg = replay_next_arg(kern, arg);
      }
    }
  }
#endif

  // The next case is L2 striping in one command queue.  The hope is two stripe groups plus one stripe are small enough to fit in the L2 caches
  // for the fftMiddleIn, tailSquare, and fftMiddleOut kernels.
  // Sadly, testing thusfar shows the extra overhead of more kernel launches and events/syncs outweighs the benefit of more L2 cache hits.

  else if (!multi_q && l2_striping) {
    for (int cache_group = 1; cache_group <= NUM_CACHE_GROUPS; ++cache_group) {
      // Check for irrelevant cache group
      if (inactiveCacheGroup(cache_group)) continue;

      // Allow larger caches to do several L2 stripes in a single "stripe group" at a time (increases occupancy, reduces kernel launch costs).
      u32 stripe_group_size = l2_striping;

//For now, only support INPLACE with its 16x16 transpose      

      // Loop over all the stripes for this data type.  Let's define a stripe as one "column" of data processed by fftMiddleIn producing 16*MIDDLE tailSquare lines.
      // Since tailSquare operates on Hermetian pairs of lines, fftMiddleIn alternates operating on low stripes and high stripes.
      u32 num_stripes = fft.shape.width / 16;
      u32 num_stripe_groups = num_stripes / stripe_group_size;
      for (u32 i = 0; i < num_stripe_groups / 2; ++i) {
        bool last_block = (i == num_stripe_groups / 2 - 1);
        // Base_lo refers to the starting fftMiddleIn column number (x coordinate).  When fftMiddleIn uses a 16x16 transpose, base_lo advances 16 at a time.
        // Base_hi refers to the fftMiddleIn column number that outputs the lines needed for Hermetian matching in tailSquare.
        u32 base_lo = i * stripe_group_size * 16;
        // u32 base_hi = (num_stripe_groups - 1 - i) * stripe_group_size * 16;

        // Iterate over the recorded kernels.  Execute each.
        int arg = 0;
        for (auto kern : recorded_kernels) {

          if (kern == KMIDIN) {
            // Do midIn on base_lo and base_hi.  If MIDDLE is odd, the first midIn call must be preceeded by a midIn call to produce the special N/2 tailSquare
            // line from the WIDTH/2 column.  The last base pair must take into account that the N/2 stripe has already been done if MIDDLE is odd.
            u32 oneStripeKernelsToExecute = fft.shape.height / 16;
            u32 kernelsToExecute = 2 * stripe_group_size * oneStripeKernelsToExecute;
            if (base_lo == 0 && (fft.shape.middle & 1)) kernelsToExecute += oneStripeKernelsToExecute;
            if (last_block && (fft.shape.middle & 1)) kernelsToExecute -= oneStripeKernelsToExecute;
            replay_one(kern, cache_group, arg, &queue, base_lo, kernelsToExecute);
          }

          else if (kern == KFFTHIN) {
            // fftMiddleIn produces 2 * stripe_group_size * 16 * MIDDLE lines
            replay_one(kern, cache_group, arg, &queue, base_lo, 2 * stripe_group_size * 16 * fft.shape.middle);
          }
 
          else if (kern == KTAILSQUARE || kern == KTAILMUL || kern == KTAILMULLOW) {
            // We can now completely process lines from the lower half of the base_lo stripe group (Hermetian mates are mostly in upper half of base_hi stripe group)
            u32 half_size = (fft.shape.middle + 1) / 2;  // For base_lo, round odd middles up.
            u32 stripeGroupLines = stripe_group_size * 16;
            u32 kernelsToExecute = half_size * stripeGroupLines;
            // We can process most lines from the lower half of the base_hi stripe group (Hermetian mates are mostly in upper half of base_lo stripe group)
            // The first line in the stripe group is the only line that is not ready for base_hi tail processing.
            if (base_lo == 0) stripeGroupLines--;               // Do one fewer line for the first tail call.
            if (last_block) stripeGroupLines++;                 // Last tail call does not skip first line
            if (fft.shape.middle & 1) half_size--;              // For base_hi, round odd middles down.
            kernelsToExecute += half_size * stripeGroupLines;
            // Call the kernel
            replay_one(kern, cache_group, arg, &queue, base_lo, kernelsToExecute);
          }

          else if (kern == KMIDOUT) {
            u32 oneStripeKernelsToExecute = fft.shape.height / 16;
            u32 kernelsToExecute = 2 * stripe_group_size * oneStripeKernelsToExecute;
            // We've completely processed lines from the base_lo stripe group.
            // The first stripe in the base_hi stripe group is not ready for output.
            // The last stripe group will does an extra stripe.
            if (base_lo == 0) kernelsToExecute -= oneStripeKernelsToExecute;  // Do one fewer stripe for the first midOut call
            if (last_block) kernelsToExecute += oneStripeKernelsToExecute;    // Last midOut call does not skip first stripe
            replay_one(kern, cache_group, arg, &queue, base_lo, kernelsToExecute);
          }

          // Skip other kernels (KFFTW)
          else {}

          // Advance argument index
          arg = replay_next_arg(kern, arg);
        }
      }

      // Iterate over the recorded kernels.  Execute any not already executed (KFFTW).
      // FFTW cannot benefit from L2 striping, it can only benefit from datatype grouping.
      int arg = 0;
      for (auto kern : recorded_kernels) {
        if (kern == KFFTW) replay_one(kern, cache_group, arg, &queue);
        arg = replay_next_arg(kern, arg);
      }
    }
  }

  // The last case is L2 striping in two command queues.  The hope is four stripe groups plus two stripes are small enough to fit in the L2 caches for the
  // fftMiddleIn, tailSquare, and fftMiddleOut kernels.  The hope also is the dual queue approach hides the overhead introduced by more kernel launches.
  // Sadly, testing thusfar shows the extra overhead of more kernel launches and events/syncs outweighs the benefit of more L2 cache hits.

  else if (multi_q && l2_striping) {
    Queue *queues[2] = {&queue, &auxQueues[0]};
    EventHolder midInEvents[2];
    EventHolder tailEvents[2];

    splitQueue();
    for (int cache_group = 1; cache_group <= NUM_CACHE_GROUPS; ++cache_group) {
      bool has_unexecuted_kernels = false;

      // Check for irrelevant cache group
      if (inactiveCacheGroup(cache_group)) continue;

      // Allow larger caches to do several L2 stripes at a time (increases occupancy, reduces kernel launch costs).
      u32 stripe_group_size = l2_striping;

//For now, only support INPLACE with its 16x16 transpose      

      // Loop over all the stripes for this data type.  Let's define a stripe as one "column" of data processed by fftMiddleIn producing 16*MIDDLE tailSquare lines.
      // Since tailSquare operates on Hermetian pairs of lines, fftMiddleIn alternates operating on low stripes and high stripes.
      u32 num_stripes = fft.shape.width / 16;
      u32 num_stripe_groups = num_stripes / stripe_group_size;
      for (u32 i = 0; i < num_stripe_groups / 2; ++i) {
        int q = (i & 1);               // Index into which command queue to use
        u32 i_within_queue = i >> 1;
        bool last_i_within_queue = (i_within_queue == num_stripe_groups / 4 - 1);
        // Base_lo refers to the starting fftMiddleIn column number (x coordinate).  When fftMiddleIn uses a 16x16 transpose, base_lo advances 16 at a time.
        // Base_hi refers to the fftMiddleIn column number that outputs the lines needed for Hermetian matching in tailSquare.
        u32 base_lo = (i_within_queue + q * num_stripe_groups / 4) * stripe_group_size * 16;
        //u32 base_hi = fft.shape.width - stripe_group_size * 16 - base_lo;

        // Iterate over the recorded kernels.  Execute each.
        int arg = 0;
        for (auto kern : recorded_kernels) {

          if (kern == KMIDIN) {
            // Do midIn on base_lo and base_hi.  If MIDDLE is odd, the first midIn call must be preceeded by a midIn call to produce the special N/2 tailSquare
            // line from the WIDTH/2 column.  The first midIn call in the second command queue also must be preceeded by a midIn call to produce one L2 stripe.
            // The last base_hi in each queue must take into account pre-read stripes in the other queue.
#define q_requires_preread(q) (((q) == 0 && fft.shape.middle & 1) || ((q) == 1))
            u32 oneStripeKernelsToExecute = fft.shape.height / 16;
            u32 kernelsToExecute = 2 * stripe_group_size * oneStripeKernelsToExecute;
            if (i_within_queue == 0 && q_requires_preread(q)) kernelsToExecute += oneStripeKernelsToExecute;
            if (last_i_within_queue && q_requires_preread(!q)) kernelsToExecute -= oneStripeKernelsToExecute;
            replay_one(kern, cache_group, arg, queues[q], base_lo, kernelsToExecute);
            if (i_within_queue == 0 && q_requires_preread(q)) midInEvents[q] = queues[q]->createSyncEvent();
            // The last block must sync with a pre-read from the other command queue
            // BUG - if block is both i_within_queue == 0 and last_block_in_queue, then midInEvents[1] does not exist!
            // Sanity checking L2_STRIPING setting to a max of WIDTH/128 at startup eliminates this bug.
            if (last_i_within_queue && q_requires_preread(!q)) queues[q]->waitForSyncEvent(&midInEvents[!q]);
          }

          else if (kern == KFFTHIN) {
            // fftMiddleIn produces 2 * stripe_group_size * 16 * MIDDLE lines
            replay_one(kern, cache_group, arg, queues[q], base_lo, 2 * stripe_group_size * 16 * fft.shape.middle);
          }

          else if (kern == KTAILSQUARE || kern == KTAILMUL || kern == KTAILMULLOW) {
            // We can now completely process lines from the lower half of the base_lo stripe group (Hermetian mates are mostly in upper half of base_hi stripe group)
            u32 half_size = (fft.shape.middle + 1) / 2;    // For base_lo, round odd middles up.
            u32 stripeGroupLines = stripe_group_size * 16;
            u32 kernelsToExecute = half_size * stripeGroupLines;
            // We can process most lines from the lower half of the base_hi stripe group (Hermetian mates are mostly in upper half of base_lo stripe group)
            // The first line in the stripe group is the only line that is not ready for base_hi tail processing.
            if (i == 0) stripeGroupLines--;                         // Do one fewer line for the first tail call.
            if (last_i_within_queue && q == 1) stripeGroupLines++;  // Last tail call does not skip first line
            half_size = fft.shape.middle / 2;                       // For base_hi, round odd middles down.
            kernelsToExecute += half_size * stripeGroupLines;
            // Call the kernel
            replay_one(kern, cache_group, arg, queues[q], base_lo, kernelsToExecute);
            if (i_within_queue == 0 && q_requires_preread(q)) tailEvents[q] = queues[q]->createSyncEvent();
// We could eliminate two events and syncs by having the last midIn wait on the first tailsquare.  It exposes a little less parallellism, but perhaps that is irrelevant.
// For even middles, we only save one event and sync.
            // The last block must sync with the first tailSquare in the other command queue
            // BUG - if block is both i_within_queue == 0 and last_block_in_queue, then tailEvents[1] does not exist!
            // Sanity checking L2_STRIPING setting to a max of WIDTH/128 at startup eliminates this bug.
            if (last_i_within_queue && q_requires_preread(!q)) queues[q]->waitForSyncEvent(&tailEvents[!q]);
          }

          else if (kern == KMIDOUT) {
            u32 oneStripeKernelsToExecute = fft.shape.height / 16;
            u32 kernelsToExecute = 2 * stripe_group_size * oneStripeKernelsToExecute;
            // We've completely processed lines from the base_lo stripe group.
            // The first stripe in the base_hi stripe group is not ready for output.
            // The last stripe group does an extra stripe.
            if (i_within_queue == 0) kernelsToExecute -= oneStripeKernelsToExecute;  // Do one fewer stripe for the first midOut call
            if (last_i_within_queue) kernelsToExecute += oneStripeKernelsToExecute;  // Last midOut call does not skip first stripe
            replay_one(kern, cache_group, arg, queues[q], base_lo, kernelsToExecute);
          }

          // Skip other kernels (KFFTW)
          else
            has_unexecuted_kernels = true;

          // Advance argument index
          arg = replay_next_arg(kern, arg);
        }
      }

      // Iterate over the recorded kernels again.  Execute any not already executed (KFFTW).
      // FFTW cannot benefit from L2 striping, it can only benefit from datatype grouping.
      if (has_unexecuted_kernels) {
        int arg = 0;
        mergeQueue();
        for (auto kern : recorded_kernels) {
          if (kern == KFFTW) replay_one(kern, cache_group, arg, &queue);
          arg = replay_next_arg(kern, arg);
        }
        splitQueue();
      }
    }
    mergeQueue();
  }

  // Empty the recorded kernels queue
  recorded_kernels.clear();
  recorded_kernel_args.clear();
}

// Replay one recorded kernel on the specified queue, with specified base and kernelsToExecute.  Two dimensional work groups are supported for some kernels.
// A kernelsToExecuteX of zero is permitted - used for the default kernelsToExecute (a.k.a. workSize) set at kernel creation that operates on all the FFT data.
void Gpu::replay_one(enum BOTTOM_HALF_KERNELS kern, int cache_group, int arg, Queue *q, int base, int kernelsToExecuteX, int kernelsToExecuteY) {

  // Call the appropriate kernel
  if (kern == KMIDIN) {
    Buffer<double> const *buf = recorded_kernel_args[arg++];
    // If not in place, the input is from the scratch buffer
    Buffer<double> const *in = in_place ? buf : &buf3;
    Buffer<double> const *out = buf;
    if (cache_group == 1) { kfftMidIn.setQueue(q); kfftMidIn.setKernelsToExecute(kernelsToExecuteX); kfftMidIn(*out, *in, base); }
    if (cache_group == 2) {
      if (fused31 && replay_square_pass) { kFusedMidTail31.setQueue(q); kFusedMidTail31(buf3, *in, base); }  // tail output -> scratch (in-place hazard)
      else { kfftMidInGF31.setQueue(q); kfftMidInGF31.setKernelsToExecute(kernelsToExecuteX); kfftMidInGF31(*out, *in, base); }
    }
    if (cache_group == 3) { kfftMidInGF61.setQueue(q); kfftMidInGF61.setKernelsToExecute(kernelsToExecuteX); kfftMidInGF61(*out, *in, base); }
    if (cache_group == 4) { kfftMidInR0.setQueue(q); kfftMidInR0.setKernelsToExecute(kernelsToExecuteX); kfftMidInR0(*out, *in, base); }
    if (cache_group == 5) { kfftMidInR1.setQueue(q); kfftMidInR1.setKernelsToExecute(kernelsToExecuteX); kfftMidInR1(*out, *in, base); }
  }

  if (kern == KFFTHIN) {
    Buffer<double> const *out = recorded_kernel_args[arg++];
    Buffer<double> const *in = recorded_kernel_args[arg++];
    if (cache_group == 1) { kfftHin.setQueue(q); kfftHin.setKernelsToExecute(kernelsToExecuteX); kfftHin(*out, *in, base); }
    if (cache_group == 2) { kfftHinGF31.setQueue(q); kfftHinGF31.setKernelsToExecute(kernelsToExecuteX); kfftHinGF31(*out, *((fused31 && replay_square_pass) ? &buf3 : in), base); }
    if (cache_group == 3) { kfftHinGF61.setQueue(q); kfftHinGF61.setKernelsToExecute(kernelsToExecuteX); kfftHinGF61(*out, *in, base); }
    if (cache_group == 4) { kfftHinR0.setQueue(q); kfftHinR0.setKernelsToExecute(kernelsToExecuteX); kfftHinR0(*out, *in, base); }
    if (cache_group == 5) { kfftHinR1.setQueue(q); kfftHinR1.setKernelsToExecute(kernelsToExecuteX); kfftHinR1(*out, *in, base); }
  }

  if (kern == KTAILSQUARE) {
    Buffer<double> const *buf = recorded_kernel_args[arg++];
    // If not in place, the output is to the scratch buffer
    Buffer<double> const *in = buf;
    Buffer<double> const *out = in_place ? buf : &buf3;
    if (!tail_single_kernel && base == 0) {
      if (cache_group == 1) { ktailSquareZero.setQueue(q); ktailSquareZero(*out, *in); }
      if (cache_group == 2 && !(fused31 && replay_square_pass)) { ktailSquareZeroGF31.setQueue(q); ktailSquareZeroGF31(*out, *in); }
      if (cache_group == 3) { ktailSquareZeroGF61.setQueue(q); ktailSquareZeroGF61(*out, *in); }
      if (cache_group == 4) { ktailSquareZeroR0.setQueue(q); ktailSquareZeroR0(*out, *in); }
      if (cache_group == 5) { ktailSquareZeroR1.setQueue(q); ktailSquareZeroR1(*out, *in); }
      if (kernelsToExecuteX) kernelsToExecuteX--;
    }
    if (cache_group == 1) { ktailSquare.setQueue(q); ktailSquare.setKernelsToExecute(kernelsToExecuteX, kernelsToExecuteY); ktailSquare(*out, *in, base); }
    if (cache_group == 2 && !(fused31 == 1 && replay_square_pass)) { ktailSquareGF31.setQueue(q); ktailSquareGF31.setKernelsToExecute(kernelsToExecuteX, kernelsToExecuteY); ktailSquareGF31(*out, *((fused31 == 2 && replay_square_pass) ? &buf3 : in), base); }
    if (cache_group == 3) { ktailSquareGF61.setQueue(q); ktailSquareGF61.setKernelsToExecute(kernelsToExecuteX, kernelsToExecuteY); ktailSquareGF61(*out, *in, base); }
    if (cache_group == 4) { ktailSquareR0.setQueue(q); ktailSquareR0.setKernelsToExecute(kernelsToExecuteX, kernelsToExecuteY); ktailSquareR0(*out, *in, base); }
    if (cache_group == 5) { ktailSquareR1.setQueue(q); ktailSquareR1.setKernelsToExecute(kernelsToExecuteX, kernelsToExecuteY); ktailSquareR1(*out, *in, base); }
  }

  if (kern == KTAILMUL) {
    Buffer<double> const *buf = recorded_kernel_args[arg++];
    Buffer<double> const *in2 = recorded_kernel_args[arg++];
    // If not in place, the output is to the scratch buffer
    Buffer<double> const *in1 = buf;
    Buffer<double> const *out = in_place ? buf : &buf3;
    if (cache_group == 1) { ktailMul.setQueue(q); ktailMul.setKernelsToExecute(kernelsToExecuteX, kernelsToExecuteY); ktailMul(*out, *in1, *in2, base); }
    if (cache_group == 2) { ktailMulGF31.setQueue(q); ktailMulGF31.setKernelsToExecute(kernelsToExecuteX, kernelsToExecuteY); ktailMulGF31(*out, *in1, *in2, base); }
    if (cache_group == 3) { ktailMulGF61.setQueue(q); ktailMulGF61.setKernelsToExecute(kernelsToExecuteX, kernelsToExecuteY); ktailMulGF61(*out, *in1, *in2, base); }
    if (cache_group == 4) { ktailMulR0.setQueue(q); ktailMulR0.setKernelsToExecute(kernelsToExecuteX, kernelsToExecuteY); ktailMulR0(*out, *in1, *in2, base); }
    if (cache_group == 5) { ktailMulR1.setQueue(q); ktailMulR1.setKernelsToExecute(kernelsToExecuteX, kernelsToExecuteY); ktailMulR1(*out, *in1, *in2, base); }
  }

  if (kern == KTAILMULLOW) {
    Buffer<double> const *buf = recorded_kernel_args[arg++];
    Buffer<double> const *in2 = recorded_kernel_args[arg++];
    // If not in place, the output is to the scratch buffer
    Buffer<double> const *in1 = buf;
    Buffer<double> const *out = in_place ? buf : &buf3;
    if (cache_group == 1) { ktailMulLow.setQueue(q); ktailMulLow.setKernelsToExecute(kernelsToExecuteX, kernelsToExecuteY); ktailMulLow(*out, *in1, *in2, base); }
    if (cache_group == 2) { ktailMulLowGF31.setQueue(q); ktailMulLowGF31.setKernelsToExecute(kernelsToExecuteX, kernelsToExecuteY); ktailMulLowGF31(*out, *in1, *in2, base); }
    if (cache_group == 3) { ktailMulLowGF61.setQueue(q); ktailMulLowGF61.setKernelsToExecute(kernelsToExecuteX, kernelsToExecuteY); ktailMulLowGF61(*out, *in1, *in2, base); }
    if (cache_group == 4) { ktailMulLowR0.setQueue(q); ktailMulLowR0.setKernelsToExecute(kernelsToExecuteX, kernelsToExecuteY); ktailMulLowR0(*out, *in1, *in2, base); }
    if (cache_group == 5) { ktailMulLowR1.setQueue(q); ktailMulLowR1.setKernelsToExecute(kernelsToExecuteX, kernelsToExecuteY); ktailMulLowR1(*out, *in1, *in2, base); }
  }

  if (kern == KMIDOUT) {
    Buffer<double> const *buf = recorded_kernel_args[arg++];
    // If not in place, the input is from the scratch buffer
    Buffer<double> const *in = in_place ? buf : &buf3;
    Buffer<double> const *out = buf;
    if (cache_group == 1) { kfftMidOut.setQueue(q); kfftMidOut.setKernelsToExecute(kernelsToExecuteX); kfftMidOut(*out, *in, base); }
    if (cache_group == 2) { kfftMidOutGF31.setQueue(q); kfftMidOutGF31.setKernelsToExecute(kernelsToExecuteX); kfftMidOutGF31(*out, *((fused31 == 1 && replay_square_pass) ? &buf3 : in), base); }  // FUSED31 tail wrote to scratch
    if (cache_group == 3) { kfftMidOutGF61.setQueue(q); kfftMidOutGF61.setKernelsToExecute(kernelsToExecuteX); kfftMidOutGF61(*out, *in, base); }
    if (cache_group == 4) { kfftMidOutR0.setQueue(q); kfftMidOutR0.setKernelsToExecute(kernelsToExecuteX); kfftMidOutR0(*out, *in, base); }
    if (cache_group == 5) { kfftMidOutR1.setQueue(q); kfftMidOutR1.setKernelsToExecute(kernelsToExecuteX); kfftMidOutR1(*out, *in, base); }
  }

  if (kern == KFFTW) {
    Buffer<double> const *out = recorded_kernel_args[arg++];
    Buffer<double> const *in = recorded_kernel_args[arg++];
    if (cache_group == 1) { kfftW.setQueue(q); kfftW(*out, *in); }
    if (cache_group == 2) { kfftWGF31.setQueue(q); kfftWGF31(*out, *in); }
    if (cache_group == 3) { kfftWGF61.setQueue(q); kfftWGF61(*out, *in); }
    if (cache_group == 4) { kfftWR0.setQueue(q); kfftWR0(*out, *in); }
    if (cache_group == 5) { kfftWR1.setQueue(q); kfftWR1(*out, *in); }
  }

  // Detached M31 width edge.  Both kernels belong to the GF31 cache group, so
  // in MULTI_Q mode they land on the main queue where the M31 bottom half runs:
  // the inverse width (KDETACHA) executes after fftMiddleOutGF31 and overlaps
  // the longer M61 chain; the forward width (KDETACHB) executes right after the
  // fused carry that produced its flat input, before the next fftMiddleInGF31.
  if (kern == KDETACHA) {
    Buffer<double> const *in = recorded_kernel_args[arg++];
    if (cache_group == 2) { kfftWGF31.setQueue(q); kfftWGF31(bufDetach31, *in); }
  }

  if (kern == KDETACHB) {
    Buffer<double> const *out = recorded_kernel_args[arg++];
    if (cache_group == 2) { kfftWOut31.setQueue(q); kfftWOut31(*out, bufDetach31); }
  }
}

// Advance the index into the array of kernel arguments
int Gpu::replay_next_arg(enum BOTTOM_HALF_KERNELS kern, int arg) {

  if (kern == KMIDIN || kern == KTAILSQUARE || kern == KMIDOUT || kern == KDETACHA || kern == KDETACHB) {
    return arg + 1;
  }

  else { //if (kern == KFFTHIN || kern == KTAILMUL || kern == KTAILMULLOW || kern == KFFTW) {
    return arg + 2;
  }
}

// Call the appropriate kernels to support hybrid FFTs and NTTs

void Gpu::fftP(Buffer<double>& buf, Buffer<Word>& in) {
  // Work around a troublesome oddball case.  ModMul calls fftP and fftMidIn on one multiplication argument.  If !in_place, fftP writes to buf3, and fftMidIn is queued.
  // Modmul then calls fftP on the other multiplication argument.  If we don't replay now, fftP overwrite buf3.
  replay();  
  // If not in place, instead write the output to the scratch buffer
  Buffer<double>  const*out = in_place ? &buf : &buf3;
  bool const packedRiesel = fft.shape.fft_type == FFT31R2 && fft.NTT_GF61;
  bool const fourRiesel = packedRiesel && fft.NTT_RIESEL;
  bool const goodThomas = fourRiesel &&
                          (args.value("GOOD_THOMAS3", 0) ||
                           args.value("GOOD_THOMAS7", 0));
  if (goodThomas && args.value("MULTI_Q", 0)) {
    splitQueue();
    kfftP.setQueue(&queue);
    kfftPRPair.setQueue(&auxQueues[0]);
    kfftPR1.setQueue(&auxQueues[1]);
    kfftP(*out, in);
    kfftPRPair(*out, in);
    kfftPR1(*out, in);
    mergeQueue();
  } else if (goodThomas) {
    kfftP(*out, in);
    kfftPRPair(*out, in);
    kfftPR1(*out, in);
  } else if (fourRiesel && args.value("MULTI_Q", 0)) {
    splitQueue();
    kfftP.setQueue(&queue);
    kfftPR0.setQueue(&queue);
    kfftPRPair.setQueue(&auxQueues[0]);
    kfftPR1.setQueue(&auxQueues[1]);
    kfftP(*out, in);
    kfftPR0(*out, in);
    kfftPRPair(*out, in);
    kfftPR1(*out, in);
    mergeQueue();
  } else if (fourRiesel) {
    kfftP(*out, in);
    kfftPRPair(*out, in);
    kfftPR0(*out, in);
    kfftPR1(*out, in);
  } else if (packedRiesel && args.value("MULTI_Q", 0)) {
    splitQueue();
    kfftP31R2.setQueue(&queue);
    kfftPRPair.setQueue(&auxQueues[0]);
    kfftP31R2(*out, in);
    kfftPRPair(*out, in);
    mergeQueue();
  } else if (packedRiesel) {
    kfftP31R2(*out, in);
    kfftPRPair(*out, in);
  } else if (fft.NTT_RIESEL && args.value("MULTI_Q", 0)) {
    splitQueue();
    kfftP.setQueue(&queue);
    kfftPR0.setQueue(&auxQueues[0]);
    kfftPR1.setQueue(&auxQueues[1]);
    kfftP(*out, in);
    kfftPR0(*out, in);
    kfftPR1(*out, in);
    mergeQueue();
  } else {
    kfftP(*out, in);
    if (fft.NTT_RIESEL) {
      kfftPR0(*out, in);
      kfftPR1(*out, in);
    }
  }
}

void Gpu::fftMidIn(Buffer<double>& buf) {
  // Record this call for later playback
  recorded_kernels.push_back(KMIDIN);
  recorded_kernel_args.push_back(&buf);
}

void Gpu::fftHin(Buffer<double>& out, Buffer<double>& in) {
  // Record this call for later playback
  recorded_kernels.push_back(KFFTHIN);
  recorded_kernel_args.push_back(&out);
  recorded_kernel_args.push_back(&in);
}

void Gpu::tailSquare(Buffer<double>& buf) {
  // Record this call for later playback
  recorded_kernels.push_back(KTAILSQUARE);
  recorded_kernel_args.push_back(&buf);
}

void Gpu::tailMul(Buffer<double>& buf, Buffer<double>& in2) {
  // Record this call for later playback
  recorded_kernels.push_back(KTAILMUL);
  recorded_kernel_args.push_back(&buf);
  recorded_kernel_args.push_back(&in2);
}

void Gpu::tailMulLow(Buffer<double>& buf, Buffer<double>& in2) {
  // Record this call for later playback
  recorded_kernels.push_back(KTAILMULLOW);
  recorded_kernel_args.push_back(&buf);
  recorded_kernel_args.push_back(&in2);
}

void Gpu::fftMidOut(Buffer<double>& buf) {
  // Record this call for later playback
  recorded_kernels.push_back(KMIDOUT);
  recorded_kernel_args.push_back(&buf);
}

void Gpu::fftW(Buffer<double>& out, Buffer<double>& in) {
  // Record this call for later playback
  recorded_kernels.push_back(KFFTW);
  recorded_kernel_args.push_back(&out);
  recorded_kernel_args.push_back(&in);
  // This kernel always ends the "bottom half".  Replay the recorded kernel calls.
  endBottomHalf();
}

void Gpu::carryA(Buffer<Word>& out, Buffer<double>& in) {
  assert(roePos <= ROE_SIZE);
  roePos < wantROE ? kCarryAROE(out, in, roePos++)
                   : kCarryA(out, in, updateCarryPos(1 << 2));
}

void Gpu::carryAParity(Buffer<Word>& out, Buffer<double>& in) {
  assert(parityCorrect);
  waitParity();
  assert(roePos <= ROE_SIZE);
  Buffer<u32>& expected = parityPrepared || parityPacked ? bufParityExpected : *parityIn;
  kCarryAParity.setFixedArgs(6, expected);
  kCarryAROEParity.setFixedArgs(6, expected);
  roePos < wantROE ? kCarryAROEParity(out, in, roePos++)
                   : kCarryAParity(out, in, updateCarryPos(1 << 2));
}

void Gpu::carryM(Buffer<Word>& out, Buffer<double>& in) {
  assert(roePos <= ROE_SIZE);
  roePos < wantROE ? kCarryMROE(out, in, roePos++)
                   : kCarryM(out, in, updateCarryPos(1 << 3));
}

void Gpu::carryLL(Buffer<Word>& out, Buffer<double>& in) {
  kCarryLL(out, in, updateCarryPos(1 << 2));
}

void Gpu::carryFused(Buffer<double>& buf) {
  // Detached M31 edge: replay the exact M31 inverse width at the end of the
  // M31 bottom half, where it hides under the longer M61 chain and fills
  // bufDetach31 for the detached fused carry below.
  if (detachedM31) { recorded_kernels.push_back(KDETACHA); recorded_kernel_args.push_back(&buf); }
  // This kernel always ends the "bottom half".  Replay the recorded kernel calls.
  endBottomHalf();
  // Only foldP reads bufFolded.  Do not serialize the whole shortened
  // convolution here: its private middle/tail/output work can overlap this
  // carry and is consumed later with finer-grained readiness.
  if (foldTransformEnabled) foldTransform->protectInput(queue);
  // Like fftP, if not in place write the output to the scratch buffer
  Buffer<double> const *in = &buf;
  Buffer<double> const *out = in_place ? &buf : &buf3;
  assert(roePos <= ROE_SIZE);
  if (parityCorrect) {
    waitParity();
    Buffer<u32>& expected = parityPrepared || parityPacked ? bufParityExpected : *parityIn;
    kCarryFused.setFixedArgs(9, expected, *parityOut);
    kCarryFusedROE.setFixedArgs(9, expected, *parityOut);
    roePos < wantROE ? kCarryFusedROE(*out, *in, roePos++)
                     : kCarryFused(*out, *in, updateCarryPos(1 << 0));
    swap(parityIn, parityOut);
    if (foldValidate && !foldValidated) {
      bufFoldMismatches.zero();
      kFoldValidate(bufFoldWords, bufFolded, bufFoldMismatches);
      u32 const mismatches = bufFoldMismatches.read(1)[0];
      if (mismatches != 0) {
        throw std::runtime_error("folded M31 carry output disagrees with standalone fold at " +
                                 std::to_string(mismatches) + " bins");
      }
      log("Folded M31 carry output validated at all %u bins\n", N / 16);
      foldValidated = true;
    }
    if (parityPrepared || parityPacked) prepareParity();
  } else {
    roePos < wantROE ? kCarryFusedROE(*out, *in, roePos++)
                     : kCarryFused(*out, *in, updateCarryPos(1 << 0));
  }
  if (foldTransformEnabled) foldTransform->launch(queue, bufFolded);
  // Detached M31 edge: the fused carry wrote flat pre-width M31 residues to
  // bufDetach31.  Replay the forward width first in the next M31 bottom half,
  // ahead of fftMiddleInGF31 which consumes its output.
  if (detachedM31 == 1) { recorded_kernels.push_back(KDETACHB); recorded_kernel_args.push_back(&buf); }
}

void Gpu::carryFusedMul(Buffer<double>& buf) {
  if (detachedM31) { recorded_kernels.push_back(KDETACHA); recorded_kernel_args.push_back(&buf); }
  // This kernel always ends the "bottom half".  Replay the recorded kernel calls.
  endBottomHalf();
  // Like fftP, if not in place write the output to the scratch buffer
  Buffer<double> const *in = &buf;
  Buffer<double> const *out = in_place ? &buf : &buf3;
  assert(roePos <= ROE_SIZE);
  roePos < wantROE ? kCarryFusedMulROE(*out, *in, roePos++)
                   : kCarryFusedMul(*out, *in, updateCarryPos(1 << 1));
  if (detachedM31 == 1) { recorded_kernels.push_back(KDETACHB); recorded_kernel_args.push_back(&buf); }
}

void Gpu::carryFusedLL(Buffer<double>& buf) {
  if (detachedM31) { recorded_kernels.push_back(KDETACHA); recorded_kernel_args.push_back(&buf); }
  // This kernel always ends the "bottom half".  Replay the recorded kernel calls.
  endBottomHalf();
  // Like fftP, if not in place write the output to the scratch buffer
  Buffer<double> const *in = &buf;
  Buffer<double> const *out = in_place ? &buf : &buf3;
  if (parityCorrect) {
    waitParity();
    Buffer<u32>& expected = parityPrepared || parityPacked ? bufParityExpected : *parityIn;
    kCarryFusedLL.setFixedArgs(9, expected, *parityOut);
    kCarryFusedLL(*out, *in, updateCarryPos(1 << 0));
    swap(parityIn, parityOut);
    if (parityPrepared || parityPacked) prepareParity();
  } else {
    kCarryFusedLL(*out, *in, updateCarryPos(1 << 0));
  }
  if (detachedM31 == 1) { recorded_kernels.push_back(KDETACHB); recorded_kernel_args.push_back(&buf); }
}

void Gpu::carryMiddle(Buffer<Word>& packed, Buffer<double>& buf) {
  assert(useMiddleCarry);
  endBottomHalf();
  assert(roePos <= ROE_SIZE);
  if (roePos < wantROE) {
    kCarryMiddleOut(buf, packed, bufCarry, roePos++);
  } else {
    kCarryMiddleOut(buf, packed, bufCarry, updateCarryPos(1 << 0));
  }
  kCarryMiddleIn(buf, packed, bufCarry);
}

void Gpu::carrySplit(Buffer<Word>& packed, Buffer<double>& buf) {
  assert(useSplitCarry);
  endBottomHalf();
  assert(roePos <= ROE_SIZE);
  if (roePos < wantROE) {
    kCarrySplitOut(buf, packed, bufCarry, roePos++);
  } else {
    kCarrySplitOut(buf, packed, bufCarry, updateCarryPos(1 << 0));
  }
  kCarrySplitIn(buf, packed, bufCarry);
}


#if 0
void Gpu::measureTransferSpeed() {
  u32 SIZE_MB = 16;
  vector<double> data(SIZE_MB * 1024 * 1024, 1);
  Buffer<double> buf{profile.make("DMA"), &queue, SIZE};

  Timer t;
  for (int i = 0; i < 4; ++i) {
    buf.write(data);
    log("buffer Write : %f GB/s\n", double(SIZE / 1024 / 1024) * sizeof(double) / (1024 * t.reset()));
  }

  for (int i = 0; i < 4; ++i) {
    buf.read(data);
    // queue.finish();
    log("buffer READ : %f GB/s\n", double(SIZE / 1024 / 1024) * sizeof(double) / (1024 * t.reset()));
  }

  queue.finish();
}
#endif

u32 Gpu::updateCarryPos(u32 bit) {
  return (statsBits & bit) && (carryPos < CARRY_SIZE) ? carryPos++ : carryPos;
}

vector<Buffer<Word>> Gpu::makeBufVector(u32 size) {
  vector<Buffer<Word>> r;
  r.reserve(size);
for (u32 i = 0; i < size; ++i) { r.emplace_back(timeBufVect, &queue, N); }
  return r;
}

pair<RoeInfo, RoeInfo> Gpu::readROE() {
  assert(roePos <= ROE_SIZE);
  if (roePos) {
    vector<float> roe = bufROE.read(roePos + 2);
    assert(roe.size() == roePos + 2);
    // Experimental FP32+M61 diagnostics store an integer count in each ROE
    // slot.  Convert the raw atomic-add result to a numeric float before the
    // normal statistics machinery consumes it.  Normal ROE collection is
    // bit-for-bit unchanged when ROE_COUNT is not requested.
    if (args.value("ROE_COUNT", 0)) {
      u64 sum = 0;
      u32 countMax = 0;
      u32 valid = 0;
      u32 nonzero = 0;
      u32 otherKernelSlots = 0;
      for (size_t i = 2; i < roe.size(); ++i) {
        u32 count;
        static_assert(sizeof(count) == sizeof(roe[i]));
        memcpy(&count, &roe[i], sizeof(count));
        // Non-fused carry kernels still write an IEEE-754 ROE value to the
        // shared diagnostic buffer at block boundaries.  Its bit pattern is
        // far above the maximum possible coefficient-pair count; exclude
        // those sparse slots from this experiment.
        if (count <= N / 2) {
          sum += count;
          countMax = max(countMax, count);
          valid++;
          nonzero += count != 0;
          roe[i] = float(count);
        } else {
          otherKernelSlots++;
          roe[i] = 0;
        }
      }
      log("ROE_COUNT >= %.3f: %u valid iterations, %u nonzero, max %u, mean %.3f; %u non-fused slots excluded\n",
          args.value("ROE_COUNT", 0) * .001, valid, nonzero, countMax,
          valid ? double(sum) / valid : 0.0, otherKernelSlots);
    }
    // Split the roe buffer into two.  One for squarings and one for multiplications.  This is likely overkill as the multiplication ROE is not used - though
    // it could be useful for debugging (in which case we could support getting roe for squarings or multipplications, but not both).
    auto [squareRoe, mulRoe] = split(roe, mulRoePos);
    // Delete first two used to calculate roePos on the GPU.  Do this after splitting the vector (mulRoePos recorded indices in "+ 2" format).
    u32 squareRoeSize = u32(squareRoe.size()) - 2;
    roe[0] = squareRoe[squareRoeSize];
    roe[1] = squareRoe[squareRoeSize+1];
    squareRoe.resize(squareRoeSize);
    // Clear the ROE buffer and mulRoePos vector
    bufROE.zero(roePos + 2);
    roePos = 0;
    mulRoePos.clear();
    return {roeStat(squareRoe), roeStat(mulRoe)};
  } else {
    return {};
  }
}

RoeInfo Gpu::readCarryStats() {
  assert(carryPos <= CARRY_SIZE);
  if (carryPos == 0) { return {}; }
  vector<float> carry = bufStatsCarry.read(carryPos + 2);
  assert(carry.size() == carryPos + 2);
  // Delete first two used to calculate carryPos on the GPU.
  carry[0] = carry[carryPos];
  carry[1] = carry[carryPos+1];
  carry.resize(carryPos);
  // Clear the GPU buffer
  bufStatsCarry.zero(carryPos + 2);
  carryPos = 0;

  RoeInfo ret = roeStat(carry);

#if 0
  log("%s\n", ret.toString().c_str());

  std::sort(carry.begin(), carry.end());
  File fo = File::openAppend("carry.txt");
  auto it = carry.begin();
  u32 n = carry.size();
  u32 c = 0;
  for (int i=0; i < 500; ++i) {
    double y = 0.23 + (0.48 - 0.23) / 500 * i;
    while (it < carry.end() && *it < y) {
      ++c;
      ++it;
    }
    fo.printf("%f %f\n", y, c / double(n));
  }

  // for (auto x : carry) { fo.printf("%f\n", x); }
  fo.printf("\n\n");
#endif

  return ret;
}

template<typename T>
static bool isAllZero(vector<T> v) { return std::all_of(v.begin(), v.end(), [](T x) { return x == 0;}); }

// Read from GPU, verifying the transfer with a sum, and retry on failure.
vector<Word> Gpu::readChecked(Buffer<Word>& buf) {
  for (int nRetry = 0; nRetry < 3; ++nRetry) {
    bufSumOut.zero();
    sum64(bufSumOut, N, buf);
    vector<u64> expectedVect(1);
    bufSumOut.readAsync(expectedVect);

    vector<Word> data = readOut(buf);
    u64 hostSum = 0;
    for (auto it = data.begin(), end = data.end(); it < end; ++it) hostSum += u64(*it);

    u64 const gpuSum = expectedVect[0];
    if (hostSum == gpuSum) {
      // A buffer containing all-zero is exceptional, so mark that through the empty vector.
      if (gpuSum == 0 && isAllZero(data)) {
        log("Read ZERO\n");
        return {};
      }
      return data;
    }

    log("GPU read failed: %016" PRIx64 " (gpu) != %016" PRIx64 " (host)\n", gpuSum, hostSum);
  }
  throw "GPU persistent read errors";
}

Words Gpu::readAndCompress(Buffer<Word>& buf)  { return compactBits(readChecked(buf), E); }
vector<u32> Gpu::readCheck() { return readAndCompress(bufCheck); }
vector<u32> Gpu::readData() { return readAndCompress(bufData); }

// ioA := ioA * inB; inB must be the output of fftMidIn; inB is preserved
void Gpu::mul(Buffer<Word>& ioA, Buffer<double>& inB, Buffer<double>& tmp1, bool mul3) {
  fftP(tmp1, ioA);
  fftMidIn(tmp1);
  tailMul(tmp1, inB);
  fftMidOut(tmp1);
  fftW(buf3, tmp1);

  // Register the current ROE pos as multiplication (vs. a squaring)
  if (mulRoePos.empty() || mulRoePos.back() < roePos) { mulRoePos.push_back(roePos + 2); }

  if (mul3) { carryM(ioA, buf3); } else { carryA(ioA, buf3); }
  carryB(ioA);
}

// ioA := ioA * inB; inB will end up in buf1 in the LEAD_MIDDLE state
void Gpu::modMul(Buffer<Word>& ioA, Buffer<Word>& inB, bool mul3) {
  modMul(ioA, inB, LEAD_NONE, mul3);
};

// ioA := ioA * inB; inB will end up in buf1 in the LEAD_MIDDLE state
void Gpu::modMul(Buffer<Word>& ioA, Buffer<Word>& inB, enum LEAD_TYPE leadInB, bool mul3) {
  if (leadInB == LEAD_NONE) fftP(buf1, inB);
  if (leadInB != LEAD_MIDDLE) fftMidIn(buf1);
  mul(ioA, buf1, buf2, mul3);
};

void Gpu::writeState(u64 k, const vector<u32>& check, u32 blockSize) {
  assert(blockSize > 0);
  writeIn(bufCheck, check);

  bufData << bufCheck;
  bufAux  << bufCheck;

  if (k) {  // Only verify bufData that was read in from a save file
    u32 n;
    for (n = 1; blockSize % (2 * n) == 0; n *= 2) {
      squareLoop(bufData, 0, n);
      modMul(bufData, bufAux);
      bufAux << bufData;
    }

    assert((n & (n - 1)) == 0);
    assert(blockSize % n == 0);

    blockSize /= n;
    assert(blockSize >= 2);

    for (u32 i = 0; i < blockSize - 2; ++i) {
      squareLoop(bufData, 0, n);
      modMul(bufData, bufAux);
    }

    squareLoop(bufData, 0, n);
  }
  // The middle-3/7 bring-up has not yet generalized tailMul's Hermitian
  // channel pairing.  GT9 has a complete multiply path and uses normal
  // initialization below.
  if (!k && (args.value("GOOD_THOMAS3", 0) ||
             args.value("GOOD_THOMAS7", 0)))
    writeIn(bufData, makeWords(E, 3));
  else modMul(bufData, bufAux, true);
}

bool Gpu::doCheck(u32 blockSize) {
  squareLoop(bufAux, bufCheck, 0, blockSize, true);
  modMul(bufCheck, bufData);
  return isEqual(bufCheck, bufAux);
}

void Gpu::logTimeKernels() {
  auto prof = profile.get();
  u64 total = 0;
  for (const TimeInfo* p : prof) { total += p->times[2]; }
  if (!total) { return; } // no profile
  
  char buf[256];
  // snprintf(buf, sizeof(buf), "Profile:\n ");

  string s = "Profile:\n";
  for (const TimeInfo* p : prof) {
    u32 const n = p->n;
    assert(n);
    double const f = 1e-3 / n;
    double const percent = 100.0 / total * p->times[2];
    if (!args.verbose && percent < 0.2) { break; }
    snprintf(buf, sizeof(buf),
             args.verbose ? "%s %5.2f%% %-18s : %6.1f us/call x %5d calls  (%.3f %.0f)\n"
                          : "%s %5.2f%% %-18s %6.1f x%6d  %.3f %.0f\n",
             logContext().c_str(),
             percent, p->name.c_str(), p->times[2] * f, n, p->times[0] * (f * 1e-3), p->times[1] * (f * 1e-3));
    s += buf;
  }
  log("%s", s.c_str());
  // log("Total time %.3fs\n", total * 1e-9);
  profile.reset();
}

vector<Word> Gpu::readWords(Buffer<Word> &buf) {
  // GPU is returning either 4-byte or 8-byte integers.  C++ code is expecting 8-byte integers.  Handle the "no conversion" case.
  if (fft.WordSize == 8) return buf.read();
  // Convert 32-bit GPU Words into 64-bit C++ Words
  vector<Word> GPUdata = buf.read();
  vector<Word> CPUdata;
  CPUdata.resize(GPUdata.size() * 2);
  for (u32 i = 0; i < GPUdata.size(); ++i) {
    CPUdata[2*i] = (i32) GPUdata[i];
    CPUdata[2*i+1] = (GPUdata[i] >> 32);
  }
  return CPUdata;
}

void Gpu::writeWords(Buffer<Word>& buf, vector<Word> &words) {
  // GPU is expecting either 4-byte or 8-byte integers.  C++ code is using 8-byte integers.  Handle the "no conversion" case.
  if (fft.WordSize == 8) buf.write(words);
  // Convert 64-bit C++ Words into 32-bit GPU Words
  else {
    vector<Word> GPUdata;
    GPUdata.resize(words.size() / 2);
    assert((words.size() & 1) == 0);
    for (u32 i = 0; i < words.size(); i += 2) {
      GPUdata[i/2] = ((i64) words[i+1] << 32) | (u32) words[i];
    }
    buf.write(GPUdata);
  }
}

vector<Word> Gpu::readOut(Buffer<Word> &buf) {
  transpOut(bufAux, buf);
  return readWords(bufAux);
}

void Gpu::writeIn(Buffer<Word>& buf, const vector<u32>& words) { writeIn(buf, expandBits(words, N, E)); }

void Gpu::writeIn(Buffer<Word>& buf, vector<Word>&& words) {
  writeWords(bufAux, words);
  transpIn(buf, bufAux);
}

Words Gpu::expExp2(const Words& A, u32 n) {
  u32 const logStep   = 10000;
  u32 const blockSize = 100;
  
  writeIn(bufData, A);
  IterationTimer timer{0};
  u32 k = 0;
  while (k < n) {
    u32 const its = std::min(blockSize, n - k);
    squareLoop(bufData, 0, its);
    k += its;
    queue.finish();
    if (k % logStep == 0) {
      float const secsPerIt = timer.reset(k);
      log("%u / %u, %s us/it\n", k, n, formatSecsPerIter(secsPerIt).c_str());
    }
  }
  return readData();
}

// A:= A^h * B
void Gpu::expMul(Buffer<Word>& A, u64 h, Buffer<Word>& B) {
  exponentiate(A, h);
  modMul(A, B);
}

// return A^x * B
Words Gpu::expMul(const Words& A, u64 h, const Words& B, bool doSquareB) {
  writeIn(bufCheck, B);
  if (doSquareB) { square(bufCheck); }

  writeIn(bufData, A);
  expMul(bufData, h, bufCheck);
  return readData();
}

static bool testBit(u64 x, int bit) { return x & (u64(1) << bit); }

// See "left-to-right binary exponentiation" on wikipedia
void Gpu::exponentiate(Buffer<Word>& bufInOut, u64 exp) {
  if (exp == 0) {
    bufInOut.set(1);
  } else if (exp > 1) {
    fftP(buf1, bufInOut);
    fftMidIn(buf1);
    fftHin(buf2, buf1); // save fully FFTed "base" to buf2
    bool midInAlreadyDone = true;

    int p = 63;
    while (!testBit(exp, p)) { --p; }

    for (--p; ; --p) {
      if (!midInAlreadyDone) fftMidIn(buf1);
      tailSquare(buf1);
      fftMidOut(buf1);
      midInAlreadyDone = false;

      if (testBit(exp, p)) {
        doCarry(buf1, bufInOut);
        fftMidIn(buf1);
        tailMulLow(buf1, buf2);
        fftMidOut(buf1);
      }

      if (!p) { break; }

      doCarry(buf1, bufInOut);
    }

    fftW(buf3, buf1);
    carryA(bufInOut, buf3);
    carryB(bufInOut);
  }
}

// does either carryFused() or the expanded version depending on useLongCarry
void Gpu::doCarry(Buffer<double>& in, Buffer<Word>& wordBuf) {
  if (useLongCarry) {
    fftW(buf3, in);
    carryA(wordBuf, buf3);
    carryB(wordBuf);
    fftP(in, wordBuf);
  } else {
    carryFused(in);
  }
}

// Use buf1 (and buf23 if not in place) to do a single squaring.
void Gpu::square(Buffer<Word>& out, Buffer<Word>& in, enum LEAD_TYPE leadIn, enum LEAD_TYPE leadOut, bool doMul3, bool doLL) {
  assert(leadOut != LEAD_MIDDLE || useMiddleCarry);
  // LL does not do Mul3
  assert(!(doMul3 && doLL));

  // Use CUDA graphs for some common squarings
  // NOTE: assumes that if doLL is set, it will always be set
  bool graph_recording = false;
  Graph *graph = NULL;
  if (use_graphs && (&out == &bufData || &out == &bufAux) && &in == &out && leadIn == LEAD_WIDTH && leadOut == LEAD_WIDTH && !doMul3) {
    // We have one graph for ROE and one for no-ROE and one for bufData and one for bufAux
    bool roe = (roePos < wantROE);
    bool srcData = (&out == &bufData);
    graph = &graph_square[2 * roe + srcData];
    // Execute an already recorded graph
    if (graph->isRecorded()) {
      graph->launch(&queue);
      queue.incSquareCount();
      if (roe) roePos++;   // WARNING: If we ever graph Gpu::Mul, we'll need to also maintain mulRoePos vector.
      return;
    }
    // Otherwise, record a new graph
    graph->beginRecording(&queue);
    graph_recording = true;
  }

  // In place FFTs use buf1.  Not in place FFTs also use buf3.
  // If leadIn is LEAD_NONE, in contains the input data, squaring starts at fftP
  // If leadIn is LEAD_WIDTH, buf1 (or buf3 if not in place) contains the input data, squaring starts at fftMidIn
  // If leadIn is LEAD_MIDDLE, buf1 contains the input data, squaring starts at tailSquare
  // If leadOut is LEAD_WIDTH, then buf1 (or buf3 if not in place) will contain the output of carryFused -- to be used as input to the next squaring.
  if (leadIn == LEAD_NONE) {
    if (parityCorrect) {
      parityInit(*parityIn, in);
      if (parityPrepared || parityPacked) prepareParity();
    }
    fftP(buf1, in);
  }
  if (leadIn != LEAD_MIDDLE) fftMidIn(buf1);
  tailSquare(buf1);

  // The experimental middle/carry path consumes tail output directly and
  // leaves the next iteration at the tail-input boundary.
  if (leadOut == LEAD_MIDDLE) {
    assert(!doLL && !doMul3);
    carryMiddle(out, buf1);
  } else {
    fftMidOut(buf1);

    // If leadOut is not allowed then we cannot use the faster carryFused kernel
    if (leadOut == LEAD_NONE) {
      fftW(buf3, buf1);
      if (!doLL && !doMul3) {
        parityCorrect ? carryAParity(out, buf3) : carryA(out, buf3);
      } else if (doLL) {
        carryLL(out, buf3);
      } else {
        carryM(out, buf3);
      }
      carryB(out);
    }

    // Use CarryFused
    else {
      assert(!useLongCarry);
      assert(!doMul3);
      if (doLL) {
        carryFusedLL(buf1);
      } else if (useSplitCarry) {
        carrySplit(out, buf1);
      } else {
        carryFused(buf1);
      }
    }
  }

  // End CUDA graph recording (and execute the just recorded graph)
  if (graph_recording) {
    graph->endRecording(&queue);
    graph->launch(&queue);
  }
}

u64 Gpu::squareLoop(Buffer<Word>& out, Buffer<Word>& in, u64 from, u64 to, bool doTailMul3) {
  assert(from < to);
  enum LEAD_TYPE leadIn = LEAD_NONE;
  for (u64 k = from; k < to; ++k) {
    enum LEAD_TYPE const leadOut = useLongCarry || (k == to - 1) ? LEAD_NONE :
                                   useMiddleCarry ? LEAD_MIDDLE : LEAD_WIDTH;
    square(out, (k==from) ? in : out, leadIn, leadOut, doTailMul3 && (k == to - 1));
    leadIn = leadOut;
  }
  return to;
}

bool Gpu::isEqual(Buffer<Word>& in1, Buffer<Word>& in2) {
  kernIsEqual(in1, in2);
  int isEq = 0;
  bufTrue.read(&isEq, 1);
  if (!isEq) { bufTrue.write({1}); }
  return isEq;
}

u64 Gpu::bufResidue(Buffer<Word> &buf) {
  readResidue(bufSmallOut, buf);
  vector<Word> words = readWords(bufSmallOut);

  int carry = 0;
  for (int i = 0; i < 32; ++i) {
    u32 const len = bitlen(N, E, N - 32 + i);
    i64 const w = (i64) words[i] + carry;
    carry = (int) (w >> len);
  }

  u64 res = 0;
  int hasBits = 0;
  for (int k = 0; k < 32 && hasBits < 64; ++k) {
    u32 const len = bitlen(N, E, k);
    i64 const tmp = (i64) words[32 + k] + carry;
    carry = (int) (tmp >> len);
    u64 const w = tmp - ((i64) carry << len);
    assert(w < (1ULL << len));
    res += w << hasBits;
    hasBits += len;
  }
  return res;
}

static string formatETA(u32 secs) {
  u32 const etaMins = (secs + 30) / 60;
  int const days  = etaMins / (24 * 60);
  int const hours = etaMins / 60 % 24;
  int const mins  = etaMins % 60;
  char buf[64];
  if (days) {
    snprintf(buf, sizeof(buf), "%dd %02d:%02d", days, hours, mins);
  } else {
    snprintf(buf, sizeof(buf), "%02d:%02d", hours, mins);
  }
  return string(buf);  
}

static string getETA(u64 step, u64 total, float secsPerStep) {
  u32 const etaSecs = max(0u, u32((total - step) * secsPerStep));
  return formatETA(etaSecs);
}

string RoeInfo::toString() const {
  if (!N) { return {}; }

  char buf[256];
  snprintf(buf, sizeof(buf), "Z(%u)=%.1f Max %f mean %f sd %f (%f, %f)",
           N, z(.5f), max, mean, sd, gumbelMiu, gumbelBeta);
  return buf;
}

static string makeLogStr(const string& status, u64 k, u64 res, float secsPerIt, u64 nIters) {
  char buf[256];
  
  snprintf(buf, sizeof(buf), "%2s %9" PRIu64 " %016" PRIx64 " %s ETA %s; ",
           status.c_str(), k, res, /* k / float(nIters) * 100, */
           formatSecsPerIter(secsPerIt).c_str(), getETA(k, nIters, secsPerIt).c_str());
  return buf;
}

void Gpu::doBigLog(u64 k, u64 res, bool checkOK, float secsPerIt, u64 nIters, u32 nErrors) {
  auto [roeSq, roeMul] = readROE();
  double const z = roeSq.z();
  zAvg.update(z, roeSq.N);
  if (roeSq.max > 0.005)
    log("%sZ=%.0f (avg %.1f), ROEmax=%.3f, ROEavg=%.3f. %s\n", makeLogStr(checkOK ? "OK" : "EE", k, res, secsPerIt, nIters).c_str(),
        z, zAvg.avg(), roeSq.max, roeSq.mean, (nErrors ? " "s + to_string(nErrors) + " errors"s : ""s).c_str());
  else
    log("%sZ=%.0f (avg %.1f) %s\n", makeLogStr(checkOK ? "OK" : "EE", k, res, secsPerIt, nIters).c_str(),
        z, zAvg.avg(), (nErrors ? " "s + to_string(nErrors) + " errors"s : ""s).c_str());

  if (roeSq.N > 2 && (z < 6 || (fft.shape.fft_type == FFT64 && z < 20))) {
    log("Danger ROE! Z=%.1f is too small, increase precision or FFT size!\n", z);
  }

  // Unless ROE log is not explicitly requested, measure only a few iterations to minimize overhead
  wantROE = args.logROE ? ROE_SIZE : 400;

  RoeInfo const carryStats = readCarryStats();
  if (carryStats.N > 2) {
    if (args.value("SPIN_STATS", 0)) {
      log("Carry readiness loops: max %.0f mean %.3f sd %.3f over %u iterations\n",
          carryStats.max, carryStats.mean, carryStats.sd, carryStats.N);
    } else {
      u32 const m = u32(ldexp(carryStats.max, 32));
      double const z = carryStats.z();
      log("Carry: %x Z(%u)=%.1f\n", m, carryStats.N, z);
    }
  }
}

bool Gpu::equals9(const Words& a) {
  if (a[0] != 9) { return false; }
  for (auto it = next(a.begin()); it != a.end(); ++it) { if (*it) { return false; }}
  return true;
}

[[maybe_unused]] static int ulps(double a, double b) {
  if (a == 0 && b == 0) { return 0; }

  u64 const aa = as<u64>(a);
  u64 const bb = as<u64>(b);
  bool const sameSign = (aa >> 63) == (bb >> 63);
  int const delta = int(sameSign ? bb - aa : bb + aa);
  return delta;
}

[[maybe_unused]] static double trigNorm(double c, double s) {
  double const c2 = c * c;
  double const err = fma(c, c, -c2);
  double const norm = c2 + fma(s, s, err);
  return norm;
}

void Gpu::selftestTrig() {

#if FFT_FP64
  const u32 n = hN / 8;
  testTrig(buf1);
  vector<double> trig = buf1.read(n * 2);
  int sup = 0, sdown = 0;
  int cup = 0, cdown = 0;
  int oneUp = 0, oneDown = 0;
  for (u32 k = 0; k < n; ++k) {
    double c = trig[2*k];
    double s = trig[2*k + 1];

#if 0
    auto [refCos, refSin] = root1(hN, k);
#else
    long double angle = M_PIl * k / (hN/2);
    double refSin = sinl(angle);
    double refCos = cosl(angle);
#endif

    if (s > refSin) { ++sup; }
    if (s < refSin) { ++sdown; }
    if (c > refCos) { ++cup; }
    if (c < refCos) { ++cdown; }
    
    double norm = trigNorm(c, s);

    if (norm < 1.0) { ++oneDown; }
    if (norm > 1.0) { ++oneUp; }
  }

  log("TRIG sin(): imperfect %d / %d (%.2f%%), balance %d\n",
      sup + sdown, n, (sup + sdown) * 100.0 / n, sup - sdown);
  log("TRIG cos(): imperfect %d / %d (%.2f%%), balance %d\n",
      cup + cdown, n, (cup + cdown) * 100.0 / n, cup - cdown);
  log("TRIG norm: up %d, down %d\n", oneUp, oneDown);
#endif

  if (isAmdGpu(shared.context->deviceId())) {
    vector<string> WHATS {"V_NOP", "V_ADD_I32", "V_FMA_F32", "V_ADD_F64", "V_FMA_F64", "V_MUL_F64", "V_MAD_U64_U32"};
    for (int w = 0; std::cmp_less(w, WHATS.size()); ++w) {
      const int what = w;
      testTime(what, bufCarry);
      vector<i64> const times = bufCarry.read(4096 * 2);
      [[maybe_unused]] i64 const prev = 0;
      u64 min = -1;
      u64 sum = 0;
      for (i64 const x : times) {
        
#if 0
        if (x != prev) {
          log("%4d : %ld\n", i, x);
          prev = x;
        }
#endif
        if (x > 0 && std::cmp_less(x, min)) { min = x; }
        if (x > 0) { sum += x; }
      }
      log("%-15s : %.2f cycles latency; time min: %d; avg %.0f\n",
          WHATS[w].c_str(), double(min - 40) / 48, int(min), double(sum) / times.size());
    }
  }
}

static u32 mod3(const std::vector<u32> &words) {
  u32 r = 0;
  // uses the fact that 2**32 % 3 == 1.
  for (u32 const w : words) { r += w % 3; }
  return r % 3;
}

static void doDiv3(u64 E, Words& words) {
  u32 r = (3 - mod3(words)) % 3;
  assert(r < 3);
  int const topBits = E % 32;
  assert(topBits > 0 && topBits < 32);
  {
    u64 const w = (u64(r) << topBits) + words.back();
    words.back() = u32(w / 3);
    r = w % 3;
  }
  for (auto it = words.rbegin() + 1, end = words.rend(); it != end; ++it) {
    u64 const w = (u64(r) << 32) + *it;
    *it = u32(w / 3);
    r = w % 3;
  }
}

void Gpu::doDiv9(u64 E, Words& words) {
  doDiv3(E, words);
  doDiv3(E, words);
}

fs::path Gpu::saveProof(const Args& args, ProofSet& proofSet) {
  bool problem_proof = false;
  for ( ; ; ) {
    for (int retry = 0; retry == 0 || (retry == 1 && !problem_proof); ++retry) {
      auto [proof, hashes] = proofSet.computeProof(this);
      fs::path const tmpFile = proof.file(args.proofToVerifyDir);
      proof.save(tmpFile);
            
      fs::path proofFile = proof.file(args.proofResultDir);

      bool const ok = Proof::load(tmpFile).verify(this, hashes);
      log("Proof '%s' verification %s\n", tmpFile.string().c_str(), ok ? "OK" : "FAILED");
      if (ok) {
        fancyRename(tmpFile, proofFile);
        log("Proof '%s' generated\n", proofFile.string().c_str());
        return proofFile;
      }
    }
    problem_proof = true;
    proofSet.reducePower();
    if (proofSet.power < 4) break;
  }
  throw "bad proof generation";
}

PRPState Gpu::loadPRP(Saver<PRPState>& saver) {
  for (int nTries = 0; nTries < 2; ++nTries) {
    if (nTries) {
      saver.dropMostRecent();    // Try an earlier savefile
    }

    PRPState state = saver.load();
    writeState(state.k, state.check, state.blockSize);
    u64 const res = dataResidue();

    if (res == state.res64) {
      log("OK %9" PRIu64 " on-load: blockSize %d, %016" PRIx64 "\n", state.k, state.blockSize, res);
      return state;
      // return {loaded.k, loaded.blockSize, loaded.nErrors};
    }

    log("EE %9" PRIu64 " on-load: %016" PRIx64 " vs. %016" PRIx64 "\n", state.k, res, state.res64);

    // The older GT3/GT7 scaffolds still lack a complete persisted-state
    // boundary.  GT9 is exact and must pass the ordinary on-load check.
    if (args.value("GOOD_THOMAS3", 0) || args.value("GOOD_THOMAS7", 0)) {
      log("Good-Thomas prototype: deferring on-load residue validation until the dedicated carry boundary is installed\n");
      return state;
    }

    if (!state.k) { break; }  // We failed on PRP start
  }

  throw "Error on load";
}

u32 Gpu::getProofPower(u64 k) {
  u32 const power = ProofSet::effectivePower(E, args.getProofPow(E), k);

  if (power != args.getProofPow(E)) {
    log("Proof using power %u (vs %u)\n", power, args.getProofPow(E));
  }

  if (!power) {
    log("Proof generation disabled!\n");
  } else {
    log("Proof of power %u requires about %.1fGB of disk space\n", power, ProofSet::diskUsageGB(E, power));
  }
  return power;
}

tuple<bool, RoeInfo> Gpu::measureCarry() {
  u32 blockSize{}, iters{}, warmup{};

  blockSize = 200;
  iters = 2000;
  warmup = 50;

  assert(iters % blockSize == 0);

  u32 k = 0;
  PRPState const state{.exponent=E, .k=0, .blockSize=blockSize, .res64=3, .check=makeWords(E, 1), .nErrors=0};
  writeState(state.k, state.check, state.blockSize);
  {
    u64 const res = dataResidue();
    if (res != state.res64) {
      log("residue expected %016" PRIx64 " found %016" PRIx64 "\n", state.res64, res);
    }
    assert(res == state.res64);
  }

  enum LEAD_TYPE leadIn = LEAD_NONE;
  modMul(bufCheck, bufData, leadIn);
  leadIn = LEAD_MIDDLE;

  enum LEAD_TYPE const leadOut = useLongCarry ? LEAD_NONE :
                                 useMiddleCarry ? LEAD_MIDDLE : LEAD_WIDTH;
  square(bufData, bufData, leadIn, leadOut);
  leadIn = leadOut;
  ++k;

  while (k < warmup) {
    square(bufData, bufData, leadIn, leadOut);
    leadIn = leadOut;
    ++k;
  }

  readCarryStats(); // ignore the warm-up iterations

  if (Signal::stopRequested()) { throw "stop requested"; }

  while (true) {
    while (k % blockSize < blockSize-1) {
      square(bufData, bufData, leadIn, leadOut);
      leadIn = leadOut;
      ++k;
    }
    square(bufData, bufData, leadIn, LEAD_NONE);
    leadIn = LEAD_NONE;
    ++k;

    if (k >= iters) { break; }

    modMul(bufCheck, bufData, leadIn);
    leadIn = LEAD_MIDDLE;
    if (Signal::stopRequested()) { throw "stop requested"; }
  }

  [[maybe_unused]] u64 const res = dataResidue();
  if (Signal::stopRequested()) { throw "stop requested"; }

  bool const ok = doCheck(blockSize);
  auto stats = readCarryStats();

  // log("%s %016" PRIx64 " %s\n", ok ? "OK" : "EE", res, roe.toString(statsBits).c_str());
  return {ok, stats};
}

tuple<bool, u64, RoeInfo, RoeInfo> Gpu::measureROE(bool  /*quick*/) {
  u32 blockSize{}, iters{}, warmup{};

  {
    blockSize = 200;
    iters = 2000;
    warmup = 50;
  }

  assert(iters % blockSize == 0);

  wantROE = ROE_SIZE; // should be large enough to capture fully this measureROE()

  u32 k = 0;
  PRPState const state{.exponent=E, .k=0, .blockSize=blockSize, .res64=3, .check=makeWords(E, 1), .nErrors=0};
  writeState(state.k, state.check, state.blockSize);
  {
    u64 const res = dataResidue();
    if (res != state.res64) {
      log("residue expected %016" PRIx64 " found %016" PRIx64 "\n", state.res64, res);
    }
    assert(res == state.res64);
  }

  enum LEAD_TYPE leadIn = LEAD_NONE;
  modMul(bufCheck, bufData, leadIn);
  leadIn = LEAD_MIDDLE;

  enum LEAD_TYPE const leadOut = useLongCarry ? LEAD_NONE :
                                 useMiddleCarry ? LEAD_MIDDLE : LEAD_WIDTH;
  square(bufData, bufData, leadIn, leadOut);
  leadIn = leadOut;
  ++k;

  while (k < warmup) {
    square(bufData, bufData, leadIn, leadOut);
    leadIn = leadOut;
    ++k;
  }

  readROE(); // ignore the warm-up iterations

  if (Signal::stopRequested()) { throw "stop requested"; }

  while (true) {
    while (k % blockSize < blockSize-1) {
      square(bufData, bufData, leadIn, leadOut);
      leadIn = leadOut;
      ++k;
    }
    square(bufData, bufData, leadIn, LEAD_NONE);
    leadIn = LEAD_NONE;
    ++k;

    if (k >= iters) { break; }

    modMul(bufCheck, bufData, leadIn);
    leadIn = LEAD_MIDDLE;
    if (Signal::stopRequested()) { throw "stop requested"; }
  }

  [[maybe_unused]] u64 const res = dataResidue();
  if (Signal::stopRequested()) { throw "stop requested"; }

  bool const ok = doCheck(blockSize);
  auto roes = readROE();

  wantROE = 0;
  // log("%s %016" PRIx64 " %s\n", ok ? "OK" : "EE", res, roe.toString(statsBits).c_str());
  return {ok, res, roes.first, roes.second};
}

double Gpu::timePRP(int quick) {        // Quick varies from 1 (slowest, longest) to 10 (quickest, shortest)
  u32 blockSize{}, iters{}, warmup{};

  if (quick == 10)     iters =   400, blockSize = 200;
  else if (quick == 9) iters =   600, blockSize = 300;
  else if (quick == 8) iters =   900, blockSize = 300;
  else if (quick == 7) iters =  1200, blockSize = 400;
  else if (quick == 6) iters =  1800, blockSize = 600;
  else if (quick == 5) iters =  3000, blockSize = 1000;
  else if (quick == 4) iters =  5000, blockSize = 1000;
  else if (quick == 3) iters =  8000, blockSize = 1000;
  else if (quick == 2) iters = 12000, blockSize = 1000;
  else if (quick == 1) iters = 20000, blockSize = 1000;
  warmup = 20;

  assert(iters % blockSize == 0);

  u32 k = 0;
  PRPState const state{.exponent=E, .k=0, .blockSize=blockSize, .res64=3, .check=makeWords(E, 1), .nErrors=0};
  writeState(state.k, state.check, state.blockSize);
  assert(dataResidue() == state.res64);

  enum LEAD_TYPE leadIn = LEAD_NONE;
  modMul(bufCheck, bufData, leadIn);
  leadIn = LEAD_MIDDLE;

  enum LEAD_TYPE const leadOut = useLongCarry ? LEAD_NONE :
                                 useMiddleCarry ? LEAD_MIDDLE : LEAD_WIDTH;
  square(bufData, bufData, leadIn, leadOut);
  leadIn = leadOut;
  ++k;

  while (k < warmup) {
    square(bufData, bufData, leadIn, leadOut);
    leadIn = leadOut;
    ++k;
  }
  queue.finish();
  if (Signal::stopRequested()) { throw "stop requested"; }

  Timer t;
  queue.setSquareTime(0);     // Busy wait on nVidia to get the most accurate timings while tuning
  while (true) {
    while (k % blockSize < blockSize-1) {
      square(bufData, bufData, leadIn, leadOut);
      leadIn = leadOut;
      ++k;
    }
    square(bufData, bufData, leadIn, LEAD_NONE);
    leadIn = LEAD_NONE;
    ++k;

    if (k >= iters) { break; }

    modMul(bufCheck, bufData, leadIn);
    leadIn = LEAD_MIDDLE;
    if (Signal::stopRequested()) { throw "stop requested"; }
  }
  queue.finish();
  double secsPerIt = t.reset() / (iters - warmup);

  if (Signal::stopRequested()) { throw "stop requested"; }

  u64 const res = dataResidue();
  bool const ok = doCheck(blockSize);
  if (!ok) {
    log("Error %016" PRIx64 "\n", res);
    secsPerIt = 0.1; // a large value to mark the error
  }
  return secsPerIt * 1e6;
}

PRPResult Gpu::isPrimePRP([[maybe_unused]] const Task& task) {
  assert(E == task.exponent);

  // This timer is used to measure total elapsed time to be written to the savefile.
  Timer elapsedTimer;

  u32 nErrors = 0;
  int nSeqErrors = 0;
  u64 lastFailedRes64 = 0;
  u32 const logStep = args.logStep;

 reload:
  elapsedTimer.reset();
  u32 blockSize{};
  u64 k{};
  double elapsedBefore = 0;

  {
    PRPState const state = loadPRP(*getSaver());
    nErrors = std::max(nErrors, state.nErrors);
    blockSize = state.blockSize;
    k = state.k;
    elapsedBefore = state.elapsed;
  }

  assert(blockSize > 0 && logStep % blockSize == 0);

  u32 const checkStep = checkStepForErrors(blockSize, nErrors);
  assert(checkStep % logStep == 0);

  u32 const power = getProofPower(k);
  
  ProofSet proofSet{E, power};

  bool isPrime = false;

  u64 finalRes64 = 0;
  vector<u32> res2048;

  // We extract the res64 at kEnd.
  // For M=2^E-1, residue "type-3" == 3^(M+1), and residue "type-1" == type-3 / 9,
  // See http://www.mersenneforum.org/showpost.php?p=468378&postcount=209
  // For both type-1 and type-3 we need to do E squarings (as M+1==2^E).
  const u64 kEnd = E;
  assert(k < kEnd);

  // We continue beyound kEnd: to the next multiple of blockSize, to do a check there
  u64 const kEndEnd = roundUp(kEnd, blockSize);

  bool skipNextCheckUpdate = false;

  u64 persistK = proofSet.next(k);
  enum LEAD_TYPE leadIn = LEAD_NONE;

  assert(k % blockSize == 0);
  assert(checkStep % blockSize == 0);

  const u64 startK = k;
  IterationTimer iterationTimer{k};

  wantROE = 0; // skip the initial iterations

  while (true) {
    assert(k < kEndEnd);

    if (!wantROE && k - startK > 30) { wantROE = args.logROE ? ROE_SIZE : 2'000; }

    if (skipNextCheckUpdate) {
      skipNextCheckUpdate = false;
    } else if (k % blockSize == 0) {
      modMul(bufCheck, bufData, leadIn);
      leadIn = LEAD_MIDDLE;
    }

    ++k; // !! early inc

    bool const doStop = (k % blockSize == 0) && (Signal::stopRequested() || (args.iters && k - startK >= args.iters));
    bool const doCheck = doStop || (k % checkStep == 0) || (k >= kEndEnd) || (k - startK == 2 * blockSize);
    bool const doLog = k % logStep == 0;
    enum LEAD_TYPE const leadOut =
      doCheck || doLog || k == persistK || k == kEnd || useLongCarry ? LEAD_NONE :
      useMiddleCarry ? LEAD_MIDDLE : LEAD_WIDTH;

    if (doStop) { log("Stopping, please wait..\n"); }

    square(bufData, bufData, leadIn, leadOut, false);
    leadIn = leadOut;

    if (k == persistK) {
      vector<Word> const rawData = readChecked(bufData);
      if (rawData.empty()) {
        log("Data error ZERO\n");
        ++nErrors;
        goto reload;
      }
      (*background)([=, E=this->E] { ProofSet::save(E, power, k, compactBits(rawData, E)); });
      persistK = proofSet.next(k);
    }

    if (k == kEnd) {
      Words words = readData();
      isPrime = equals9(words);
      doDiv9(E, words);
      finalRes64 = residue(words);
      res2048.clear();
      assert(words.size() >= 64);
      res2048.insert(res2048.end(), words.begin(), std::next(words.begin(), 64));
      log("%s %8" PRIu64 " / %" PRIu64 ", %s\n", isPrime ? "PP" : "CC", kEnd, E, hex(finalRes64).c_str());
    }

    if (!doCheck && !doLog) continue;

    u64 const res = dataResidue();
    float const secsPerIt = iterationTimer.reset(k);
    queue.setSquareTime((int) (secsPerIt * 1'000'000));

    vector<Word> rawCheck = readChecked(bufCheck);
    if (rawCheck.empty()) {
      ++nErrors;
      log("%9" PRIu64 " %016" PRIx64 " read NULL check\n", k, res);
      if (++nSeqErrors > 2) { throw "sequential errors"; }
      goto reload;
    }

    if (!doCheck) {
      (*background)([=, this] {
        getSaver()->saveUnverified({.exponent=E, .k=k, .blockSize=blockSize, .res64=res, .check=compactBits(rawCheck, E), .nErrors=nErrors,
                                    .elapsed=elapsedBefore + elapsedTimer.at()});
      });

      log("   %9" PRIu64 " %016" PRIx64 " %s\n", k, res, formatSecsPerIter(secsPerIt).c_str());
      RoeInfo const carryStats = readCarryStats();
      if (carryStats.N) {
        if (args.value("SPIN_STATS", 0)) {
          log("Carry readiness loops: max %.0f mean %.3f sd %.3f over %u iterations\n",
              carryStats.max, carryStats.mean, carryStats.sd, carryStats.N);
        } else {
          u32 const m = u32(ldexp(carryStats.max, 32));
          double const z = carryStats.z();
          log("Carry: %x Z(%u)=%.1f\n", m, carryStats.N, z);
        }
      }
    } else {
      bool const ok = this->doCheck(blockSize);
      [[maybe_unused]] float const secsCheck = iterationTimer.reset(k);

      if (ok) {
        nSeqErrors = 0;
        // lastFailedRes64 = 0;
        skipNextCheckUpdate = true;

        if (k < kEnd) {
          (*background)([=, this, rawCheck = std::move(rawCheck)] {
            getSaver()->save({.exponent=E, .k=k, .blockSize=blockSize, .res64=res, .check=compactBits(rawCheck, E), .nErrors=nErrors, .elapsed=elapsedBefore + elapsedTimer.at()});
          });
        }

        doBigLog(k, res, ok, secsPerIt, kEndEnd, nErrors);

        if (k >= kEndEnd) {
          fs::path const proofFile = saveProof(args, proofSet);
          return {.isPrime=isPrime, .res64=finalRes64, .nErrors=nErrors, .proofPath=proofFile.string(), .res2048=toHex(res2048)};
        }        
      } else {
        ++nErrors;
        doBigLog(k, res, ok, secsPerIt, kEndEnd, nErrors);
        if (args.value("GOOD_THOMAS3", 0) || args.value("GOOD_THOMAS7", 0) ||
            args.value("GOOD_THOMAS9", 0))
          logTimeKernels();
        if (++nSeqErrors > 2) {
          log("%d sequential errors, will stop.\n", nSeqErrors);
          throw "too many errors";
        }
        if (res == lastFailedRes64) {
          log("Consistent error %016" PRIx64 ", will stop.\n", res);
          throw "consistent error";
        }
        lastFailedRes64 = res;
        if (!doStop) { goto reload; }
      }

      logTimeKernels();

      if (doStop) {
        queue.finish();
        throw "stop requested";
      }

      iterationTimer.reset(k);
    }
  }
}

LLResult Gpu::isPrimeLL([[maybe_unused]] const Task& task) {
  assert(E == task.exponent);
  wantROE = 0;

  Timer elapsedTimer;

  Saver<LLState> saver{E, 1000, args.nSavefiles};

  reload:
  elapsedTimer.reset();

  u64 startK = 0;
  double elapsedBefore = 0;
  {
    LLState state = saver.load();

    elapsedBefore = state.elapsed;
    startK = state.k;
    u64 const expectedRes = (u64(state.data[1]) << 32) | state.data[0];
    writeIn(bufData, state.data);
    u64 const res = dataResidue();
    if (res != expectedRes) { throw "Invalid savefile (res64)"; }
    assert(res == expectedRes);
    log("LL loaded @ %" PRIu64 " : %016" PRIx64 "\n", startK, res);
  }

  IterationTimer iterationTimer{startK};

  u64 k = startK;
  u64 const kEnd = E - 2;
  enum LEAD_TYPE leadIn = LEAD_NONE;

  while (true) {
    ++k;
    bool doStop = (k >= kEnd) || (args.iters && k - startK >= args.iters);

    if (Signal::stopRequested()) {
      doStop = true;
      log("Stopping, please wait..\n");
    }

    bool const doLog = (k % args.logStep == 0) || doStop;
    enum LEAD_TYPE const leadOut = doLog || useLongCarry ? LEAD_NONE :
                                   useMiddleCarry ? LEAD_MIDDLE : LEAD_WIDTH;

    squareLL(bufData, leadIn, leadOut);
    leadIn = leadOut;

    if (!doLog) continue;

    u64 res64 = 0;
    auto data = readData();
    bool const isAllZero = data.empty();

    if (isAllZero) {
      if (k < kEnd) {
        log("Error: early ZERO @ %" PRIu64 "\n", k);
        if (doStop) {
          throw "stop requested";
        } 
          goto reload;
       
      }
      res64 = 0;
    } else {
      assert(data.size() >= 2);
      res64 = (u64(data[1]) << 32) | data[0];
      saver.save({.exponent=E, .k=k, .data=std::move(data), .elapsed=elapsedBefore + elapsedTimer.at()});
    }

    float const secsPerIt = iterationTimer.reset(k);
    queue.setSquareTime((int) (secsPerIt * 1'000'000));
    log("%9" PRIu64 " %016" PRIx64 " %s ETA %s\n", k, res64, formatSecsPerIter(secsPerIt).c_str(), getETA(k, kEnd, secsPerIt).c_str());

    if (k >= kEnd) { return {.isPrime=isAllZero, .res64=res64}; }

    if (doStop) { throw "stop requested"; }
  }
}

array<u64, 4> Gpu::isCERT(const Task& task) {
  assert(E == task.exponent);
  wantROE = 0;

  // Get CERT start value
  char fname[32];
  sprintf(fname, "M%" PRIu64 ".cert", E);

// AutoPrimenet.py does not add the cert entry to worktodo.txt until it has successfully downloaded the .cert file.

  { // Enclosing this code in braces ensures the file will be closed by the File destructor.  The later file deletion requires the file be closed in Windows.
    File fi = File::openReadThrow(fname);
    u32 const nBytes = u32((E - 1) / 8 + 1);
    Words const B = fi.readBytesLE(nBytes);
    writeIn(bufData, B);
  }

  Timer elapsedTimer;

  elapsedTimer.reset();

  u32 const startK = 0;

  IterationTimer iterationTimer{startK};

  u32 k = 0;
  u32 const kEnd = task.squarings;
  enum LEAD_TYPE leadIn = LEAD_NONE;

  while (true) {
    ++k;
    bool doStop = (k >= kEnd);

    if (Signal::stopRequested()) {
      doStop = true;
      log("Stopping, please wait..\n");
    }

    bool const doLog = (k % 100'000 == 0) || doStop;
    enum LEAD_TYPE const leadOut = doLog || useLongCarry ? LEAD_NONE :
                                   useMiddleCarry ? LEAD_MIDDLE : LEAD_WIDTH;

    squareCERT(bufData, leadIn, leadOut);
    leadIn = leadOut;

    if (!doLog) continue;

    Words data = readData();
    assert(data.size() >= 2);
    u64 const res64 = (u64(data[1]) << 32) | data[0];

    float const secsPerIt = iterationTimer.reset(k);
    queue.setSquareTime((int) (secsPerIt * 1'000'000));
    log("%7u / %7u %016" PRIx64 " %s ETA %s\n", k, kEnd, res64, formatSecsPerIter(secsPerIt).c_str(), getETA(k, kEnd, secsPerIt).c_str());

    if (k >= kEnd) {
      fs::remove (fname);
      return std::move(SHA3{}.update(data.data(), u32((E-1)/8+1))).finish();
    }

    if (doStop) { throw "stop requested"; }
  }
}


void Gpu::clear(bool isPRP) {
  if (isPRP) {
    Saver<PRPState>::clear(E);
  } else {
    Saver<LLState>::clear(E);
  }
}

Saver<PRPState> *Gpu::getSaver() {
  if (!saver) { saver = make_unique<Saver<PRPState>>(E, args.blockSize, args.nSavefiles); }
  return saver.get();
}
