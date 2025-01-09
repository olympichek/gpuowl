// Copyright Mihai Preda

#include "TrigBufCache.h"

// Trial balloon.  I'm hoping we can create one #define that indicates a GPU's preference for doing extra DP work in exchange for
// less memory accesses.  I can envision a Titan V (slow memory, 1:2 SP:DP ratio) setting this value to a high value, whereas a
// consumer grade nVidia GPU with 1:32 or 1:64 SP:DP ratio would set this to zero.
#define PREFER_DP_TO_MEM        2               // Excellent DP GPU such as Titan V or Radeon VII Pro.
//#define PREFER_DP_TO_MEM      1               // Good DP GPU.  Tuned for Radeon VII.
//#define PREFER_DP_TO_MEM      0               // Poor DP GPU.  A typical consumer grade GPU.

// Klunky defines for single-wide vs. double-wide tailSquare
// Clean this up once we determine which options to make user visible
#define SINGLE_WIDE             0       // Old single-wide tailSquare vs. new double-wide tailSquare
#define SINGLE_KERNEL           0       // Implement tailSquare in a single kernel vs. two kernels

#define SAVE_ONE_MORE_MUL 1             // I'll probably make this the only option -- still studying if rocm optimizer is making a mess of this in carryfused

#define _USE_MATH_DEFINES
#include <cmath>

#ifndef M_PIl
#define M_PIl 3.141592653589793238462643383279502884L
#endif

#ifndef M_PI
#define M_PI 3.1415926535897931
#endif

static_assert(sizeof(double2) == 16, "size double2");

// For small angles, return "fancy" cos - 1 for increased precision
double2 root1Fancy(u32 N, u32 k) {
  assert(!(N&7));
  assert(k < N);
  assert(k < N/4);

  long double angle = M_PIl * k / (N / 2);
  return {double(cosl(angle) - 1), double(sinl(angle))};
}

static double trigNorm(double c, double s) { return c * c + s * s; }
static double trigError(double c, double s) { return abs(trigNorm(c, s) - 1.0); }

// Round trig long double to double as to satisfy c^2 + s^2 == 1 as best as possible
static double2 roundTrig(long double lc, long double ls) {
  double c1 = lc;
  double c2 = nexttoward(c1, lc);
  double s1 = ls;
  double s2 = nexttoward(s1, ls);

  double c = c1;
  double s = s1;
  for (double tryC : {c1, c2}) {
    for (double tryS : {s1, s2}) {
      if (trigError(tryC, tryS) < trigError(c, s)) {
        c = tryC;
        s = tryS;
      }
    }
  }
  return {c, s};
}

// Returns the primitive root of unity of order N, to the power k.
double2 root1(u32 N, u32 k) {
  assert(k < N);
  if (k >= N/2) {
    auto [c, s] = root1(N, k - N/2);
    return {-c, -s};
  } else if (k > N/4) {
    auto [c, s] = root1(N, N/2 - k);
    return {-c, s};
  } else if (k > N/8) {
    auto [c, s] = root1(N, N/4 - k);
    return {s, c};
  } else {
    assert(k <= N/8);

    long double angle = M_PIl * k / (N / 2);

#if 1
    return roundTrig(cosl(angle), sinl(angle));
#else
    double c = cos(double(angle)), s = sin(double(angle));
    if ((c * c + s * s == 1.0)) {
      return {c, s};
    } else {
      return {double(cosl(angle)), double(sinl(angle))};
    }
#endif
  }
}

namespace {
static const constexpr bool LOG_TRIG_ALLOC = false;

vector<double2> genSmallTrig(u32 size, u32 radix) {
  if (LOG_TRIG_ALLOC) { log("genSmallTrig(%u, %u)\n", size, radix); }

  vector<double2> tab;
// old fft_WIDTH
  for (u32 line = 1; line < radix; ++line) {
    for (u32 col = 0; col < size / radix; ++col) {
      tab.push_back(radix / line >= 8 ? root1Fancy(size, col * line) : root1(size, col * line));
    }
  }
  tab.resize(size);

// New fft_WIDTH
  // Epsilon value, 2^-250, should have an exact representation as a double
  const double epsilon = 5.5271478752604445602472651921923E-76;  // Protect against divide by zero
  // Sine/cosine values for first fft8
//TO DO: explore using long doubles through the division (though dividing by double(cosine) may make sense)
  for (u32 line = 1; line < radix-1; line += 2) {  // Output pairs of lines
    for (u32 col = 0; col < size / radix; ++col) {
      double2 rootA = root1(size, col * line); rootA.first += epsilon;
      double2 rootB = root1(size, col * (line + 1)); rootB.first += epsilon;
      double2 sine_over_cosines = {rootA.second / rootA.first, rootB.second / rootB.first};
      tab.push_back(sine_over_cosines);
    }
  }
  for (u32 line = radix-1; line < radix; ++line) {
    for (u32 col = 0; col < size / radix; col += 2) {
      double2 rootA = root1(size, col * line); rootA.first += epsilon;
      double2 rootB = root1(size, (col + 1) * line); rootB.first += epsilon;
      double2 sine_over_cosines = {rootA.second / rootA.first, rootB.second / rootB.first};
      tab.push_back(sine_over_cosines);
    }
  }
  // Sine/cosine values for second fft8
  for (u32 line = 0; line < radix; ++line) {
    for (u32 col = 0; col < size / radix / 8; col += 2) {
      double2 rootA = root1(size, 8 * col * line); rootA.first += epsilon;
      double2 rootB = root1(size, 8 * (col + 1) * line); rootB.first += epsilon;
      double2 sine_over_cosines = {rootA.second / rootA.first, rootB.second / rootB.first};
      tab.push_back(sine_over_cosines);
    }
  }
  // Cosine values for first fft8 (output in post-shufl order)
//TODO: Examine why when sine is 0.0 cosine is not 1.0 or -1.0 (printf is outputting 0.999... and -0.999...)
  for (u32 col = 0; col < size / radix; ++col) {  // Output pairs of lines, col 0..7 is line0, col 8..15 is line1, etc.
    if (col & 8) continue;                        // Skip line1, line3, etc.
    for (u32 line = 0; line < radix; ++line) {
      double2 rootA = root1(size, col * line); rootA.first += epsilon;
      double2 rootB = root1(size, (col + 8) * line); rootB.first += epsilon;
#if SAVE_ONE_MORE_MUL
      if (col / 8 == 2) { // Compute cosine3 / cosine1
      rootB.first /= root1(size, (col + 8 - 16) * line).first + epsilon;
      }
      if (col / 8 == 4) { // Compute cosine5 / cosine1
      rootB.first /= root1(size, (col + 8 - 32) * line).first + epsilon;
      }
#endif
      if (col / 8 == 6) { // Compute cosine6 / cosine2
        rootA.first /= root1(size, (col - 32) * line).first + epsilon;
      }
      if (col / 8 == 6) { // Compute cosine7 / cosine3
        rootB.first /= root1(size, (col + 8 - 32) * line).first + epsilon;
      }
      double2 cosines = {rootA.first, rootB.first};
      tab.push_back(cosines);
    }
  }
  // Cosine values for second fft8 (output in post-shufl order)
  for (u32 col = 0; col < size / radix / 8; ++col) {
    for (u32 line = 0; line < radix; line += 2) {
      double2 rootA = root1(size, 8 * col * line); rootA.first += epsilon;
      double2 rootB = root1(size, 8 * col * (line + 1)); rootB.first += epsilon;
#if SAVE_ONE_MORE_MUL
      if (col == 3) { // Compute cosine3 / cosine1
        rootA.first /= root1(size, 8 * (col - 2) * line).first + epsilon;
        rootB.first /= root1(size, 8 * (col - 2) * (line + 1)).first + epsilon;
      }
      if (col == 5) { // Compute cosine5 / cosine1
        rootA.first /= root1(size, 8 * (col - 4) * line).first + epsilon;
        rootB.first /= root1(size, 8 * (col - 4) * (line + 1)).first + epsilon;
      }
#endif
      if (col == 6) { // Compute cosine6 / cosine2
        rootA.first /= root1(size, 8 * (col - 4) * line).first + epsilon;
        rootB.first /= root1(size, 8 * (col - 4) * (line + 1)).first + epsilon;
      }
      if (col == 7) { // Compute cosine7 / cosine3
        rootA.first /= root1(size, 8 * (col - 4) * line).first + epsilon;
        rootB.first /= root1(size, 8 * (col - 4) * (line + 1)).first + epsilon;
      }
      double2 cosines = {rootA.first, rootB.first};
      tab.push_back(cosines);
    }
  }

  return tab;
}

// Generate the small trig values for fft_HEIGHT plus optionally trig values used in pairSq.
vector<double2> genSmallTrigCombo(u32 width, u32 middle, u32 size, u32 radix) {
  if (LOG_TRIG_ALLOC) { log("genSmallTrigCombo(%u, %u)\n", size, radix); }

  vector<double2> tab;
// old fft_HEIGHT
  for (u32 line = 1; line < radix; ++line) {
    for (u32 col = 0; col < size / radix; ++col) {
      tab.push_back(radix / line >= 8 ? root1Fancy(size, col * line) : root1(size, col * line));
    }
  }
  tab.resize(size);

// New fft_HEIGHT
  // Epsilon value, 2^-250, should have an exact representation as a double
  const double epsilon = 5.5271478752604445602472651921923E-76;  // Protect against divide by zero
  // Sine/cosine values for first fft8
//TO DO: explore using long doubles through the division (though dividing by double(cosine) may make sense)
  for (u32 line = 1; line < radix-1; line += 2) {  // Output pairs of lines
    for (u32 col = 0; col < size / radix; ++col) {
      double2 rootA = root1(size, col * line); rootA.first += epsilon;
      double2 rootB = root1(size, col * (line + 1)); rootB.first += epsilon;
      double2 sine_over_cosines = {rootA.second / rootA.first, rootB.second / rootB.first};
      tab.push_back(sine_over_cosines);
    }
  }
  for (u32 line = radix-1; line < radix; ++line) {
    for (u32 col = 0; col < size / radix; col += 2) {
      double2 rootA = root1(size, col * line); rootA.first += epsilon;
      double2 rootB = root1(size, (col + 1) * line); rootB.first += epsilon;
      double2 sine_over_cosines = {rootA.second / rootA.first, rootB.second / rootB.first};
      tab.push_back(sine_over_cosines);
    }
  }
  // Sine/cosine values for second fft8
  for (u32 line = 0; line < radix; ++line) {
    for (u32 col = 0; col < size / radix / 8; col += 2) {
      double2 rootA = root1(size, 8 * col * line); rootA.first += epsilon;
      double2 rootB = root1(size, 8 * (col + 1) * line); rootB.first += epsilon;
      double2 sine_over_cosines = {rootA.second / rootA.first, rootB.second / rootB.first};
      tab.push_back(sine_over_cosines);
    }
  }
  // Cosine values for first fft8 (output in post-shufl order)
//TODO: Examine why when sine is 0.0 cosine is not 1.0 or -1.0 (printf is outputting 0.999... and -0.999...)
  for (u32 col = 0; col < size / radix; ++col) {  // Output pairs of lines, col 0..7 is line0, col 8..15 is line1, etc.
    if (col & 8) continue;                        // Skip line1, line3, etc.
    for (u32 line = 0; line < radix; ++line) {
      double2 rootA = root1(size, col * line); rootA.first += epsilon;
      double2 rootB = root1(size, (col + 8) * line); rootB.first += epsilon;
#if SAVE_ONE_MORE_MUL
      if (col / 8 == 2) { // Compute cosine3 / cosine1
      rootB.first /= root1(size, (col + 8 - 16) * line).first + epsilon;
      }
      if (col / 8 == 4) { // Compute cosine5 / cosine1
      rootB.first /= root1(size, (col + 8 - 32) * line).first + epsilon;
      }
#endif
      if (col / 8 == 6) { // Compute cosine6 / cosine2
        rootA.first /= root1(size, (col - 32) * line).first + epsilon;
      }
      if (col / 8 == 6) { // Compute cosine7 / cosine3
        rootB.first /= root1(size, (col + 8 - 32) * line).first + epsilon;
      }
      double2 cosines = {rootA.first, rootB.first};
      tab.push_back(cosines);
    }
  }
  // Cosine values for second fft8 (output in post-shufl order)
  for (u32 col = 0; col < size / radix / 8; ++col) {
    for (u32 line = 0; line < radix; line += 2) {
      double2 rootA = root1(size, 8 * col * line); rootA.first += epsilon;
      double2 rootB = root1(size, 8 * col * (line + 1)); rootB.first += epsilon;
#if SAVE_ONE_MORE_MUL
      if (col == 3) { // Compute cosine3 / cosine1
        rootA.first /= root1(size, 8 * (col - 2) * line).first + epsilon;
        rootB.first /= root1(size, 8 * (col - 2) * (line + 1)).first + epsilon;
      }
      if (col == 5) { // Compute cosine5 / cosine1
        rootA.first /= root1(size, 8 * (col - 4) * line).first + epsilon;
        rootB.first /= root1(size, 8 * (col - 4) * (line + 1)).first + epsilon;
      }
#endif
      if (col == 6) { // Compute cosine6 / cosine2
        rootA.first /= root1(size, 8 * (col - 4) * line).first + epsilon;
        rootB.first /= root1(size, 8 * (col - 4) * (line + 1)).first + epsilon;
      }
      if (col == 7) { // Compute cosine7 / cosine3
        rootA.first /= root1(size, 8 * (col - 4) * line).first + epsilon;
        rootB.first /= root1(size, 8 * (col - 4) * (line + 1)).first + epsilon;
      }
      double2 cosines = {rootA.first, rootB.first};
      tab.push_back(cosines);
    }
  }
  tab.resize(size*4);

  // From tailSquare pre-calculate some or all of these:  T2 trig = slowTrig_N(line + H * lowMe, ND / NH * 2);
#if PREFER_DP_TO_MEM == 2             // No pre-computed trig values
#elif PREFER_DP_TO_MEM == 1           // best option on a Radeon VII.
  u32 height = size;
  // Output line 0 trig values to be read by every u,v pair of lines
  for (u32 me = 0; me < height / radix; ++me) {
    tab.push_back(root1(width * middle * height, width * middle * me));
  }
  // Output the one or two T2 multipliers to be read by one u,v pair of lines
  for (u32 line = 0; line < width * middle / 2; ++line) {
    tab.push_back(root1Fancy(width * middle * height, line));
    if (!SINGLE_WIDE) tab.push_back(root1Fancy(width * middle * height, width * middle - line));
  }
#else
  u32 height = size;
  for (u32 u = 0; u < width * middle / 2; ++u) {
    for (u32 v = 0; v < (SINGLE_WIDE ? 1 : 2); ++v) {
      u32 line = (v == 0) ? u : width * middle - u;
      for (u32 me = 0; me < height / radix; ++me) {
        tab.push_back(root1(width * middle * height, line + width * middle * me));
      }
    }
  }
#endif

  return tab;
}

// starting from a MIDDLE of 5 we consider angles in [0, 2Pi/MIDDLE] as worth storing with the
// cos-1 "fancy" trick.
#define SHARP_MIDDLE 5

vector<double2> genMiddleTrig(u32 smallH, u32 middle, u32 width) {
  if (LOG_TRIG_ALLOC) { log("genMiddleTrig(%u, %u, %u)\n", smallH, middle, width); }
  vector<double2> tab;
  if (middle == 1) {
    tab.resize(1);
  } else {
    if (middle < SHARP_MIDDLE) {
      for (u32 k = 0; k < smallH; ++k) { tab.push_back(root1(smallH * middle, k)); }
      for (u32 k = 0; k < width; ++k)  { tab.push_back(root1(middle * width, k)); }
    } else {
      for (u32 k = 0; k < smallH; ++k) { tab.push_back(root1Fancy(smallH * middle, k)); }
      for (u32 k = 0; k < width; ++k)  { tab.push_back(root1Fancy(middle * width, k)); }
    }
  }
  return tab;
}

} // namespace

TrigBufCache::~TrigBufCache() = default;

TrigPtr TrigBufCache::smallTrig(u32 W, u32 nW) {
  lock_guard lock{mut};
  auto& m = small;
  decay_t<decltype(m)>::key_type key{W, nW};

  TrigPtr p{};
  auto it = m.find(key);
  if (it == m.end() || !(p = it->second.lock())) {
    p = make_shared<TrigBuf>(context, genSmallTrig(W, nW));
    m[key] = p;
    smallCache.add(p);
  }
  return p;
}

TrigPtr TrigBufCache::smallTrigCombo(u32 width, u32 middle, u32 W, u32 nW) {
  lock_guard lock{mut};
  auto& m = small;
#if PREFER_DP_TO_MEM == 2             // No pre-computed trig values
  decay_t<decltype(m)>::key_type key{W, nW};
#else
  // Hack so that width 512 and height 512 don't share the same buffer.  Width could share the height buffer since it is a subset of the combo height buffer.
  // Also add in middle so that 512:15:512 does not share buffer with 512:14:512.
  decay_t<decltype(m)>::key_type key{W, nW+middle};
#endif

  TrigPtr p{};
  auto it = m.find(key);
  if (it == m.end() || !(p = it->second.lock())) {
    p = make_shared<TrigBuf>(context, genSmallTrigCombo(width, middle, W, nW));
    m[key] = p;
    smallCache.add(p);
  }
  return p;
}

TrigPtr TrigBufCache::middleTrig(u32 SMALL_H, u32 MIDDLE, u32 width) {
  lock_guard lock{mut};
  auto& m = middle;
  decay_t<decltype(m)>::key_type key{SMALL_H, MIDDLE, width};

  TrigPtr p{};
  auto it = m.find(key);
  if (it == m.end() || !(p = it->second.lock())) {
    p = make_shared<TrigBuf>(context, genMiddleTrig(SMALL_H, MIDDLE, width));
    m[key] = p;
    middleCache.add(p);
  }
  return p;
}
