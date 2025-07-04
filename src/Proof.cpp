// Copyright (C) Mihai Preda.

#include "Proof.h"
#include "Sha3Hash.h"
#include "MD5.h"
#include "Gpu.h"

#include <algorithm>
#include <vector>
#include <string>
#include <cassert>
#include <filesystem>
#include <cinttypes>
#include <charconv>
#include <gmpxx.h>

#if defined(__BYTE_ORDER__) && __BYTE_ORDER__ != __ORDER_LITTLE_ENDIAN__
#error Byte order must be Little Endian
#endif


namespace proof {

array<u64, 4> hashWords(u32 E, const Words& words) {
  return std::move(SHA3{}.update(words.data(), (E-1)/8+1)).finish();
}

array<u64, 4> hashWords(u32 E, array<u64, 4> prefix, const Words& words) {
  return std::move(SHA3{}.update(prefix).update(words.data(), (E-1)/8+1)).finish();
}

string fileHash(const fs::path& filePath) {
  File fi = File::openReadThrow(filePath);
  char buf[64 * 1024];
  MD5 h;
  u32 size = 0;
  while ((size = fi.readUpTo(buf, sizeof(buf)))) { h.update(buf, size); }
  return std::move(h).finish();
}

ProofInfo getInfo(const fs::path& proofFile) {
  string hash = proof::fileHash(proofFile);
  File fi = File::openReadThrow(proofFile);
  u32 E = 0, power = 0;
  vector<string> knownFactors{};
  char number [2048] = {0};
  char c = 0;
  if (fi.scanf(Proof::HEADER_v2, &power, number, &c) != 3 || c != '\n') {
    log("Proof file '%s' has invalid header\n", proofFile.string().c_str());
    throw "Invalid proof header";
  }
  try {
    E = proof::mersenneFromString(number, knownFactors);
  } catch (string &e) {
    log("Proof file '%s' has invalid header\n", proofFile.string().c_str());
    log(e.c_str());
    throw "Invalid proof header";
  }
  return {power, E, knownFactors, hash};
}

// Mersenne number example: M124647911
// Mersenne cofactor example: M18178631/36357263/145429049/8411216206439
string mersenneToString(u32 E, const vector<string>& knownFactors) {
  string result = "M" + to_string(E);
  for (const auto& factor : knownFactors) {
    result += "/" + factor;
  }
  return result;
}

u32 mersenneFromString(const string& number, vector<string>& knownFactors) {
  if (number.empty() || number[0] != 'M') {
    throw string("Invalid Mersenne number format: must start with M");
  }
  
  // Split by '/' to separate exponent from factors
  vector<string> parts = split(number.substr(1), '/'); // Remove 'M' prefix and split
  
  if (parts.empty()) {
    throw string("Invalid Mersenne number format: no exponent found");
  }
  
  // First part is the exponent
  u32 exponent = 0;
  const char* expStr = parts[0].c_str();
  const char* end = expStr + parts[0].size();
  auto [ptr, ec] = from_chars(expStr, end, exponent, 10);
  
  if (ptr != end || ec != std::errc{}) {
    throw string("Invalid exponent: " + parts[0]);
  }
  
  // Remaining parts are known factors
  knownFactors.clear();
  for (size_t i = 1; i < parts.size(); ++i) {
    if (!parts[i].empty()) {
      // Validate factor is numeric
      try {
        mpz_class test{parts[i]};
        if (test <= 0) {
          throw string("Invalid factor: '" + parts[i] + "' (not positive)");
        }
        knownFactors.push_back(parts[i]);
      } catch (const std::invalid_argument&) {
        throw string("Invalid factor: '" + parts[i] + "' (not numeric)");
      }
    }
  }
  
  return exponent;
}
}

// ---- Proof ----

fs::path Proof::file(const fs::path& proofDir) const {
  string strE = to_string(E);
  u32 power = middles.size();
  return proofDir / (strE + '-' + to_string(power) + ".proof");  
}

void Proof::save(const fs::path& proofFile) const {
  File fo = File::openWrite(proofFile);
  u32 power = middles.size();
  string number = proof::mersenneToString(E, knownFactors);
  fo.printf(HEADER_v2, power, number.c_str(), '\n');
  fo.write(B.data(), (E-1)/8+1);
  for (const Words& w : middles) { fo.write(w.data(), (E-1)/8+1); }
}

Proof Proof::load(const fs::path& path) {
  File fi = File::openReadThrow(path);
  u32 E = 0, power = 0;
  vector<string> knownFactors{};
  char number [2048] = {0};
  char c = 0;
  if (fi.scanf(HEADER_v2, &power, number, &c) != 3 || c != '\n') {
    log("Proof file '%s' has invalid header\n", path.string().c_str());
    throw "Invalid proof header";
  }
  try {
    E = proof::mersenneFromString(number, knownFactors);
  } catch (string &e) {
    log("Proof file '%s' has invalid header\n", path.string().c_str());
    log(e.c_str());
    throw "Invalid proof header";
  }
  u32 nBytes = (E - 1) / 8 + 1;
  Words B = fi.readBytesLE(nBytes);
  vector<Words> middles;
  for (u32 i = 0; i < power; ++i) { middles.push_back(fi.readBytesLE(nBytes)); }
  return {E, knownFactors, B, middles};
}

bool Proof::verify(Gpu *gpu, const vector<u64>& hashes) const {
  // log("B         %016" PRIx64 "\n", res64(B));
  // for (u32 i = 0; i < middles.size(); ++i) { log("Middle[%u] %016" PRIx64 "\n", i, res64(middles[i])); }
  
  u32 power = middles.size();
  assert(power > 0);

  bool isPrime = (B == makeWords(E, 9));

  Words A{makeWords(E, 3)};
  Words B{this->B};
  
  auto hash = proof::hashWords(E, B);

  u32 span = E;
  for (u32 i = 0; i < power; ++i, span = (span + 1) / 2) {
    const Words& M = middles[i];
    hash = proof::hashWords(E, hash, M);
    u64 h = hash[0];
    
    if (hashes.size() > i && h != hashes.at(i)) {
      log("proof [%u] : hash expected %016" PRIx64 " != %016" PRIx64 "\n", i, hashes[i], h);
      return false;
    }

    bool doSquareB = span % 2;
    B = gpu->expMul(M, h, B, doSquareB);
    A = gpu->expMul(A, h, M, false);

    if (gpu->args.verbose) { log("proof [%u] : A %016" PRIx64 ", B %016" PRIx64 ", h %016" PRIx64 "\n", i, res64(A), res64(B), h); }
  }
    
  log("proof verification: doing %d iterations\n", span);
  A = gpu->expExp2(A, span);

  bool ok = (A == B);
  if (ok) {
    log("proof: %u proved %s\n", E, isPrime ? "probable prime" : "composite");
  } else {
    log("proof: invalid (%016" PRIx64 " expected %016" PRIx64 ")\n", res64(A), res64(B));
  }
  return ok;
}

// ---- ProofSet ----

ProofSet::ProofSet(u32 E, const vector<string>& knownFactors, u32 power)
  : E{E}, knownFactors{knownFactors}, power{power} {
  
  assert(E & 1); // E is supposed to be prime
  if (power <= 0 || power > 12) {
    log("Invalid proof power: %u\n", power);
    throw "Invalid proof power";
  }

  fs::create_directories(proofPath(E));

  vector<u32> spans;
  for (u32 span = (E + 1) / 2; spans.size() < power; span = (span + 1) / 2) { spans.push_back(span); }

  points.push_back(0);
  for (u32 p = 0, span = (E + 1) / 2; p < power; ++p, span = (span + 1) / 2) {
    for (u32 i = 0, end = points.size(); i < end; ++i) {
      points.push_back(points[i] + span);
    }
  }

  assert(points.size() == (1u << power));
  assert(points.front() == 0);

  points.front() = E;
  std::sort(points.begin(), points.end());

  assert(points.size() == (1u << power));
  assert(points.back() == E);

  points.push_back(u32(-1)); // guard element
  cacheIt = points.begin();

  for ([[maybe_unused]] u32 p : points) {
    assert(p > E || isInPoints(E, power, p));
  }
}

bool ProofSet::isInPoints(u32 E, u32 power, u32 k) {
  if (k == E) { return true; } // special-case E
  u32 start = 0;
  for (u32 p = 0, span = (E + 1) / 2; p < power; ++p, span = (span + 1) / 2) {
    assert(k >= start);
    if (k > start + span) {
      start += span;
    } else if (k == start + span) {
      return true;
    }
  }
  return false;
}

bool ProofSet::canDo(u32 E, const vector<string>& knownFactors, u32 power, u32 currentK) {
  assert(power > 0 && power <= 12);
  return ProofSet{E, knownFactors, power}.isValidTo(currentK);
}

u32 ProofSet::bestPower(u32 E) {
  // Best proof powers assuming no disk space concern.
  // We increment power by 1 for each fourfold increase of the exponent.
  // The values below produce power=10 at wavefront, and power=11 at 100Mdigits:
  // power=10 from 60M to 240M, power=11 from 240M up.

  assert(E > 0);
  // log2(x)/2 is log4(x)
  int power = 10 + floor(log2(E / 60e6) / 2);
  assert(power >= 2);
  return power;
}

double ProofSet::diskUsageGB(u32 E, u32 power) {
  //  -3 because convert exponent bits to bytes
  // -30 because convert bytes to GB
  // +power because needs 2^power residues for proof generation
  // and 5% on top
  return power ? ldexp(E, -33 + int(power)) * 1.05 : 0.0;
}

u32 ProofSet::effectivePower(u32 E, const vector<string>& knownFactors, u32 power, u32 currentK) {
  for (u32 p = power; p > 0; --p) {
    // log("validating proof residues for power %u\n", p);
    if (canDo(E, knownFactors, p, currentK)) { return p; }
  }
  return 0;
}
    
bool ProofSet::fileExists(u32 k) const {
  return File::size(proofPath(E) / to_string(k)) == i64(E / 32 + 2) * 4;
}

bool ProofSet::isValidTo(u32 limitK) const {
  auto it = upper_bound(points.begin(), points.end(), limitK);

  if (it == points.begin()) {
    return true;
  }

  --it;
  try {
    load(*it);
  } catch (...) {
    return false;
  }

  while (it != points.begin()) {
    if (!fileExists(*--it)) { return false; }
  }
  
  return true;
}

u32 ProofSet::next(u32 k) const {
  if (*cacheIt <= k || (cacheIt > points.begin() && *prev(cacheIt) > k)) {
    cacheIt = upper_bound(points.begin(), points.end(), k);
  }
  return *cacheIt;
}

void ProofSet::save(u32 E, u32 power, u32 k, const Words& words) {
  assert(k && k <= E);
  assert(isInPoints(E, power, k));

  File::openWrite(proofPath(E) / to_string(k)).writeChecked(words);
  assert(load(E, power, k) == words);
}

Words ProofSet::load(u32 E, u32 power, u32 k) {
  assert(k && k <= E);
  assert(isInPoints(E, power, k));
  return File::openReadThrow(proofPath(E) / to_string(k)).readChecked<u32>(E/32 + 1);
}

std::pair<Proof, vector<u64>> ProofSet::computeProof(Gpu *gpu) const {
  Words B = load(E);
  Words A = makeWords(E, 3);

  vector<Words> middles;
  vector<u64> hashes;

  auto hash = proof::hashWords(E, B);

  vector<Buffer<i32>> bufVect = gpu->makeBufVector(power);

  for (u32 p = 0; p < power; ++p) {
    auto bufIt = bufVect.begin();
    assert(p == hashes.size());
    u32 s = (1u << (power - p - 1));
    for (u32 i = 0; i < (1u << p); ++i) {
      Words w = load(points[s * (i * 2 + 1) - 1]);
      gpu->writeIn(*bufIt++, w);
      for (u32 k = 0; i & (1u << k); ++k) {
        assert(k <= p - 1);
        --bufIt;
        u64 h = hashes[p - 1 - k];
        gpu->expMul(*(bufIt - 1), h, *bufIt);
      }
    }
    assert(bufIt == bufVect.begin() + 1);
    Words w = gpu->readAndCompress(bufVect.front());
    if (w.empty()) { throw "Read ZERO during proof generation"; }
    middles.push_back(w);
    hash = proof::hashWords(E, hash, middles.back());
    hashes.push_back(hash[0]);

    log("proof [%u] : M %016" PRIx64 ", h %016" PRIx64 "\n", p, res64(middles.back()), hashes.back());
  }
  return {Proof{E, knownFactors, std::move(B), std::move(middles)}, hashes};
}
