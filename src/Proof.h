// Copyright (C) Mihai Preda.

#pragma once

#include "File.h"
#include "common.h"

namespace fs = std::filesystem;

class Gpu;

struct ProofInfo {
  u32 power;
  u32 exp;
  string md5;
};

namespace proof {

array<u64, 4> hashWords(u32 E, const Words& words);

array<u64, 4> hashWords(u32 E, array<u64, 4> prefix, const Words& words);

string fileHash(const fs::path& filePath);

ProofInfo getInfo(const fs::path& proofFile);

}

class Proof {  
public:
  const u32 E;
  const Words B;
  const vector<Words> middles;

  /*Example header:
    PRP PROOF\n
    VERSION=2\n
    HASHSIZE=64\n
    POWER=8\n
    NUMBER=M216091\n
  */
  static const constexpr char* HEADER_v2 = "PRP PROOF\nVERSION=2\nHASHSIZE=64\nPOWER=%u\nNUMBER=M%u%c";

  static Proof load(const fs::path& path);

  void save(const fs::path& proofResultDir) const;

  fs::path file(const fs::path& proofDir) const;
  
  bool verify(Gpu *gpu, const vector<u64>& hashes = {}) const;
};

class ProofSet {
public:
  const u32 E;
  const u32 power;
  const u32 instance;
  
private:  
  vector<u32> points;  
  
  bool isValidTo(u32 limitK) const;

  static bool canDo(u32 E, u32 power, u32 currentK, u32 instance);

  mutable decltype(points)::const_iterator cacheIt{};

  bool fileExists(u32 k) const;

  static fs::path proofPath(u32 E, u32 instance) {
    fs::path worker = "worker-" + to_string(instance);
    return fs::path(worker / to_string(E)) / "proof";
  }
public:
  
  static u32 bestPower(u32 E);
  static u32 effectivePower(u32 E, u32 power, u32 currentK, u32 instance);
  static double diskUsageGB(u32 E, u32 power);
  static bool isInPoints(u32 E, u32 power, u32 k);
  
  ProofSet(u32 E, u32 power, u32 instance);
    
  u32 next(u32 k) const;

  static void save(u32 E, u32 power, u32 k, const Words& words, u32 instance);
  static Words load(u32 E, u32 power, u32 k, u32 instance);
        
  void save(u32 k, const Words& words) const { return save(E, power, k, words, instance); }
  Words load(u32 k) const { return load(E, power, k, instance); }


  std::pair<Proof, vector<u64>> computeProof(Gpu *gpu) const;
};
