// Copyright (C) Mihai Preda

#pragma once

#include "Args.h"
#include "common.h"
#include "GpuCommon.h"

#include <string>
#include <vector>

class Args;
class Result;
class Context;
class Queue;
class TrigBufCache;

class Task {
public:
  enum Kind {PRP, VERIFY, LL, CERT};

  Kind kind;
  u32 exponent;
  string AID;  // Assignment ID
  string line; // the verbatim worktodo line, used in deleteTask().
  u32 squarings;  // For CERTs

  std::vector<string> knownFactors; // For PRP on cofactors
  int residueType = 1;  // Default Type 1, Type 5 for cofactors
  bool isCofactor() const { return !knownFactors.empty(); }

  string verifyPath; // For Verify
  void execute(GpuCommon shared, Queue* q, u32 instance);

  void writeResultPRP(const Args&, bool isPrime, u64 res64, const std::string& res2048, u32 fftSize, u32 nErrors, const fs::path& proofPath) const;
  void writeResultLL(const Args&, bool isPrime, u64 res64, u32 fftSize) const;
  void writeResultCERT(const Args&, array <u64, 4> hash, u32 squarings, u32 fftSize) const;
};
