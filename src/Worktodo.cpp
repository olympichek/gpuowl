// Copyright (C) Mihai Preda.

#include "Worktodo.h"

#include "Task.h"
#include "File.h"
#include "common.h"
#include "Args.h"
#include "fs.h"

#include <cassert>
#include <string>
#include <optional>
#include <charconv>
#include <gmpxx.h>

namespace {

bool isHex(const string& s) {
  u32 dummy{};
  const char *end = s.c_str() + s.size();
  auto[ptr, ec] = std::from_chars(s.c_str(), end, dummy, 16);
  return (ptr == end);
}

// Split string respecting quoted sections (for PRP-CF assignment parsing)
vector<string> splitRespectingQuotes(const string& s, char delim) {
  vector<string> result;
  string current;
  bool inQuotes = false;
  
  for (char c : s) {
    if (c == '"') {
      inQuotes = !inQuotes;
      current += c;
    } else if (c == delim && !inQuotes) {
      result.push_back(current);
      current.clear();
    } else {
      current += c;
    }
  }
  if (!current.empty()) {
    result.push_back(current);
  }
  return result;
}

// Validate that factors are actually factors of the Mersenne number 2^p - 1
bool validateMersenneFactors(u64 exponent, const std::vector<string>& factors) {
  try {
    // Compute 2^p - 1
    mpz_class mersenne = (mpz_class{1} << exponent) - 1;
    
    // Check each factor divides the Mersenne number
    for (const auto& factorStr : factors) {
      if (factorStr.empty()) continue;
      
      mpz_class factor{factorStr};
      if (factor <= 1) {
        log("Failed to parse PRP-CF assignment: factor '%s' must be > 1\n", factorStr.c_str());
        return false;
      }
      
      if (mersenne % factor != 0) {
        log("Failed to parse PRP-CF assignment: '%s' is not a factor of M%lu\n", factorStr.c_str(), exponent);
        return false;
      }
    }
    
    return true;
    
  } catch (const std::exception& e) {
    log("Failed to parse PRP-CF assignment: %s\n", e.what());
    return false;
  }
}

// Parse comma-separated factors from quoted string like "36357263,145429049,8411216206439"
std::vector<string> parseFactors(const string& factorStr) {
  string trimmed = rstripNewline(factorStr);
  
  std::vector<string> factors;
  if (trimmed.size() >= 2 && trimmed.front() == '"' && trimmed.back() == '"') {
    string content = trimmed.substr(1, trimmed.size() - 2); // Remove quotes
    factors = split(content, ',');
    // Validate factors are numeric (using GMP for large factor support)
    for (const auto& factor : factors) {
      if (factor.empty()) continue;
      try {
        mpz_class test{factor};
        if (test <= 0) {
          log("parseFactors: invalid factor '%s' (not positive)\n", factor.c_str());
          return {}; // Factor must be positive
        }
      } catch (const std::invalid_argument&) {
        log("parseFactors: invalid factor '%s' (not numeric)\n", factor.c_str());
        return {}; // Invalid factor format
      }
    }
  }
  return factors;
}

// Examples:
// PRP=FEEE9DCD59A0855711265C1165C4C693,1,2,124647911,-1,77,0
// PRP=D01D05DD3394CFF8887960999DC0D9EE,1,2,18178631,-1,99,2,"36357263,145429049,8411216206439"
// DoubleCheck=E0F583710728343C61643028FBDBA0FB,70198703,75,1
// Cert=B2EE67DC0A514753E488794C9DD6F6BD,1,2,82997591,-1,162105
std::optional<Task> parse(const std::string& line) {
  if (line.empty() || line[0] == '#') { return {}; }

  vector<string> topParts = split(line, '=');

  bool isPRP = false;
  bool isLL = false;
  bool isCERT = false;

  if (topParts.size() == 2) {
    string kind = topParts.front();
    if (kind == "PRP" || kind == "PRPDC") {
      isPRP = true;
    } else if (kind == "Test" || kind == "DoubleCheck") {
      isLL = true;
    } else if (kind == "Cert") {
      isCERT = true;
    }
  }

  if (isPRP || isLL) {
    vector<string> parts = splitRespectingQuotes(topParts.back(), ',');
    if (!parts.empty() && (parts.front() == "N/A" || parts.front().empty())) {
      parts.erase(parts.begin()); // skip empty AID
    }

    string AID;
    if (!parts.empty() && parts.front().size() == 32 && isHex(parts.front())) {
      AID = parts.front();
      parts.erase(parts.begin());
    }

    string s = (parts.size() >= 4 && parts[0] == "1" && parts[1] == "2" && parts[3] == "-1") ? parts[2]
      : (!parts.empty() ? parts[0] : "");

    const char *end = s.c_str() + s.size();
    u64 exp{};
    auto [ptr, _] = from_chars(s.c_str(), end, exp, 10);
    if (ptr != end) { exp = 0; }
    
    if (exp > 1000) {
      Task task{isPRP ? Task::PRP : Task::LL, u32(exp), AID, line, 0};
      
      // Check for cofactor format: PRP=AID,1,2,exponent,-1,how_far_factored,tests_saved,"factors"
      if (isPRP && parts.size() >= 7) {
        const string& lastPart = parts.back();
        auto factors = parseFactors(lastPart);
        if (!factors.empty() && validateMersenneFactors(exp, factors)) {
          task.knownFactors = factors;
          task.residueType = 5; // Type 5 for cofactor tests
        }
        else {
          log("Skipping PRP-CF assignment with invalid factors: \"%s\"\n", rstripNewline(line).c_str());
          return {}; // Skip this assignment
        }
      }
      
      return task;
    }
  }
  if (isCERT) {
    vector<string> parts = split(topParts.back(), ',');
    if (!parts.empty() && parts.front().size() == 32 && isHex(parts.front())) {
      string AID;
      AID = parts.front();
      parts.erase(parts.begin());

      if (parts.size() == 5 && parts[0] == "1" && parts[1] == "2" && parts[3] == "-1") {
	string s = parts[2];
	const char *end = s.c_str() + s.size();
	u64 exp{0};
	from_chars(s.c_str(), end, exp, 10);
	s = parts[4];
	end = s.c_str() + s.size();
	u64 squarings{0};
	from_chars(s.c_str(), end, squarings, 10);
//printf ("Exec cert %d %d \n", (int) exp, (int) squarings);
	if (exp > 1000 && squarings > 100) { return {{Task::CERT, u32(exp), AID, line, u32(squarings) }}; }
      }
    }
  }
  log("worktodo.txt line ignored: \"%s\"\n", rstripNewline(line).c_str());
  return {};
}

// Among the valid tasks from fileName, return the "best" which means the smallest CERT, or otherwise the exponent PRP/LL
static std::optional<Task> bestTask(const fs::path& fileName) {
  optional<Task> best;
  for (const string& line : File::openRead(fileName)) {
    optional<Task> task = parse(line);
    if (task && (!best
                 || (best->kind != Task::CERT && task->kind == Task::CERT)
                 || ((best->kind != Task::CERT || task->kind == Task::CERT) && task->exponent < best->exponent))) {
      best = task;
    }
  }
  return best;
}

string workName(i32 instance) { return "worktodo-" + to_string(instance) + ".txt"; }

optional<Task> getWork(Args& args, i32 instance) {
  fs::path localWork = workName(instance);

  // Try to get a task from the local worktodo-<N> file.
  if (optional<Task> task = bestTask(localWork)) { return task; }

  if (args.masterDir.empty()) { return {}; }

  fs::path worktodo = args.masterDir / "worktodo.txt";

  /*
    We need to aquire a task from the global worktodo.txt, and "atomically"
    add the task to the local worktodo-N.txt and remove it from worktodo.txt

    Below we call the global worktodo.txt "global worktodo", and worktodo-N.txt "local worktodo".

    We want to avoid filesystem-based locking, so we approximate it this way:
    1. read the file-size of the global worktodo
    2. read one task from the global worktodo
    3. append the task to the local worktodo
    4. write the new content of the global worktodo without the task to a temporary file
    5. compare the size of the global worktodo with its initial size (as an heuristic to detect modifications to it)
    6a. if the size is not changed, rename the temporary file to global worktodo and done
    6b. if the size is changed (i.e. global worktodo was modified in the meantime):
       7. remove the task from the local worktodo (undo the local task add)
       8. start again (from step 1)
  */

  for (int retry = 0; retry < 2; ++retry) {
    u64 initialSize = fileSize(worktodo);
    if (!initialSize) { return {}; }

    optional<Task> task = bestTask(worktodo);
    if (!task) { return {}; }

    string workLine = task->line;
    File::append(localWork, workLine);

    if (deleteLine(worktodo, workLine, initialSize)) {
      return task;
    }

    // Undo add to local worktodo. Attempt twice.
    bool found = deleteLine(localWork, workLine) || deleteLine(localWork, workLine);
    assert(found);
    if (!found) { return {}; }
  }

  log("Could not extract a task from '%s'\n", worktodo.string().c_str());
  // must be tough luck to be preempted twice while mutating the global worktodo
  assert(false);
  return {};
}

} // namespace

std::optional<Task> Worktodo::getTask(Args &args, i32 instance) {
  if (instance == 0) {
    if (args.prpExp) {
      u32 exp = args.prpExp;
      args.prpExp = 0;
      return Task{Task::PRP, exp};
    } else if (args.llExp) {
      u32 exp = args.llExp;
      args.llExp = 0;
      return Task{Task::LL, exp};
    } else if (!args.verifyPath.empty()) {
      auto path = args.verifyPath;
      args.verifyPath.clear();
      return Task{.kind=Task::VERIFY, .verifyPath=path};
    }
  }
  return getWork(args, instance);
}

bool Worktodo::deleteTask(const Task &task, i32 instance) {
  // Some tasks don't originate in worktodo.txt and thus don't need deleting.
  if (task.line.empty()) { return true; }
  return deleteLine(workName(instance), task.line);
}
