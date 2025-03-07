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

namespace {

bool isHex(const string& s) {
  u32 dummy{};
  const char *end = s.c_str() + s.size();
  auto[ptr, ec] = std::from_chars(s.c_str(), end, dummy, 16);
  return (ptr == end);
}

// Examples:
// PRP=FEEE9DCD59A0855711265C1165C4C693,1,2,124647911,-1,77,0
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
    vector<string> parts = split(topParts.back(), ',');
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
    if (exp > 1000) { return {{isPRP ? Task::PRP : Task::LL, u32(exp), AID, line, 0}}; }
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
