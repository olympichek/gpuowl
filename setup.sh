#!/usr/bin/env bash
# PRPLL speedup-campaign bootstrap for a fresh GPU box.
#
# Usage:
#   ./setup.sh [target-dir]          # bootstrap (default target: /workspace/gpuowl)
#   ./setup.sh pack                  # build prpll-campaign-artifacts-YYYYMMDD.zip
#                                    # from the current checkout (for migrating on)
#
# Optional environment overrides:
#   REPO_URL        (default https://github.com/olympichek/gpuowl.git)
#   BRANCH          (default prototype/radix2-width-ownership)
#   SKIP_VALIDATE=1 # skip the 100k-iteration baseline validation run
#
# If a prpll-campaign-artifacts-*.zip sits next to this script, the previous
# box's experiment logs/configs are unpacked into the repo checkout.  The two
# campaign ledgers (sol-prpll-speedup-attempts.md, fable-prpll-speedup-
# attempts.md) travel in git and must be read before starting new experiments.
# If an AI agent continues the campaign, install
# _migration/claude-memory-*.md into its persistent memory.
set -euo pipefail

REPO_URL="${REPO_URL:-https://github.com/olympichek/gpuowl.git}"
BRANCH="${BRANCH:-prototype/radix2-width-ownership}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROD_FFT="1:512:8:512:202"
PROD_USE="INPLACE=1,LOADS=22042,STORES=21,TABMUL_CHAIN32=1,MODM31=2,MULTI_Q=1,L1CUDA=3,GRAPHS=0"

# --- pack mode: bundle this box's experiment artifacts for the next box -----
if [ "${1:-}" = "pack" ]; then
  REPO="$SCRIPT_DIR"
  [ -d "$REPO/.git" ] || { echo "ERROR: run pack from inside the repo checkout"; exit 1; }
  cd "$REPO"
  UNPUSHED=$(git rev-list --count @{upstream}..HEAD 2>/dev/null || echo "?")
  if [ "$UNPUSHED" != 0 ]; then
    echo "WARNING: $UNPUSHED unpushed commit(s) on $(git branch --show-current)."
    echo "         The new box clones from origin — push first or the ledgers will be stale."
  fi
  ZIP="$REPO/prpll-campaign-artifacts-$(date +%Y%m%d).zip"
  rm -f "$ZIP"
  # scratch run dirs (.name.XXXX/): logs/telemetry/configs only — no checkpoint
  # state (.prp), proofs, or kernel caches.
  find . -maxdepth 2 \
      \( -path ./.git -o -name kernel-cache -o -name proof -o -name proof-tmp \
         -o -name '1[0-9]*' \) -prune -o \
      -type f \( -name '*.out' -o -name '*.csv' -o -name 'gpuowl*.log' \
         -o -name 'tune.txt' -o -name 'config.txt' -o -name 'clocks.txt' \) -path './.*' -print \
    | zip -q "$ZIP" -@
  [ -d _migration ] && zip -qr "$ZIP" _migration
  echo "Packed: $ZIP ($(du -h "$ZIP" | cut -f1))"
  echo "Copy setup.sh + this zip to the new box and run ./setup.sh there."
  exit 0
fi

TARGET="${1:-/workspace/gpuowl}"
echo "== PRPLL campaign setup =="

# --- 1. GPU + toolchain sanity --------------------------------------------
command -v nvidia-smi >/dev/null || { echo "ERROR: nvidia-smi not found"; exit 1; }
echo "GPU: $(nvidia-smi --query-gpu=name,power.limit,power.min_limit,power.max_limit,memory.total --format=csv,noheader | head -1)"
if ! command -v nvcc >/dev/null; then
  echo "WARNING: nvcc not on PATH. The build needs the CUDA toolkit (nvrtc + headers),"
  echo "         usually under /usr/local/cuda on Vast base images."
fi
command -v g++  >/dev/null || { echo "ERROR: g++ not found (C++20 required)"; exit 1; }
command -v make >/dev/null || { echo "ERROR: make not found"; exit 1; }

# --- 2. Hardware-control probe ---------------------------------------------
# The 300 W and 600 W campaign boxes had NONE of these; each one that works
# here unlocks a recorded open item (see next-steps below and the ledgers).
echo "--- hardware-control probe ---"
CURPL=$(nvidia-smi --query-gpu=power.limit --format=csv,noheader,nounits | head -1 | cut -d. -f1)
if nvidia-smi -pl "$CURPL" >/dev/null 2>&1; then
  HAVE_PL=1; echo "power-limit control (-pl):  AVAILABLE ($(nvidia-smi --query-gpu=power.min_limit,power.max_limit --format=csv,noheader | head -1))"
else
  HAVE_PL=0; echo "power-limit control (-pl):  blocked (Insufficient Permissions)"
fi
if nvidia-smi -lgc 1500,1500 >/dev/null 2>&1; then
  nvidia-smi -rgc >/dev/null 2>&1 || true
  echo "clock locking (-lgc/-rgc):  AVAILABLE (reset to default)"
else
  echo "clock locking (-lgc/-rgc):  blocked"
fi
NCU="$(command -v ncu || ls /usr/local/cuda/bin/ncu 2>/dev/null || true)"
if [ -n "$NCU" ] && command -v nvcc >/dev/null; then
  TMPD=$(mktemp -d)
  printf '__global__ void k(){}\nint main(){k<<<1,1>>>();cudaDeviceSynchronize();return 0;}\n' > "$TMPD/t.cu"
  if (nvcc -arch=native -o "$TMPD/t" "$TMPD/t.cu" 2>/dev/null || nvcc -o "$TMPD/t" "$TMPD/t.cu" 2>/dev/null) \
     && timeout 180 "$NCU" --set basic "$TMPD/t" > "$TMPD/out" 2>&1 \
     && ! grep -q ERR_NVGPUCTRPERM "$TMPD/out"; then
    echo "NCU counters:               AVAILABLE"
  else
    echo "NCU counters:               blocked (ERR_NVGPUCTRPERM or probe failure)"
  fi
  rm -rf "$TMPD"
else
  echo "NCU counters:               probe skipped (ncu or nvcc missing); $(grep -o 'RmProfilingAdminOnly: [0-9]*' /proc/driver/nvidia/params 2>/dev/null || echo 'RmProfilingAdminOnly: ?')"
fi

# --- 3. Clone / update the repo -------------------------------------------
if [ -d "$TARGET/.git" ]; then
  echo "Repo already present at $TARGET; fetching $BRANCH"
  git -C "$TARGET" fetch origin "$BRANCH"
  git -C "$TARGET" checkout "$BRANCH"
  git -C "$TARGET" pull --ff-only origin "$BRANCH" || true
else
  git clone --branch "$BRANCH" "$REPO_URL" "$TARGET"
fi

# --- 4. Unpack previous-box artifacts (logs/configs), if provided ----------
ZIPFILE="$(ls -1 "$SCRIPT_DIR"/prpll-campaign-artifacts-*.zip 2>/dev/null | head -1 || true)"
if [ -n "$ZIPFILE" ]; then
  echo "Unpacking artifacts: $ZIPFILE"
  unzip -o -q "$ZIPFILE" -d "$TARGET"
else
  echo "No prpll-campaign-artifacts-*.zip next to setup.sh; continuing (ledgers are in git)."
fi

# --- 5. Build ---------------------------------------------------------------
echo "Building CUDA backend..."
make -C "$TARGET" CUDA=1 -j"$(nproc)"
test -x "$TARGET/build-cuda/prpll" || { echo "ERROR: build produced no prpll binary"; exit 1; }

# --- 6. Raise the power limit to the board maximum (default: max throughput) -
MAXPL=$(nvidia-smi --query-gpu=power.max_limit --format=csv,noheader,nounits | head -1 | cut -d. -f1)
if [ "$CURPL" -lt "$MAXPL" ] && [ "$HAVE_PL" = 1 ]; then
  nvidia-smi -pl "$MAXPL" >/dev/null 2>&1 && echo "Power limit raised: ${CURPL} W -> ${MAXPL} W"
else
  echo "Power limit: ${CURPL} W (board max ${MAXPL} W)"
fi

# --- 7. Baseline validation: 100k iterations, exact residues ----------------
if [ "${SKIP_VALIDATE:-0}" != "1" ]; then
  RUNDIR="$TARGET/.setup-validate"
  rm -rf "$RUNDIR"; mkdir -p "$RUNDIR"
  echo "Running 100k-iteration production validation (under a minute on a 600 W GB202)..."
  nvidia-smi --query-gpu=power.draw,clocks.sm,temperature.gpu,clocks_event_reasons.active --format=csv,noheader -l 1 \
    > "$RUNDIR/telemetry.csv" 2>/dev/null &
  SMIPID=$!
  ( cd "$RUNDIR" && "$TARGET/build-cuda/prpll" -dir . -prp 136279841 -iters 100000 \
      -fft "$PROD_FFT" -use "$PROD_USE" ) > "$RUNDIR/run.out" 2>&1 || true
  kill "$SMIPID" 2>/dev/null || true
  LOG="$RUNDIR/gpuowl-0.log"
  echo "--- validation summary ---"
  grep -E "FFT:|OK |00000 " "$LOG" | tail -8 || true
  # Log line format: date time exponent OK iter residue ...
  R2K=$(grep -E "OK +2000 " "$LOG" | awk '{print $6}' | head -1 || true)
  R100K=$(grep -E "OK +100000 " "$LOG" | awk '{print $6}' | head -1 || true)
  FAIL=0
  [ "$R2K" = "05d6515c416b83e2" ]   || { echo "RESIDUE MISMATCH at 2k:   got '${R2K:-none}' want 05d6515c416b83e2"; FAIL=1; }
  [ "$R100K" = "52775eea4730be87" ] || { echo "RESIDUE MISMATCH at 100k: got '${R100K:-none}' want 52775eea4730be87"; FAIL=1; }
  [ "$FAIL" = 0 ] && echo "Residues at 2k and 100k match the recorded production values."
  echo "--- under-load telemetry (samples above 100 W) ---"
  awk -F', ' '$1+0 > 100 {pw+=$1; ck+=$2+0; if ($2+0<mn||!n) mn=$2+0; if ($2+0>mx) mx=$2+0; n++} END {
    if (n) printf "samples=%d  mean power=%.0f W  SM clock mean=%d MHz (min %d, max %d)\n", n, pw/n, ck/n, mn, mx
    else print "no under-load samples captured (run too short for 1-s sampling)"
  }' "$RUNDIR/telemetry.csv"
  echo "References: 300 W Max-Q 201.5 us @1552 MHz; 600 W WSE 137-138 us cool /"
  echo "  148.9 us thermal-steady @2210 MHz.  Model: t(f) ~ 6.0 us + 315842/f_MHz."
  echo "  (5M-iteration reference residue: 7eca65291732df02.)"
  [ "$FAIL" = 0 ] || exit 1
fi

cat <<'EON'

== Setup complete ==
Read the campaign ledgers BEFORE any new experiment:
    sol-prpll-speedup-attempts.md    - registry of ~95 rejected designs
    fable-prpll-speedup-attempts.md  - 600 W phase, scaling model, gate closures

The software/tune space is closed (see "Updated campaign boundary (600 W)").
What THIS box can newly unlock depends on the hardware-control probe above:

  1. power-limit control (-pl): run the direct power sweep the 300/600 W boxes
     could not: for PL in 150..max step 50, fresh-dir 200k+ production run to
     thermal steady state with 1-s telemetry; record us/iter, sustained clock,
     us*W (energy/iter).  Deliverables: perf(PL) curve vs the P^0.49 model and
     the perf/W knee (fleet operating point).  Then map two-worker aggregate
     (-prps 136279841,136279879 -workers 2) across the same PLs to find its
     break-even (needs f_2w/f_1w >= 0.914; it was 0.840 at 600 W).
  2. clock locking (-lgc): validate t(f) = 6.0 + 315842/f directly at fixed
     clocks (decouples the V/f governor); measure iso-clock power of one vs
     two workers to separate work-density from voltage effects.
  3. NCU counters: profile the M61 tail/middle kernels for stall reasons
     during production co-run — the last open software route (an invisible
     intra-kernel stall).  Start: ncu --set full on tailSquareGF61 /
     fftMiddleOutGF61 with production -use flags.
Record everything in fable-prpll-speedup-attempts.md (matched alternating
runs, residues at every checkpoint; commit with "Record ..." messages).
EON
