#!/usr/bin/env bash
# PRPLL speedup-campaign bootstrap for a fresh GPU box (e.g. 600 W RTX PRO 6000).
#
# Usage:
#   ./setup.sh [target-dir]          # default target: /workspace/gpuowl
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
set -euo pipefail

REPO_URL="${REPO_URL:-https://github.com/olympichek/gpuowl.git}"
BRANCH="${BRANCH:-prototype/radix2-width-ownership}"
TARGET="${1:-/workspace/gpuowl}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROD_FFT="1:512:8:512:202"
PROD_USE="INPLACE=1,LOADS=22042,STORES=21,TABMUL_CHAIN32=1,MODM31=2,MULTI_Q=1,L1CUDA=3,GRAPHS=0"

echo "== PRPLL campaign setup =="

# --- 1. GPU + toolchain sanity --------------------------------------------
command -v nvidia-smi >/dev/null || { echo "ERROR: nvidia-smi not found"; exit 1; }
echo "GPU: $(nvidia-smi --query-gpu=name,power.limit,power.max_limit,memory.total --format=csv,noheader | head -1)"
if ! command -v nvcc >/dev/null; then
  echo "WARNING: nvcc not on PATH. The build needs the CUDA toolkit (nvrtc + headers),"
  echo "         usually under /usr/local/cuda on Vast base images."
fi
command -v g++  >/dev/null || { echo "ERROR: g++ not found (C++20 required)"; exit 1; }
command -v make >/dev/null || { echo "ERROR: make not found"; exit 1; }

# --- 2. Clone / update the repo -------------------------------------------
if [ -d "$TARGET/.git" ]; then
  echo "Repo already present at $TARGET; fetching $BRANCH"
  git -C "$TARGET" fetch origin "$BRANCH"
  git -C "$TARGET" checkout "$BRANCH"
  git -C "$TARGET" pull --ff-only origin "$BRANCH" || true
else
  git clone --branch "$BRANCH" "$REPO_URL" "$TARGET"
fi

# --- 3. Unpack previous-box artifacts (logs/configs), if provided ----------
ZIPFILE="$(ls -1 "$SCRIPT_DIR"/prpll-campaign-artifacts-*.zip 2>/dev/null | head -1 || true)"
if [ -n "$ZIPFILE" ]; then
  echo "Unpacking artifacts: $ZIPFILE"
  unzip -o -q "$ZIPFILE" -d "$TARGET"
else
  echo "No prpll-campaign-artifacts-*.zip next to setup.sh; continuing (ledgers are in git)."
fi

# --- 4. Build ---------------------------------------------------------------
echo "Building CUDA backend..."
make -C "$TARGET" CUDA=1 -j"$(nproc)"
test -x "$TARGET/build-cuda/prpll" || { echo "ERROR: build produced no prpll binary"; exit 1; }

# --- 5. Try to raise the power limit to the board maximum -------------------
MAXPL=$(nvidia-smi --query-gpu=power.max_limit --format=csv,noheader,nounits | head -1 | cut -d. -f1)
CURPL=$(nvidia-smi --query-gpu=power.limit     --format=csv,noheader,nounits | head -1 | cut -d. -f1)
if [ "$CURPL" -lt "$MAXPL" ]; then
  if nvidia-smi -pl "$MAXPL" >/dev/null 2>&1; then
    echo "Power limit raised: ${CURPL} W -> ${MAXPL} W"
  else
    echo "NOTE: cannot raise power limit (${CURPL}/${MAXPL} W) in this container; runs stay at ${CURPL} W."
  fi
else
  echo "Power limit already at board maximum: ${CURPL} W"
fi

# --- 6. Baseline validation: 100k iterations, exact residues ----------------
if [ "${SKIP_VALIDATE:-0}" != "1" ]; then
  RUNDIR="$TARGET/.setup-validate"
  rm -rf "$RUNDIR"; mkdir -p "$RUNDIR"
  echo "Running 100k-iteration production validation (a few minutes)..."
  nvidia-smi --query-gpu=power.draw,clocks.sm,temperature.gpu --format=csv,noheader -l 2 \
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
  echo "--- under-load telemetry (samples above 150 W) ---"
  awk -F', ' '$1+0 > 150 {pw+=$1; ck+=$2+0; if ($2+0<mn||!n) mn=$2+0; if ($2+0>mx) mx=$2+0; n++} END {
    if (n) printf "samples=%d  mean power=%.0f W  SM clock mean=%d MHz (min %d, max %d)\n", n, pw/n, ck/n, mn, mx
    else print "no under-load samples captured (run too short for 2-s sampling)"
  }' "$RUNDIR/telemetry.csv"
  echo "Reference on the Max-Q 300 W box: ~201.5 us/iteration at ~1552 MHz sustained."
  [ "$FAIL" = 0 ] || exit 1
fi

cat <<'EON'

== Setup complete ==
Next steps on a higher-power board:
  1. Compare the us/iteration above against the 201.5-us Max-Q reference and
     note the sustained SM clock from .setup-validate/telemetry.csv.
  2. Re-tune: the incumbent -use line was tuned at 300 W / 1552 MHz and the
     optimum may shift at higher clocks:
        ./build-cuda/prpll -dir <fresh-dir> -prp 136279841 -tune ntt
     Verify any suggested change with matched alternating 50k runs (fresh
     directory per run; residues must match at every 10k checkpoint:
     52316d51aa52e6b7 @10k, 6a5c8b8989125413 @20k, 9139db3046e846d4 @30k,
     b597cca6031938ff @40k, cce5a14b7c17aecb @50k).
  3. Read the campaign ledgers BEFORE attempting new speedup ideas:
        sol-prpll-speedup-attempts.md    - registry of ~95 rejected designs
        fable-prpll-speedup-attempts.md  - saturation model, latest state
EON
