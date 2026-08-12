#!/usr/bin/env bash
# Post-E1 GPU queue:
#  (1) warm reruns of iadd3_lo/hi at lgc2100 (first-run cold-start anomaly check)
#  (2) production prpll draw at lgc2100 / lgc2400 / unlocked-600W (reference row)
set -euo pipefail
cd /home/ubuntu/gpuowl
CC=/tmp/claude-1001/-home-ubuntu-gpuowl/29ee1069-0039-4be4-8ffe-bebd253e2d4d/scratchpad/clockcost
OUT=e1-clockcost
FFT="1:512:8:512:202"
USE="INPLACE=1,LOADS=22042,STORES=21,TABMUL_CHAIN32=1,MODM31=2,MULTI_Q=1,L1CUDA=3,GRAPHS=0"

sudo nvidia-smi -pl 600 >/dev/null
sudo nvidia-smi -lgc 2100,2100 >/dev/null
"$CC" header > "$OUT/warmcheck.csv"
"$CC" lop3_lo 30 >> "$OUT/warmcheck.csv"   # warm the GPU first
for r in iadd3_lo iadd3_hi lop3_lo; do "$CC" "$r" 60 >> "$OUT/warmcheck.csv"; done

prodrun() {
  local tag="$1"
  local D
  D=$(mktemp -d ".e1prod${tag}.XXXXXX")
  nvidia-smi --query-gpu=timestamp,power.draw,clocks.sm,temperature.gpu \
    --format=csv,noheader -lms 500 > "$OUT/prod-$tag-telemetry.csv" &
  local TPID=$!
  LD_LIBRARY_PATH=/usr/local/cuda-13.2/lib64 ./build-cuda/prpll -dir "$D" -prp 136279841 \
    -iters 400000 -fft "$FFT" -use "$USE" > "$D/run.out" 2>&1 || echo "PROD FAIL $tag"
  kill "$TPID" 2>/dev/null || true; wait "$TPID" 2>/dev/null || true
  grep -H "OK    400000\|OK    100000" "$D/run.out" || true
}

prodrun lgc2100
sudo nvidia-smi -lgc 2400,2400 >/dev/null
prodrun lgc2400
sudo nvidia-smi -rgc >/dev/null
prodrun unlocked600
echo "== e1post complete =="
