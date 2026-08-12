#!/usr/bin/env bash
# 1) -lmc 405 probe: the only non-idle alternative memory state; expect a
#    large loss (middles are DRAM-latency-exposed); closes the -lmc knob.
# 2) TAIL_TRIGS61 generate-vs-load: audit-derived cheap test.  Default 0
#    reads all M61 tail trigs from memory; 1/2 compute progressively more
#    from scratch.  Alternating fresh 100k runs at 600 W.
set -euo pipefail
cd /home/ubuntu/gpuowl
FFT="1:512:8:512:202"
BASE="INPLACE=1,LOADS=22042,STORES=21,TABMUL_CHAIN32=1,MODM31=2,MULTI_Q=1,L1CUDA=3,GRAPHS=0"
sudo nvidia-smi -pl 600 >/dev/null

D=$(mktemp -d ".lmc405.XXXXXX")
sudo nvidia-smi -lmc 405,405 >/dev/null
sleep 2
nvidia-smi --query-gpu=timestamp,power.draw,clocks.sm,clocks.mem,temperature.gpu \
  --format=csv,noheader -lms 1000 > "$D/telemetry.csv" &
TPID=$!
./build-cuda/prpll -dir "$D" -prp 136279841 -iters 20000 -fft "$FFT" -use "$BASE" \
  > "$D/run.out" 2>&1 || echo "FAIL lmc405"
kill "$TPID" 2>/dev/null || true; wait "$TPID" 2>/dev/null || true
sudo nvidia-smi -rmc >/dev/null 2>&1 || sudo nvidia-smi -lmc 12481,12481 >/dev/null 2>&1 || true
echo "== done lmc405: $D"

for round in 1 2 3; do
  for arm in C T1 T2; do
    case $arm in
      C)  USE="$BASE";;
      T1) USE="$BASE,TAIL_TRIGS61=1";;
      T2) USE="$BASE,TAIL_TRIGS61=2";;
    esac
    D=$(mktemp -d ".trig${arm}r${round}.XXXXXX")
    ./build-cuda/prpll -dir "$D" -prp 136279841 -iters 100000 -fft "$FFT" -use "$USE" \
      > "$D/run.out" 2>&1 || echo "FAIL trig $arm r$round"
  done
  echo "== trig round $round done"
done
echo "== lmctrig complete"
