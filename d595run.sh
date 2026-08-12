#!/usr/bin/env bash
# Driver-595.84 perf battery.  ECC-off/580 baselines: 100k@600W tail
# 150.2-151.1 us; locked-2100 171.7; ladder fit c=4.98 k=347,700.
set -euo pipefail
cd /home/ubuntu/gpuowl
FFT="1:512:8:512:202"
USE="INPLACE=1,LOADS=22042,STORES=21,TABMUL_CHAIN32=1,MODM31=2,MULTI_Q=1,L1CUDA=3,GRAPHS=0"
sudo nvidia-smi -pl 600 >/dev/null
for i in 1 2 3; do
  D=$(mktemp -d ".d595r${i}.XXXXXX")
  ./build-cuda/prpll -dir "$D" -prp 136279841 -iters 100000 -fft "$FFT" -use "$USE" \
    > "$D/run.out" 2>&1 || echo "FAIL d595 r$i"
done
for F in 2100 1500 1700 1900 2300; do
  sudo nvidia-smi -lgc "$F,$F" >/dev/null; sleep 2
  D=$(mktemp -d ".d595lgc${F}.XXXXXX")
  nvidia-smi --query-gpu=timestamp,power.draw,clocks.sm,temperature.gpu \
    --format=csv,noheader -lms 1000 > "$D/telemetry.csv" &
  TPID=$!
  ./build-cuda/prpll -dir "$D" -prp 136279841 -iters 100000 -fft "$FFT" -use "$USE" \
    > "$D/run.out" 2>&1 || echo "FAIL d595 lgc$F"
  kill "$TPID" 2>/dev/null || true; wait "$TPID" 2>/dev/null || true
done
sudo nvidia-smi -rgc >/dev/null
echo "== d595run complete"
