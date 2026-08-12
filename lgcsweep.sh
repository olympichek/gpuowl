#!/usr/bin/env bash
# t(f) model on this box: -lgc locked-clock ladder at 600 W, 100k per point.
set -euo pipefail
cd /home/ubuntu/gpuowl
FFT="1:512:8:512:202"
USE="INPLACE=1,LOADS=22042,STORES=21,TABMUL_CHAIN32=1,MODM31=2,MULTI_Q=1,L1CUDA=3,GRAPHS=0"
sudo nvidia-smi -pl 600 >/dev/null
for F in 1500 1700 1900 2100 2300; do
  sudo nvidia-smi -lgc "$F,$F" >/dev/null
  sleep 2
  D=$(mktemp -d ".lgc${F}.XXXXXX")
  nvidia-smi --query-gpu=timestamp,power.draw,clocks.sm,temperature.gpu,clocks_event_reasons.active \
    --format=csv,noheader -lms 1000 > "$D/telemetry.csv" &
  TPID=$!
  ./build-cuda/prpll -dir "$D" -prp 136279841 -iters 100000 -fft "$FFT" -use "$USE" \
    > "$D/run.out" 2>&1 || echo "FAIL lgc $F"
  kill "$TPID" 2>/dev/null || true; wait "$TPID" 2>/dev/null || true
  echo "== done lgc $F: $D"
done
sudo nvidia-smi -rgc >/dev/null
echo "== lgcsweep complete"
