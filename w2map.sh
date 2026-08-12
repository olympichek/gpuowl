#!/usr/bin/env bash
# Two-worker break-even map across power points.  1-worker baselines come
# from the .plsw* sweep at the same limits.  200k/exponent per point.
set -euo pipefail
cd /home/ubuntu/gpuowl
FFT="1:512:8:512:202"
USE="INPLACE=1,LOADS=22042,STORES=21,TABMUL_CHAIN32=1,MODM31=2,MULTI_Q=1,L1CUDA=3,GRAPHS=0"
for PL in 600 500 400 300; do
  sudo nvidia-smi -pl "$PL" >/dev/null
  sleep 3
  D=$(mktemp -d ".w2pl${PL}.XXXXXX")
  nvidia-smi --query-gpu=timestamp,power.draw,clocks.sm,temperature.gpu,clocks_event_reasons.active \
    --format=csv,noheader -lms 1000 > "$D/telemetry.csv" &
  TPID=$!
  ./build-cuda/prpll -dir "$D" -prps 136279841,136279879 -workers 2 -iters 200000 \
    -fft "$FFT" -use "$USE" > "$D/run.out" 2>&1 || echo "FAIL w2 at $PL"
  kill "$TPID" 2>/dev/null || true; wait "$TPID" 2>/dev/null || true
  echo "== done w2 $PL: $D"
done
sudo nvidia-smi -pl 600 >/dev/null
echo "== w2map complete"
