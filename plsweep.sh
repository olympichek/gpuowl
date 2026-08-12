#!/usr/bin/env bash
# Direct perf/W sweep: nvidia-smi -pl 600..300 W, production config, 1M iters
# per point, 1-s telemetry.  Final 600b point re-checks drift/hysteresis.
set -euo pipefail
cd /home/ubuntu/gpuowl
PROD_FFT="1:512:8:512:202"
PROD_USE="INPLACE=1,LOADS=22042,STORES=21,TABMUL_CHAIN32=1,MODM31=2,MULTI_Q=1,L1CUDA=3,GRAPHS=0"
for W in 600 550 500 450 400 350 300 600b; do
  PL=${W%b}
  sudo nvidia-smi -pl "$PL" >/dev/null
  sleep 3
  D=$(mktemp -d ".plsw${W}.XXXXXX")
  cp -r .setup-validate/kernel-cache "$D/" 2>/dev/null || true
  nvidia-smi --query-gpu=timestamp,power.draw,clocks.sm,temperature.gpu,clocks_event_reasons.active \
    --format=csv,noheader -lms 1000 > "$D/telemetry.csv" &
  TPID=$!
  ./build-cuda/prpll -dir "$D" -prp 136279841 -iters 1000000 -fft "$PROD_FFT" -use "$PROD_USE" \
    > "$D/run.out" 2>&1 || echo "WARN: prpll exit $? at $W"
  kill "$TPID" 2>/dev/null || true
  wait "$TPID" 2>/dev/null || true
  echo "== done $W: $D"
done
sudo nvidia-smi -pl 600 >/dev/null
echo "== sweep complete"
