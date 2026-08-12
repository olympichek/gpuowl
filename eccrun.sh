#!/usr/bin/env bash
# ECC-off measurement: run AFTER ecc.mode.current reads Disabled.
# Baselines (ECC ON, this box): unlocked 600 W 100k tail = 152.45 +- 0.5 us
# (n=15) / 1M tail = 153.3; locked-2100 (meas 2085 MHz) = 171.75 us.
# Model target if ECC explains the +10.1% k gap: locked-2100 -> ~156.5 us.
set -euo pipefail
cd /home/ubuntu/gpuowl
FFT="1:512:8:512:202"
USE="INPLACE=1,LOADS=22042,STORES=21,TABMUL_CHAIN32=1,MODM31=2,MULTI_Q=1,L1CUDA=3,GRAPHS=0"
ECC=$(nvidia-smi --query-gpu=ecc.mode.current --format=csv,noheader | head -1 | tr -d ' ')
[ "$ECC" = "Disabled" ] || { echo "ABORT: ECC current = $ECC (need Disabled)"; exit 1; }
sudo nvidia-smi -pl 600 >/dev/null
for i in 1 2 3; do
  D=$(mktemp -d ".eccoff600r${i}.XXXXXX")
  ./build-cuda/prpll -dir "$D" -prp 136279841 -iters 100000 -fft "$FFT" -use "$USE" \
    > "$D/run.out" 2>&1 || echo "FAIL eccoff600 r$i"
done
sudo nvidia-smi -lgc 2100,2100 >/dev/null
D=$(mktemp -d ".eccoff-lgc2100.XXXXXX")
./build-cuda/prpll -dir "$D" -prp 136279841 -iters 100000 -fft "$FFT" -use "$USE" \
  > "$D/run.out" 2>&1 || echo "FAIL eccoff-lgc"
sudo nvidia-smi -rgc >/dev/null
D=$(mktemp -d ".eccoff1m.XXXXXX")
nvidia-smi --query-gpu=timestamp,power.draw,clocks.sm,temperature.gpu \
  --format=csv,noheader -lms 1000 > "$D/telemetry.csv" &
TPID=$!
./build-cuda/prpll -dir "$D" -prp 136279841 -iters 1000000 -fft "$FFT" -use "$USE" \
  > "$D/run.out" 2>&1 || echo "FAIL eccoff1m"
kill "$TPID" 2>/dev/null || true; wait "$TPID" 2>/dev/null || true
echo "== eccrun complete; analyze .eccoff* dirs (residues: 2k=05d6515c416b83e2, 100k=52775eea4730be87, 1M=52b03a7cc55e677d)"
