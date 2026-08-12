#!/usr/bin/env bash
set -euo pipefail
cd /home/ubuntu/gpuowl
FFT="1:512:8:512:202"
B="INPLACE=1,TABMUL_CHAIN32=1,MODM31=2,MULTI_Q=1,L1CUDA=3,GRAPHS=0"
declare -A ARMS=(
  [C]="$B,LOADS=22042,STORES=21"
  [L1]="$B,LOADS=20042,STORES=21"
  [L2]="$B,LOADS=10042,STORES=21"
  [S1]="$B,LOADS=22042,STORES=20"
  [S2]="$B,LOADS=22042,STORES=24"
  [MT]="$B,LOADS=22042,STORES=21,MIDDLE_IN_LDS_TRANSPOSE=1,MIDDLE_OUT_LDS_TRANSPOSE=1"
)
sudo nvidia-smi -pl 600 >/dev/null
for round in 1 2 3; do
  for a in C L1 C L2 C S1 C S2 C MT; do
    D=$(mktemp -d ".ms${a}r${round}.XXXXXX")
    ./build-cuda/prpll -dir "$D" -prp 136279841 -iters 100000 -fft "$FFT" \
      -use "${ARMS[$a]}" > "$D/run.out" 2>&1 || echo "FAIL $a r$round"
  done
  echo "== round $round"
done
echo "== microsweep complete"
