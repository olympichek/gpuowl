#!/usr/bin/env bash
set -euo pipefail
cd /home/ubuntu/gpuowl
FFT="1:512:8:512:202"
BASE="INPLACE=1,LOADS=22042,STORES=21,TABMUL_CHAIN32=1,MODM31=2,MULTI_Q=1,L1CUDA=3,GRAPHS=0"
sudo nvidia-smi -pl 600 >/dev/null
for round in 1 2 3; do
  for arm in C A1 C A2; do
    case $arm in
      C)  USE="$BASE";;
      A1) USE="$BASE,ASYNC_MID61=1";;
      A2) USE="$BASE,ASYNC_MID61=2";;
    esac
    D=$(mktemp -d ".asy${arm}r${round}.XXXXXX")
    ./build-cuda/prpll -dir "$D" -prp 136279841 -iters 100000 -fft "$FFT" -use "$USE" \
      > "$D/run.out" 2>&1 || echo "FAIL $arm r$round"
  done
  echo "== round $round done"
done
echo "== asyncab complete"
