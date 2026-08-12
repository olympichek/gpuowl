#!/usr/bin/env bash
# NVRTC 13.2 vs 13.3 A/B: alternating fresh-dir 100k runs, unlocked clocks, 600 W.
# Gates: 2k=05d6515c416b83e2, 100k=52775eea4730be87.
set -euo pipefail
cd /home/ubuntu/gpuowl
FFT="1:512:8:512:202"
USE="INPLACE=1,LOADS=22042,STORES=21,TABMUL_CHAIN32=1,MODM31=2,MULTI_Q=1,L1CUDA=3,GRAPHS=0"
sudo nvidia-smi -pl 600 >/dev/null
sudo nvidia-smi -rgc >/dev/null
for r in 1 2 3; do
  for v in 13.2 13.3; do
    D=$(mktemp -d ".nv${v/./}r${r}.XXXXXX")
    LD_LIBRARY_PATH=/usr/local/cuda-$v/lib64 ./build-cuda/prpll -dir "$D" -prp 136279841 \
      -iters 100000 -fft "$FFT" -use "$USE" > "$D/run.out" 2>&1 || echo "FAIL $v r$r"
  done
done
echo "== A/B done =="
grep -H "OK    100000" .nv13*/run.out || true
