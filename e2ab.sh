#!/usr/bin/env bash
# Production A/B: compiler-axis levers at unlocked clocks, 600 W.
#  arm1: NVRTC 13.3 with PTX version clamped to 9.2 (x3)
#  arm2: driver JIT opt levels 1,2,3 under NVRTC 13.2 (x1 each; expand if promising)
#  anchor: NVRTC 13.2 default (x2, interleaved)
# Gate: 100k residue 52775eea4730be87.
set -euo pipefail
cd /home/ubuntu/gpuowl
FFT="1:512:8:512:202"
USE="INPLACE=1,LOADS=22042,STORES=21,TABMUL_CHAIN32=1,MODM31=2,MULTI_Q=1,L1CUDA=3,GRAPHS=0"
sudo nvidia-smi -pl 600 >/dev/null
sudo nvidia-smi -rgc >/dev/null

run1() {  # tag, then extra env assignments as k=v pairs
  local tag="$1"; shift
  local D
  D=$(mktemp -d ".e2ab-${tag}.XXXXXX")
  env LD_LIBRARY_PATH=/usr/local/cuda-13.2/lib64 "$@" \
    ./build-cuda/prpll -dir "$D" -prp 136279841 -iters 100000 -fft "$FFT" -use "$USE" \
    > "$D/run.out" 2>&1 || echo "FAIL $tag"
}

run1 anchor-a
run1 v133-r1 LD_LIBRARY_PATH=/usr/local/cuda-13.3/lib64 PRPLL_PTX_VERSION=9.2
run1 jito1   PRPLL_JIT_OPT=1
run1 v133-r2 LD_LIBRARY_PATH=/usr/local/cuda-13.3/lib64 PRPLL_PTX_VERSION=9.2
run1 jito2   PRPLL_JIT_OPT=2
run1 v133-r3 LD_LIBRARY_PATH=/usr/local/cuda-13.3/lib64 PRPLL_PTX_VERSION=9.2
run1 jito3   PRPLL_JIT_OPT=3
run1 anchor-b
echo "== e2ab done =="
grep -H "OK    100000" .e2ab-*/run.out || true
