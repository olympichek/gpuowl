#!/usr/bin/env bash
# E2 capture: per-NVRTC-version PTX (prpll kernel-cache) + final SASS (driver JIT cache).
# Runs 2000 iters per version (residue gate 05d6515c416b83e2), then carves cubins
# out of the JIT cache and nvdisasm's them.
set -euo pipefail
cd /home/ubuntu/gpuowl
FFT="1:512:8:512:202"
USE="INPLACE=1,LOADS=22042,STORES=21,TABMUL_CHAIN32=1,MODM31=2,MULTI_Q=1,L1CUDA=3,GRAPHS=0"
OUT=e2-sass
mkdir -p "$OUT"

for V in 13.0 13.2 13.3; do
  TAG=${V/./}
  D=$(mktemp -d ".e2cap${TAG}.XXXXXX")
  JC="$PWD/$OUT/jit-$TAG"
  rm -rf "$JC"; mkdir -p "$JC"
  EXTRA=()
  [ "$V" = "13.3" ] && EXTRA=(PRPLL_PTX_VERSION=9.2)
  env LD_LIBRARY_PATH=/usr/local/cuda-$V/lib64 CUDA_CACHE_PATH="$JC" "${EXTRA[@]}" \
    ./build-cuda/prpll -dir "$D" -prp 136279841 -iters 2000 -fft "$FFT" -use "$USE" -cache \
    > "$D/run.out" 2>&1 || echo "RUN FAIL $V"
  grep -H "OK      2000" "$D/run.out" || echo "GATE MISSING $V ($D)"
  mkdir -p "$OUT/ptx-$TAG"
  cp "$D"/kernel-cache/* "$OUT/ptx-$TAG/" 2>/dev/null || echo "no kernel-cache for $V"
done

python3 e2carve.py "$OUT"
echo "== E2 capture complete =="
