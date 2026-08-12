#!/usr/bin/env bash
# NCU-guided occupancy experiment: per-kernel register caps (and one
# prefer-shared arm) vs production, alternating fresh 100k runs at 600 W.
# Counters basis: carryFused/tailSquareGF61 reg-limited to 20 warps/SM with
# long_scoreboard (DRAM latency) the top stall; GF61 middles reg-limited to
# 32 warps with long_sb 22-24 cycles/issue.
set -euo pipefail
cd /home/ubuntu/gpuowl
FFT="1:512:8:512:202"
BASE="INPLACE=1,LOADS=22042,STORES=21,TABMUL_CHAIN32=1,MODM31=2,MULTI_Q=1,L1CUDA=3,GRAPHS=0"
declare -A ARMS=(
  [C]="$BASE"
  [R1]="$BASE,REGCF3161=80"
  [R2]="$BASE,REGCF3161=80,REGTS61=80"
  [R3]="$BASE,REGMI61=48,REGMO61=48"
  [R4]="$BASE,REGCF3161=80,REGTS61=80,REGMI61=48,REGMO61=48"
  [R5]="INPLACE=1,LOADS=22042,STORES=21,TABMUL_CHAIN32=1,MODM31=2,MULTI_Q=1,L1CUDA=1,GRAPHS=0,REGCF3161=64,REGTS61=64"
)
ORDER=(C R1 C R2 C R3 C R4 C R5)
sudo nvidia-smi -pl 600 >/dev/null
for round in 1 2 3; do
  for a in "${ORDER[@]}"; do
    D=$(mktemp -d ".reg${a}r${round}.XXXXXX")
    ./build-cuda/prpll -dir "$D" -prp 136279841 -iters 100000 -fft "$FFT" \
      -use "${ARMS[$a]}" > "$D/run.out" 2>&1 || echo "FAIL $a r$round ($D)"
  done
  echo "== round $round done"
done
echo "== regsweep complete"
