#!/usr/bin/env bash
# E1 clock-cost table: three protocols over all clockcost runs.
#  P1a: locked 2100 MHz -> watts table (J/cycle = W/f)
#  P1b: locked 2400 MHz -> watts table (V-f scaling point)
#  P2:  300 W cap, unlocked -> sustained-MHz table (binding only for heavy mixes)
set -euo pipefail
cd /home/ubuntu/gpuowl
CC=/tmp/claude-1001/-home-ubuntu-gpuowl/29ee1069-0039-4be4-8ffe-bebd253e2d4d/scratchpad/clockcost
OUT=e1-clockcost
mkdir -p "$OUT"
SECS=60
nvidia-smi --query-gpu=driver_version,ecc.mode.current,power.limit,clocks.max.sm --format=csv > "$OUT/env.txt"

run_all() {
  local f="$1"
  "$CC" header > "$f"
  for r in $("$CC" list); do
    "$CC" "$r" "$SECS" >> "$f"
  done
}

sudo nvidia-smi -pl 600 >/dev/null
sudo nvidia-smi -lgc 2100,2100 >/dev/null
run_all "$OUT/p1a-lgc2100.csv"

sudo nvidia-smi -lgc 2400,2400 >/dev/null
run_all "$OUT/p1b-lgc2400.csv"

sudo nvidia-smi -rgc >/dev/null
sudo nvidia-smi -pl 300 >/dev/null
run_all "$OUT/p2-pl300.csv"

sudo nvidia-smi -pl 600 >/dev/null
sudo nvidia-smi -rgc >/dev/null
echo "== E1 sweep complete =="
