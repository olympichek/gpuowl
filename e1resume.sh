#!/usr/bin/env bash
# Resume E1 after the k_lds OOB fix: finish p1a's missing tail, then p1b, p2.
set -euo pipefail
cd /home/ubuntu/gpuowl
CC=/tmp/claude-1001/-home-ubuntu-gpuowl/29ee1069-0039-4be4-8ffe-bebd253e2d4d/scratchpad/clockcost
OUT=e1-clockcost
SECS=60

run_all() {
  local f="$1"
  "$CC" header > "$f"
  for r in $("$CC" list); do
    "$CC" "$r" "$SECS" >> "$f"
  done
}

sudo nvidia-smi -pl 600 >/dev/null
sudo nvidia-smi -lgc 2100,2100 >/dev/null
for r in lds_lo lds_hi mix61 dram_lo dram_hi; do
  "$CC" "$r" "$SECS" >> "$OUT/p1a-lgc2100.csv"
done

sudo nvidia-smi -lgc 2400,2400 >/dev/null
run_all "$OUT/p1b-lgc2400.csv"

sudo nvidia-smi -rgc >/dev/null
sudo nvidia-smi -pl 300 >/dev/null
run_all "$OUT/p2-pl300.csv"

sudo nvidia-smi -pl 600 >/dev/null
sudo nvidia-smi -rgc >/dev/null
echo "== E1 sweep complete (resumed) =="
