#!/usr/bin/env python3
"""E1 summary: merge the three protocol CSVs into the clock-cost table.

Usage: e1sum.py [e1-clockcost]
Outputs per run: W@2100, J/cycle@2100 (pJ), W@2400, V-f exponent proxy
(W2400/W2100 vs f-ratio 1.147), sustained MHz @300W cap (binding flag),
pJ/op@2100, temps; plus hi-lo toggle deltas per class.
"""
import csv, os, sys, collections

def load(path):
    rows = {}
    if not os.path.exists(path):
        return rows
    with open(path) as f:
        for r in csv.DictReader(f):
            rows[r["run"]] = {k: float(v) if k != "run" else v for k, v in r.items() if k != "run"}
    return rows

def main(root):
    p1a = load(os.path.join(root, "p1a-lgc2100.csv"))
    p1b = load(os.path.join(root, "p1b-lgc2400.csv"))
    p2 = load(os.path.join(root, "p2-pl300.csv"))
    runs = [r for r in p1a]
    print(f"{'run':16s} {'W@2100':>7s} {'pJ/cyc':>7s} {'W@2400':>7s} {'Wratio':>6s} {'MHz@300W':>8s} {'bind':>4s} {'pJ/op':>7s} {'T1a':>4s} {'T1b':>4s}")
    for r in runs:
        a = p1a.get(r, {})
        b = p1b.get(r, {})
        c = p2.get(r, {})
        w1 = a.get("watts_energy", 0)
        w2 = b.get("watts_energy", 0)
        f1 = a.get("median_mhz", 0)
        pjc = w1 / f1 * 1e3 / 188 if f1 else 0  # pJ per cycle per SM (chip W / f / SMs)
        wr = w2 / w1 if w1 else 0
        mhz300 = c.get("median_mhz", 0)
        wc = c.get("watts_energy", 0)
        bind = "Y" if wc >= 285 else "n"  # cap 300W: binding if drawing ~cap
        print(f"{r:16s} {w1:7.1f} {pjc:7.2f} {w2:7.1f} {wr:6.3f} {mhz300:8.0f} {bind:>4s} {a.get('pj_per_op',0):7.3f} {a.get('mean_temp',0):4.0f} {b.get('mean_temp',0):4.0f}")
    print("\nToggle deltas (hi - lo), W@2100 / W@2400 / MHz@300W:")
    seen = set()
    for r in runs:
        if r.endswith("_hi"):
            base = r[:-3] + "_lo"
            if base in p1a and base not in seen:
                seen.add(base)
                dw1 = p1a[r]["watts_energy"] - p1a[base]["watts_energy"]
                dw2 = (p1b.get(r, {}).get("watts_energy", 0) - p1b.get(base, {}).get("watts_energy", 0))
                dmh = (p2.get(r, {}).get("median_mhz", 0) - p2.get(base, {}).get("median_mhz", 0))
                print(f"  {r[:-3]:14s} {dw1:+7.1f} W  {dw2:+7.1f} W  {dmh:+6.0f} MHz")

if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "e1-clockcost")
